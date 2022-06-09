from pathlib import Path
from typing import IO, List
from multiprocessing import Queue, Process
import logging

import tqdm
import numpy as np
from meth5 import MetH5File, MethlyationValuesContainer
from meth5.sparse_matrix import SparseMethylationMatrixContainer

from pycoMeth.meth_seg.segments_csv_io import SegmentsWriterBED, SegmentsWriterBedGraph
from pycoMeth.meth_seg.segment import segment
from pycoMeth.common import pycoMethError, get_logger


def worker_segment(input_queue: Queue, output_queue: Queue, max_segments_per_window: int):
    import warnings
    
    warnings.filterwarnings("ignore")
    
    while True:
        job = input_queue.get()
        if job is None:
            break
        
        sparse_matrix, fraction = job
        llrs = sparse_matrix.met_matrix
        
        if sparse_matrix.shape[1] <= 1:
            # Too few CpG-sites. Nothing to segment.
            segmentation = np.zeros(sparse_matrix.shape[1])
            result_tuple = (
                llrs,
                segmentation,
                sparse_matrix.genomic_coord,
                sparse_matrix.genomic_coord_end,
                sparse_matrix.read_samples,
            )
        else:
            # Perform segmentation
            segmentation = segment(sparse_matrix, max_segments_per_window)
            result_tuple = (
                llrs,
                segmentation,
                sparse_matrix.genomic_coord,
                sparse_matrix.genomic_coord_end,
                sparse_matrix.read_samples,
            )
        output_queue.put((result_tuple, fraction))


def worker_output(
    output_queue: Queue,
    out_tsv_file: IO,
    out_bedgraph_filebase: str,
    chromosome: str,
    read_groups_keys: str,
    print_diff_met: bool,
    quiet: bool,
):
    writers = [SegmentsWriterBED(out_tsv_file, chromosome)]
    if out_bedgraph_filebase is not None:
        writers.append(SegmentsWriterBedGraph(out_bedgraph_filebase, chromosome))
    with tqdm.tqdm(total=100) as pbar:
        while True:
            res = output_queue.get()
            if res is None:
                break
            seg_result, fraction = res
            llrs, segments, genomic_starts, genomic_ends, samples = seg_result
            
            for writer in writers:
                writer.write_segments_llr(
                    llrs, segments, genomic_starts, genomic_ends, samples, compute_diffmet=print_diff_met
                )
            pbar.update(fraction)
        pbar.n = 100
        pbar.refresh()


def load_met_matrix(
    filename: str,
    values_container: MethlyationValuesContainer,
    read_groups_keys: List[str],
    read_groups_to_include: List[str],
) -> SparseMethylationMatrixContainer:
    met_matrix: SparseMethylationMatrixContainer = values_container.to_sparse_methylation_matrix(
        read_read_names=False, read_groups_key=read_groups_keys
    )
    
    if read_groups_keys is None:
        """Read groups are read names (read-level mode)"""
        met_matrix.read_samples = np.array([f"{filename}" for _ in met_matrix.read_names])
    else:
        if read_groups_to_include is not None:
            """Filter out only allowed read groups"""
            idx = np.array([r in read_groups_to_include for r in met_matrix.read_samples])
            met_matrix = met_matrix.get_submatrix_from_read_mask(idx)
        """Read groups are read from h5 file"""
        met_matrix.read_samples = np.array([f"{filename}_{sn}" for sn in met_matrix.read_samples])
    return met_matrix


def worker_reader(
    m5files: List[Path],
    chunk_size: int,
    chromosome: str,
    window_size: int,
    input_queue: Queue,
    chunks: List[int],
    progress_per_chunk: float,
    read_groups_keys: List[str],
    read_groups_to_include: List[str],
):
    firstfile = m5files[0]
    with MetH5File(firstfile, "r", chunk_size=chunk_size) as m5:
        chrom_container = m5[chromosome]
        
        for chunk in chunks:
            values_container = chrom_container.get_chunk(chunk)
            chunk_start, chunk_end = values_container.start, values_container.end
            met_matrix = load_met_matrix(firstfile.name, values_container, read_groups_keys, read_groups_to_include)
            
            for other_m5file in m5files[1:]:
                with MetH5File(other_m5file, "r", chunk_size=chunk_size) as other_m5:
                    other_values_container = other_m5[chromosome].get_values_in_range(chunk_start, chunk_end)
                    other_met_matrix = load_met_matrix(
                        other_m5file.name, other_values_container, read_groups_keys, read_groups_to_include
                    )
                    if met_matrix.met_matrix.shape[0] == 0:
                        # First file had no data in the requested samples
                        met_matrix = other_met_matrix
                    elif other_met_matrix.met_matrix.shape[0] > 0:
                        met_matrix = met_matrix.merge(other_met_matrix, sample_names_mode="keep")
            
            if read_groups_keys is None and len(m5files) == 1:
                met_matrix.read_samples = met_matrix.read_names
            total_sites = len(met_matrix.genomic_coord)
            num_windows = (total_sites // window_size) + 1
            progress_per_window = progress_per_chunk / num_windows
            for window_start in range(0, total_sites + 1, window_size):
                window_end = window_start + window_size
                logging.debug(f"Submitting window {window_start}-{window_end}")
                sub_matrix = met_matrix.get_submatrix(window_start, window_end)
                input_queue.put((sub_matrix, progress_per_window))


def validate_chromosome_selection(m5file: Path, chromosome: str, chunk_size: int):
    with MetH5File(m5file, "r", chunk_size=chunk_size) as m5:
        if chromosome not in m5.get_chromosomes():
            raise ValueError(f"Chromosome {chromosome} not found in m5 file.")


def validate_chunk_selection(m5file: Path, chromosome: str, chunk_size: int, chunks: List[int]):
    with MetH5File(m5file, "r", chunk_size=chunk_size) as m5:
        num_chunks = m5[chromosome].get_number_of_chunks()
        if max(chunks) >= m5[chromosome].get_number_of_chunks():
            raise ValueError(f"Chunk {max(chunks)} not in chromosome. Must be in range {0}-{num_chunks - 1}")


def Meth_Seg(
    h5_file_list: [Path],
    output_tsv_fn: str,
    chromosome: str,
    chunk_size: int = int(5e4),
    chunks: [int] = None,
    workers: int = 1,
    reader_workers: int = 1,
    progress: bool = False,
    window_size: int = 300,
    max_segments_per_window: int = 10,
    read_groups_keys: [str] = None,
    read_groups_to_include: [str] = None,
    print_diff_met: bool = False,
    output_bedgraph_fn: str = None,
    verbose: bool = False,
    quiet: bool = False,
    **kwargs,
):
    """
    Methylation segmentation method implemented as a bayesian changepoint detection algorithm
    * h5_file_list
        A list of MetH5 files containing methylation llr
    * chromosome
        The chromosome to segment
    * chunk_size
        Number of llrs per chunk - for best performance, should be a multiple of the  chunksize used in creating of the h5 files
        Default is the same as the default for creating meth5 files.
    * chunks
        List of chunk IDs or None if all chunks of the chromsome are to be segmented
    * workers
        Number of worker processes
    * reader_workers
        Number of reader worker processes
    * progress
        True if  progress bar is desired
    * output_tsv_fn
        Output TSV file
    * window_size
        Window size for segmentation in number of CpG calling sites. Default: 300.
        Increasing this increases memory requirement
    * max_segments_per_window
        Maximum number of segments per window. Should probably be somewhere between 8 and 20.
        The larger the number, the more expensive the computation.
    * read_groups_keys
        If read groups should be considered (e.g. haplotype) pass the read group key. You can provide more than one.
    * print_diff_met
        Whether output TSV file should contain methylation rate difference between samples
    * output_bedgraph_fn
        Base name for bedgraphs to be written. One bedgraph per sample/read_group will be created.
    """
    log = get_logger(name="pycoMeth_Meth_Seg", verbose=verbose, quiet=quiet)
    log.debug("Checking options and input files")
    if read_groups_to_include is not None and read_groups_keys is None:
        raise pycoMethError("read_groups_to_include defined, but missing read_groups_keys parameter")
    
    input_queue = Queue(maxsize=workers * 5)
    output_queue = Queue(maxsize=workers * 100)
    
    for m5file in h5_file_list:
        validate_chromosome_selection(m5file, chromosome, chunk_size)
    
    firstm5 = h5_file_list[0]
    if chunks is None:
        # No chunks have been provided, take all
        with MetH5File(firstm5, mode="r", chunk_size=chunk_size) as f:
            chunks = list(range(f[chromosome].get_number_of_chunks()))
    else:
        # flatten chunk list, since we allow a list of chunks or a list of chunk ranges
        # (which are converted to lists in parsing)
        chunks = [chunk for subchunks in chunks for chunk in ([subchunks] if isinstance(subchunks, int) else subchunks)]
    
    validate_chunk_selection(firstm5, chromosome, chunk_size, chunks)
    
    # sort and make unique
    chunks = sorted(list(set(chunks)))
    progress_per_chunk = 100 / len(chunks)
    
    segmentation_processes = [Process(target=worker_segment, args=(input_queue, output_queue, max_segments_per_window))]
    for p in segmentation_processes:
        p.start()
    
    reader_workers = min(reader_workers, len(chunks))
    chunk_per_process = np.array_split(chunks, reader_workers)
    reader_processes = [
        Process(
            target=worker_reader,
            args=(
                h5_file_list,
                chunk_size,
                chromosome,
                window_size,
                input_queue,
                p_chunks,
                progress_per_chunk,
                read_groups_keys,
                read_groups_to_include,
            ),
        )
        for p_chunks in chunk_per_process
    ]
    for p in reader_processes:
        p.start()
    
    output_process = Process(
        target=worker_output,
        args=(
            output_queue,
            output_tsv_fn,
            output_bedgraph_fn,
            chromosome,
            read_groups_keys,
            print_diff_met,
            not progress,
        ),
    )
    output_process.start()
    
    for p in reader_processes:
        p.join()
    
    # Deal poison pills to segmentation workers
    for p in segmentation_processes:
        input_queue.put(None)
    
    for p in segmentation_processes:
        p.join()
    
    # Deal poison pill to writer worker
    output_queue.put(None)
    output_process.join()
