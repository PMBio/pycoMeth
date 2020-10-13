# -*- coding: utf-8 -*-

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Standard library imports
import itertools
from typing import IO, List, Generator
import fileinput

# Third party imports
from tqdm import tqdm
import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu
from statsmodels.stats.multitest import multipletests
from multiprocessing import Pool
from meth5.meth5_wrapper import MetH5File

# Local imports
from pycoMeth.common import *
from pycoMeth.FileParser import FileParser
from pycoMeth.CoordGen import CoordGen, Coord

# ~~~~~~~~~~~~~~~~~ Multiprocessing Worker methods ~~~~~~~~~~~~~~~~~~~~~~~~~~#


class MethCompWorker:
    def __init__(
        self,
        h5_read_groups_key,
        sample_id_list,
        h5_file_list,
        min_diff_llr,
        min_samples,
        pvalue_method,
        min_num_reads_per_interval,
    ):
        self.h5_read_groups_key = h5_read_groups_key
        self.min_diff_llr = min_diff_llr
        self.min_samples = min_samples
        self.pvalue_method = pvalue_method
        self.min_pval = np.nextafter(float(0), float(1))
        self.min_num_reads_per_interval = min_num_reads_per_interval
        self.sample_hf_files = {}

        if h5_read_groups_key is None:
            for sample_id, h5_file in zip(sample_id_list, h5_file_list):
                hf = MetH5File(h5_file, "r")
                self.sample_hf_files[sample_id] = hf
        else:
            hf = MetH5File(h5_file_list[0], "r")
            for sample_id in sample_id_list:
                self.sample_hf_files[sample_id] = hf

    def __del__(self):
        for hf in self.sample_hf_files.values():
            try:
                hf.close()
            except:
                pass

    def compute_pvalue(self, interval, sample_llrs, sample_pos, sample_reads):

        label_list = list([k for k in sample_llrs.keys() if len(sample_llrs[k]) > 0])
        raw_llr_list = [sample_llrs[k] for k in label_list]
        raw_pos_list = [sample_pos[k] for k in label_list]

        n_samples = sum(1 for llrs in raw_llr_list if len(llrs) > 0)
        
        n_reads = [len(sample_reads[k]) for k in label_list]

        # Collect median llr
        med_llr_list = [np.median(llrs) for llrs in raw_llr_list]

        # Evaluate median llr value
        neg_med = sum(med <= -self.min_diff_llr for med in med_llr_list)
        pos_med = sum(med >= self.min_diff_llr for med in med_llr_list)
        ambiguous_med = sum(
            -self.min_diff_llr <= med <= self.min_diff_llr for med in med_llr_list
        )

        counters_to_increase = []

        if n_samples < self.min_samples:
            # Not enough samples
            comment = "Insufficient samples"
            pvalue = np.nan
        elif np.min(n_reads) < self.min_num_reads_per_interval:
            # Not enough reads in one of the samples
            comment = "Insufficient coverage"
            pvalue = np.nan
        # Sufficient samples and effect size
        elif not neg_med or not pos_med:
            comment = "Insufficient effect size"
            pvalue = np.nan

        # Sufficient samples and effect size
        else:
            comment = "Valid"
            # Run stat test
            if self.pvalue_method == "KW":
                statistics, pvalue = kruskal(*raw_llr_list)
            elif self.pvalue_method == "MW":
                statistics, pvalue = mannwhitneyu(raw_llr_list[0], raw_llr_list[1])

            # Fix and categorize p-values
            if pvalue is np.nan or pvalue is None or pvalue > 1 or pvalue < 0:
                pvalue = 1.0
                counters_to_increase.append("Sites with invalid pvalue")

            # Correct very low pvalues to minimal float size
            elif pvalue == 0:
                pvalue = self.min_pval

        # Update counters result table
        counters_to_increase.append(comment)

        res = OrderedDict()
        res["chromosome"] = interval.chr_name
        res["start"] = interval.start
        res["end"] = interval.end
        res["pvalue"] = pvalue
        res["adj_pvalue"] = np.nan
        res["n_samples"] = n_samples
        res["neg_med"] = neg_med
        res["pos_med"] = pos_med
        res["ambiguous_med"] = ambiguous_med
        res["label_list"] = list_to_str(label_list)
        res["med_llr_list"] = list_to_str(med_llr_list)
        res["raw_llr_list"] = list_to_str(raw_llr_list)
        res["comment"] = comment
        res["raw_pos_list"] = list_to_str(raw_pos_list)
        res["unique_cpg_pos"] = len(set(itertools.chain.from_iterable(raw_pos_list)))

        return res, counters_to_increase

    def __call__(self, interval):
        sample_llrs = {}
        sample_pos = {}
        sample_reads = {}

        for sample_id, hf in self.sample_hf_files.items():
            chrom_container = hf[interval.chr_name]
            interval_container = chrom_container.get_values_in_range(
                interval.start, interval.end
            )
            
            if interval_container is None:
                continue
                
            llrs = interval_container.get_llrs()[:]
            pos = interval_container.get_ranges()[:, 0]
            read_names = interval_container.get_read_names()[:]

            if self.h5_read_groups_key is not None:
                read_samples = interval_container.get_read_groups(
                    self.h5_read_groups_key
                )
                
                mask = read_samples == sample_id
                llrs = llrs[mask]
                pos = pos[mask]
                read_names = read_names[mask]
            sample_llrs[sample_id] = llrs.tolist()
            sample_pos[sample_id] = pos.tolist()
            sample_reads[sample_id] = read_names

        return self.compute_pvalue(interval, sample_llrs, sample_pos, sample_reads)


def initializer(**args):
    """Initializes a worker object at the beginning when the
    multiprocessing pool is created and puts it in the global
    namespace."""
    global worker
    worker = MethCompWorker(**args)


def worker_function(*args):
    """Calls the work function of the worker object in the global
    namespace."""
    return worker(*args)


# ~~~~~~~~~~~~~~~~~~~~~~~~CpG_Comp MAIN CLASS~~~~~~~~~~~~~~~~~~~~~~~~#


def Meth_Comp(
    h5_file_list: [str],
    ref_fasta_fn: str,
    read_group_file: str = None,
    interval_bed_fn: str = None,
    output_bed_fn: str = None,
    output_tsv_fn: str = None,
    interval_size: int = 1000,
    min_num_reads_per_interval: int = 10,
    max_missing: int = 0,
    min_diff_llr: float = 2,
    sample_id_list: [str] = None,
    pvalue_adj_method: str = "fdr_bh",
    pvalue_threshold: float = 0.01,
    only_tested_sites: bool = False,
    verbose: bool = False,
    quiet: bool = False,
    progress: bool = False,
    worker_processes: int = 4,
    
    **kwargs,
):
    """Compare methylation values for each CpG positions or intervals
    between n samples and perform a statistical test to evaluate if the
    positions are significantly different. For 2 samples a Mann_Withney
    test is performed otherwise multiples samples are compared with a
    Kruskal Wallis test. pValues are adjusted for multiple tests using
    the Benjamini & Hochberg procedure for controlling the false
    discovery rate.

    * h5_file_list
        A list of MetH5 files containing methylation llr
    * read_group_file
        Tab-delimited file assigning reads to read groups (e.g. samples or haplotypes). (optional)
    * ref_fasta_fn
        Reference file used for alignment in Fasta format (ideally already indexed with samtools faidx)
    * interval_bed_fn
        SORTED bed file containing **non-overlapping** intervals to bin CpG data into (Optional) (can be gzipped)
    * interval_size
        Size of the sliding window in which to aggregate CpG sites data from if no BED file is provided
   * min_num_reads_per_interval
        Minimum number of reads per sample per interval. The entire interval will be discarded if one sample
        does not have sufficient coverage.
    * output_bed_fn
        Path to write a summary result file in BED format (At least 1 output file is required) (can be gzipped)
    * output_tsv_fn
        Path to write an more extensive result report in TSV format (At least 1 output file is required) (can be gzipped)
    * max_missing
        Max number of missing samples to perform the test
    * min_diff_llr
        Minimal llr boundary for negative and positive median llr.
        The test if only performed if at least one sample has a median llr above (methylated) and 1 sample has a median llr below (unmethylated)
    * sample_id_list
        list of sample ids to annotate results in tsv file
    * pvalue_adj_method
        Method to use for pValue multiple test adjustment
    * pvalue_threshold
        Alpha parameter (family-wise error rate) for pValue adjustment
    * only_tested_sites
        Do not include sites that were not tested because of insufficient samples or effect size in the report
    * worker_processes
        Number of processes to be launched
    """
    # Init method
    opt_summary_dict = opt_summary(local_opt=locals())
    log = get_logger(name="pycoMeth_CpG_Comp", verbose=verbose, quiet=quiet)

    log.warning("Checking options and input files")
    log_dict(opt_summary_dict, log.debug, "Options summary")

    # At least one output file is required, otherwise it doesn't make any sense
    log.debug("Checking required output")

    if not output_bed_fn and not output_tsv_fn:
        raise pycoMethError("At least 1 output file is requires (-t or -b)")

    # Automatically define tests and maximal missing samples depending on number of files to compare
    read_sample_assignment = None
    if read_group_file is not None:
        read_sample_assignment = read_readgroups_file(read_group_file)
        read_sample_assignment = read_sample_assignment.loc[
            read_sample_assignment["group_set"] != -1
        ]
        read_sample_assignment = read_sample_assignment.to_dict()["group"]
        # Number of read groups minus the "-1" which stands for "unphased"
        sample_id_list = sorted(
            list(set(read_sample_assignment.values()).difference({-1}))
        )
    else:

        if not sample_id_list:
            sample_id_list = list(range(len(h5_file_list)))

    all_samples = len(sample_id_list)

    min_samples = all_samples - max_missing

    coordgen = CoordGen(ref_fasta_fn, verbose, quiet)
    log_list(coordgen, log.debug, "Coordinate reference summary")

    if interval_bed_fn:
        log.debug("Bed annotation generator")

        def intervals_gen_fun():
            return bed_intervals_gen(coordgen=coordgen, interval_bed_fn=interval_bed_fn)

    else:
        log.debug("Sliding window generator")

        def intervals_gen_fun():
            return sliding_intervals_gen(coordgen=coordgen, interval_size=interval_size)

    # Go through intervals once (should be cheap) to count how many we are investigating
    num_intervals = sum(1 for _ in intervals_gen_fun())
    # Recreate the intervals generator
    intervals_gen = intervals_gen_fun()

    # 3 values = Kruskal Wallis test
    if all_samples >= 3:
        pvalue_method = "KW"
        log.debug("Multiple comparison mode (Kruskal_Wallis test)")
        if min_samples < 3:
            log.debug("Automatically raise number of minimal samples to 3")
            min_samples = 3
    # 2 values = Mann_Withney test
    elif all_samples == 2:
        pvalue_method = "MW"
        log.debug("Pairwise comparison mode (Mann_Withney test)")
        if min_samples:
            log.debug("No missing samples allowed for 2 samples comparison")
            min_samples = 2
    else:
        raise pycoMethError("Meth_Comp needs at least 2 input files")

    log.warning("Opening H5 files")
    try:
        # Define StatsResults to collect valid sites and perform stats
        stats_results = StatsResults(
            pvalue_adj_method=pvalue_adj_method,
            pvalue_threshold=pvalue_threshold,
            only_tested_sites=only_tested_sites,
        )

        h5_read_groups_key = "pycometh_rg" if  read_group_file is not None else None# TODO expose parameter
        # Ensure every h5file is readable and has an index
        try:
            for h5_file in h5_file_list:
                hf = MetH5File(h5_file, "a")
                hf.create_chunk_index(force_update=False)
                if read_sample_assignment is not None:
                    hf.annotate_read_groups(
                        h5_read_groups_key,
                        read_sample_assignment,
                        exists_ok=True,
                        overwrite=False,
                    )
                hf.close()
        except:
            raise pycoMethError(
                "Unable to read/write h5 files. Must be writable to create index!"
            )

        del read_sample_assignment

        log.info("Starting asynchronous file parsing")
        with tqdm(
            total=num_intervals,
            unit=" intervals",
            unit_scale=True,
            desc="\tProgress",
            disable=not progress,
        ) as pbar:

            log.info("Launching %d worker processes" % worker_processes)

            pool = Pool(
                worker_processes,
                initializer=lambda: initializer(
                    h5_read_groups_key=h5_read_groups_key,
                    sample_id_list=sample_id_list,
                    h5_file_list=h5_file_list,
                    min_diff_llr=min_diff_llr,
                    min_samples=min_samples,
                    pvalue_method=pvalue_method,
                    min_num_reads_per_interval=min_num_reads_per_interval,
                ),
            )

            # Continue reading lines from all files
            log.debug("Starting deep parsing")
            fp_done = 0

            # Init file writer
            with Comp_Writer(
                bed_fn=output_bed_fn, tsv_fn=output_tsv_fn, verbose=verbose,
            ) as writer:

                def callback(*args):
                    result_line = stats_results.callback(*(args[0]))
                    writer.write(result_line)
                    pbar.update(1)

                abort = False

                def error_callback(err):
                    log.critical("Error in worker thread ")
                    log.critical(str(err))
                    nonlocal abort
                    abort = True

                async_results = []
                # TODO perhaps perform this in batches (e.g. submit 10k intervals,
                #  wait for all to finish, then submit the next 10k, etc...)
                #  instead of submitting every interval into the queue and then waiting
                #  for all of them to finish.
                #  That would allow for a more reasonable timeout.
                for interval in intervals_gen:
                    if abort:
                        raise pycoMethError("Aborting due to error in worker thread")
                    ar = pool.apply_async(
                        worker_function,
                        args=[interval],
                        callback=callback,
                        error_callback=error_callback,
                    )
                    async_results.append(ar)

                for i, ar in enumerate(async_results):
                    if abort:
                        break
                    ar.wait(timeout=3 * 24 * 3600)

            # Exit condition
            if not stats_results.res_list:
                log.info("No valid p-Value could be computed")

            else:
                # Convert results to dataframe and correct pvalues for multiple tests
                log.info("Adjust pvalues")
                stats_results.multitest_adjust()

                rewriter = Comp_ReWriter(
                    [f for f in (output_bed_fn, output_tsv_fn) if f is not None]
                )
                rewriter.write_adjusted_pvalues(stats_results.res_list)

    except:
        writer.abort()
        raise
    finally:
        # Print counters
        log_dict(stats_results.counter, log.info, "Results summary")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~StatsResults HELPER CLASS~~~~~~~~~~~~~~~~~~~~~~~~~~~#


class StatsResults:
    def __init__(
        self,
        pvalue_adj_method="fdr_bh",
        pvalue_threshold=0.01,
        only_tested_sites=False,
    ):
        """"""
        # Save self variables
        self.pvalue_adj_method = pvalue_adj_method
        self.pvalue_threshold = pvalue_threshold
        self.only_tested_sites = only_tested_sites

        # Init self collections
        self.res_list = []
        self.counter = Counter()

        # Get minimal non-zero float value
        self.min_pval = np.nextafter(float(0), float(1))

    # ~~~~~~~~~~~~~~MAGIC AND PROPERTY METHODS~~~~~~~~~~~~~~#

    def __repr__(self):
        return dict_to_str(self.counter)

    def __len__(self):
        return len(self.res_list)

    def __iter__(self):
        for i in self.res_list:
            yield i

    # ~~~~~~~~~~~~~~PUBLIC METHODS~~~~~~~~~~~~~~#

    def callback(self, res, counters_to_increase):
        # filter out non tested site if required
        if self.only_tested_sites and res["pvalue"] is np.nan:
            return

        for c in counters_to_increase:
            self.counter[c] += 1

        reduced_res = {
            "pvalue": res["pvalue"],
            "adj_pvalue": res["adj_pvalue"],
            "comment": res["comment"],
        }
        self.res_list.append(reduced_res)

        return res

    def multitest_adjust(self):
        """"""
        # Collect non-nan pvalues
        pvalue_idx = []
        pvalue_list = []
        for i, res in enumerate(self.res_list):
            if not np.isnan(res["pvalue"]):
                pvalue_idx.append(i)
                pvalue_list.append(res["pvalue"])

        # Adjust values
        adj_pvalue_list = multipletests(
            pvals=pvalue_list,
            alpha=self.pvalue_threshold,
            method=self.pvalue_adj_method,
        )[1]

        # add adjusted values to appropriate category
        for i, adj_pvalue in zip(pvalue_idx, adj_pvalue_list):

            # Fix and categorize p-values
            if (
                adj_pvalue is np.nan
                or adj_pvalue is None
                or adj_pvalue > 1
                or adj_pvalue < 0
            ):
                adj_pvalue = 1.0
                comment = "Non-significant pvalue"

            elif adj_pvalue <= self.pvalue_threshold:
                # Correct very low pvalues to minimal float size
                if adj_pvalue == 0:
                    adj_pvalue = self.min_pval
                # update counter if pval is still significant after adjustment
                comment = "Significant pvalue"
            else:
                comment = "Non-significant pvalue"

            # update counters and update comment and adj p-value
            self.counter[comment] += 1
            self.res_list[i]["comment"] = comment
            self.res_list[i]["adj_pvalue"] = adj_pvalue


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~Comp_Writer HELPER CLASS~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Comp_Writer:
    """Extract data for valid sites and write to BED and/or TSV file."""

    def __init__(self, bed_fn=None, tsv_fn=None, verbose=True):
        """"""
        self.log = get_logger(name="Comp_Writer", verbose=verbose)
        self.bed_fn = bed_fn
        self.tsv_fn = tsv_fn

        # Init file pointers
        self.bed_fp = self._init_bed() if bed_fn else None
        self.tsv_fp = self._init_tsv() if tsv_fn else None

        # Color score table
        self.colors = OrderedDict()
        self.colors[10] = "10,7,35"
        self.colors[9] = "32,12,74"
        self.colors[8] = "60,9,101"
        self.colors[7] = "87,15,109"
        self.colors[6] = "112,25,110"
        self.colors[5] = "137,34,105"
        self.colors[4] = "163,43,97"
        self.colors[3] = "187,55,84"
        self.colors[2] = "209,70,67"
        self.colors[1] = "230,230,230"
        self.colors[0] = "230,230,230"
        self.aborted = False

    # ~~~~~~~~~~~~~~PUBLIC METHODS~~~~~~~~~~~~~~#

    def write(self, res):
        """"""
        if self.bed_fn:
            self._write_bed(res)
        if self.tsv_fn:
            self._write_tsv(res)

    def __enter__(self):
        self.log.debug("Opening Writer")
        return self

    def __exit__(self, exception_type, exception_val, trace):
        self.log.debug("Closing Writer")
        for fp in (self.bed_fp, self.tsv_fp):
            try:
                fp.close()
                if self.aborted:
                    # There was an error - delete the partial file
                    os.remove(fp.name)
            except:
                pass

    def abort(self):
        self.aborted = True

    # ~~~~~~~~~~~~~~PRIVATE METHODS~~~~~~~~~~~~~~#
    def _init_bed(self):
        """Open BED file and write file header."""
        self.log.debug("Initialise output bed file")
        mkbasedir(self.bed_fn, exist_ok=True)
        fp = (
            gzip.open(self.bed_fn, "wt")
            if self.bed_fn.endswith(".gz")
            else open(self.bed_fn, "w")
        )
        # Write header line
        fp.write("track name=meth_comp itemRgb=On\n")
        return fp

    def _write_bed(self, res):
        """Write line to BED file."""
        # Log transform pvalue and cast to int
        if np.isnan(res["adj_pvalue"]):
            score = 0
        else:
            score = int(-np.log10(res["adj_pvalue"]))
        # Define color for bed file
        color = self.colors.get(score, self.colors[10])
        # Write line
        res_line = [
            res["chromosome"],
            res["start"],
            res["end"],
            ".",
            score,
            ".",
            res["start"],
            res["end"],
            color,
        ]
        self.bed_fp.write(str_join(res_line, sep="\t", line_end="\n"))

    def _init_tsv(self):
        """Open TSV file and write file header."""
        self.log.debug("Initialise output tsv file")
        mkbasedir(self.tsv_fn, exist_ok=True)
        fp = (
            gzip.open(self.tsv_fn, "wt")
            if self.tsv_fn.endswith(".gz")
            else open(self.tsv_fn, "w")
        )
        # Write header line

        header = [
            "chromosome",
            "start",
            "end",
            "n_samples",
            "pvalue",
            "adj_pvalue",
            "neg_med",
            "pos_med",
            "ambiguous_med",
            "unique_cpg_pos",
            "labels",
            "med_llr_list",
            "raw_llr_list",
            "raw_pos_list",
            "comment",
        ]
        fp.write(str_join(header, sep="\t", line_end="\n"))
        return fp

    def _write_tsv(self, res):
        """Write line to TSV file."""

        res_line = [
            res["chromosome"],
            res["start"],
            res["end"],
            res["n_samples"],
            res["pvalue"],
            res["adj_pvalue"],
            res["neg_med"],
            res["pos_med"],
            res["ambiguous_med"],
            res["unique_cpg_pos"],
            res["label_list"],
            res["med_llr_list"],
            res["raw_llr_list"],
            res["raw_pos_list"],
            res["comment"],
        ]
        self.tsv_fp.write(str_join(res_line, sep="\t", line_end="\n"))


class Comp_ReWriter:
    """Reads delimited files with a header line, and rewrites the file
    in-place using python 3's fileinput module.

    This way we don't have to hold the entire file in memory or copy it
    around.
    """

    def __init__(self, filenames: List[str], separators: List[str] = None):
        """
        :param filenames: The list of filenames to modify
        :param separators: The separators for each file. Can be None (default), in
                           which case tab is assumed. If provided, it must be the same
                           length as filenames
        """

        self.filenames = filenames
        if separators is None:
            self.separators = ["\t"] * len(self.filenames)
        else:
            assert len(separators) == len(self.filenames)
            self.separators = separators

    def write_adjusted_pvalues(self, res_list):
        for sep, filename in zip(self.separators, self.filenames):
            if filename is not None:
                with fileinput.input(filename, inplace=True) as fi_fp:
                    for res, line in zip((None, *res_list), fi_fp):
                        line = line.strip()
                        if res is None:
                            header = {k: i for i, k in enumerate(line.split(sep))}
                            print(line)
                        else:
                            updated_line = line.split(sep)
                            updated_line[header["adj_pvalue"]] = str(res["adj_pvalue"])
                            updated_line[header["comment"]] = res["comment"]
                            print(sep.join(updated_line))


def read_readgroups_file(readgroups_file: IO):
    """Reads file that assigns read to read groups (such as haplotypes,
    samples, clusters, etc)

    :param readgroups_file: path to the tab-separated file
    :return: pandas dataframe with columns "read_name", "group" and "group_set"
    """
    # Loading
    try:
        read_groups = pd.read_csv(
            readgroups_file,
            sep="\t",
            header=0,
            dtype={"read_name": str, "group": int, "group_set": "category"},
        )
    except Exception as e:
        raise pycoMethError("Unable to read read groups file", e)

    # Validation
    if len(read_groups.columns) == 2:
        should_colnames = ["read_name", "group"]
    elif len(read_groups.columns) == 3:
        should_colnames = ["read_name", "group", "group_set"]
    else:
        raise pycoMethError(
            "Invalid number of columns in read groups file (should be 2 or 3)"
        )

    if not all([col in read_groups.columns for col in should_colnames]):
        raise pycoMethError(
            "Invalid column names in read groups file (should be %s)"
            % should_colnames.join(", ")
        )

    # Finished validation, now add group_set column if not present
    if "group_set" not in read_groups.columns:
        read_groups["group_set"] = 1

    read_groups = read_groups.set_index("read_name")

    return read_groups


def sliding_intervals_gen(coordgen, interval_size=1000) -> Generator[Coord, None, None]:
    """Generate sliding window coordinate intervals over the entire
    reference genome provided."""
    for chr_name, chr_len in coordgen.chr_name_len.items():
        for start in range(0, chr_len, interval_size):
            end = start + interval_size if start + interval_size <= chr_len else chr_len
            yield (coordgen(chr_name, start, end))


def bed_intervals_gen(coordgen, interval_bed_fn) -> Generator[Coord, None, None]:
    """Generate coordinate intervals corresponding to the provided bed
    file."""
    with FileParser(
        fn=interval_bed_fn,
        colnames=["chrom", "start", "end"],
        dtypes={"start": int, "end": int},
        force_col_len=False,
        comment="track",
        quiet=True,
    ) as bed:

        prev_ct = None
        for line in bed:
            ct = coordgen(line.chrom, line.start, line.end)
            if prev_ct and ct < prev_ct:
                raise ValueError(
                    "Unsorted coordinate found in bed file {} found after {}. Chromosomes have to be ordered as in fasta reference file".format(
                        ct, prev_ct
                    )
                )
            prev_ct = ct
            yield (ct)
