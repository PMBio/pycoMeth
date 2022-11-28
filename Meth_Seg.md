# Meth_Seg

This sub-command implements a Bayesian changepoint detection algorithm implemented as an HMM. See our pre-print for more details on the method:
https://www.biorxiv.org/content/10.1101/2022.02.16.480699v2

Unlike the CGI-finder, this method attempts to find segments based on methylation information alone. 

Note that this method is computationally more expensive, and it is highly recommended **to parallelize on a compute cluster**.

## Usage
    pycometh Meth_Seg [-h] -i [H5_FILE_LIST [H5_FILE_LIST ...]] -c CHROMOSOME [-n [CHUNKS [CHUNKS ...]]] -t
                      OUTPUT_TSV_FN [-b OUTPUT_BEDGRAPH_FN] [-r [READ_GROUPS_KEYS [READ_GROUPS_KEYS ...]]]
                      [-s [READ_GROUPS_TO_INCLUDE [READ_GROUPS_TO_INCLUDE ...]]] [--chunk_size CHUNK_SIZE]
                      [-m MAX_SEGMENTS_PER_WINDOW] [-p WORKERS] [--reader_workers READER_WORKERS] 
                      [-w WINDOW_SIZE] [--print_diff_met]

## Arguments
 * `-i H5_FILE_LIST`: Input MetH5 files containing methylation calls
 * `-c CHROMOSOME`: Which chromosome to segment
 * `-n CHUNKS`: Which MetH5 chunks to segment, provided as list of integers (space separated). Optional. If not provided, all chunks of the chromosome will be segmented. 
 * `-t OUTPUT_TSV_FN`: Tab delimited output file
 * `-r READ_GROUPS_KEYS`: One or more read-group keys which store read-groups in the MetH5 file (space separated). Optional.
 * `-s READ_GROUPS_TO_INCLUDE`: List of read-groups which should be included. Optional. If not provided, all read-groups will be included.
 * `--chunk_size CHUNK_SIZE`: Number of methylation calls per chunk. Can deviate from the chunk size used in MetH5 creation, but it's highly recommended to at least use a multiple of that.
 * `-m MAX_SEGMENTS_PER_WINDOW`: The maximum number of segments in a window. This parameter controls segmentation granularity. Optional. Default: 10. 
 * `-p WORKERS`: Number of processes to spawn for segmentation. Optional. Default: 1
 * `--reader_workers READER_WORKERS`: Number of processes to spawn for input pre-loading. Optional. Default: 1
 * `-w WINDOW_SIZE`: Number of CpG sites per segmentation task. Optional. Default: 300
 * `--print_diff_met`: Make output TSV file include the methylation rate difference between read-groups

## Output format

The TSV output format contains one line per segments, with the following columns: 
* chromosome
* start position
* end position
* number of unique CpG call locations in the segment
* Optional, if `--print_diff_met` is provided: the methylation rate differences between read groups 

The BED output format is the minimal standard genomic [BED3](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format listing the coordinates of segments.


## Usage tips

### Read-groups

`pycometh Meth_Seg` is designed for multi-read-group (i.e. multi-sample or haplotyped or both) consensus segmentations.
That means if you have methylation calls from multiple biological samples, or with haplotype assignment, you will get one single
methylome segmentation that should fit all three groups.

There are two ways to define read-groups. One is via separate MetH5 files. Each input file will be considered as its own read-group.
The second one is via the read-group annotation in the MetH5 file. If your MetH5 file contains read-groups, you can define
under which read group key they are filed using the `-r` argument, and you can further filter which read-groups to include with the `-s`
argument. 

For example, let's say you have 2 biological samples stored in `case.m5` and `control.m5`. Let's also say both samples have been
haplotyped and haplotype assignment is stored in the MetH5 file under the read-group key `haplotype`. The read-groups are 
called `H1`, `H2` and `unphased`. Let's say you don't to include unphased reads in your segmentation. In that case, you would run:

    $ pycometh Meth_Seg -i case.m5 control.m5 -r haplotype -s H1 H2 [...other arguments...]

### Segmentation granularity

The arguments `-w` and `-m` control how granular the segmentation is supposed to be. The default values have shown to be a reasonable
granularity (about 10 segments max per 300 CpG-sites). If your device has lots of RAM you can increase the window size in order to
reduce granularity, or increase both values in order to reduce the number of windows while keeping the same granularity.

### Parallel processing

One of the advantages of the MetH5 format is that chunked data storage allows for good load balancing in parallel processing.
It is therefore recommended to run the segmentation on a compute cluster, using the `-n` argument, giving each compute job the 
same amount of chunks to compute. 

You can use the meth5 CLI to find out how many chunks you store in your MetH5 file. For example:

    $ meth5 --chunk_size 2000000 list_chunks -i Primary_cpg_ori.h5
    ------ Primary_cpg_ori.h5 -------
    | Chromosome | Number of chunks |
    | chr1       | 32               |
    | chr10      | 18               |
    | chr11      | 16               |
    | chr12      | 16               |
    | chr13      | 11               |
    <truncated>

Then you can, for example, parallelize by running 5 chunks per cluster job (writing only the relevant arguments)

    $ bsub pycometh Meth_Seg --chunk_size 200000 -c chr1 -n 1 2 3 4 5 [...other arguments...] # job 1
    $ bsub pycometh Meth_Seg --chunk_size 200000 -c chr1 -n 6 7 8 9 10 [...other arguments...] # job 2
    # etc...

and eventually merge the output.

Also, make use of the `-p` parameter to define the number of worker processes, but note that workload is split based on chunks,
so you don't gain anything from running 8 worker processes per cluster job if you only assign 5 chunks to one job. 

**The number of chunks per job is optimally a multiple of the number of worker processes** 