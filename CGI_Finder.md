# CGI_Finder

Simple method to find putative CpG islands in DNA sequences by using a sliding window and merging overlapping windows satisfying the CpG island definition. Results can be saved in bed and tsv format.

In order for a CpG island to be called, it must fulfill three criteria:

* Sufficient lenght. Threshold set my `-w` argument, default 0.
* Sufficient GC content computed as (N_C + N_G) / N. Threshold set by `-c` argument, default 0.5.
* Sufficient normalized CpG ratio, computed as (N_CG * N) / (N_C * N_G). Threshold set by `-r` argument, default 0.6.

## Usage
    pycometh CGI_Finder [-h] -f REF_FASTA_FN [-b OUTPUT_BED_FN] [-t OUTPUT_TSV_FN] [-m MERGE_GAP] 
                        [-w MIN_WIN_LEN] [-c MIN_CG_FREQ] [-r MIN_OBS_CG_RATIO] [-v] [-q] [-p]

## Arguments
 * `-f REF_FASTA_FN`: Filename to genome reference in FASTA format
 * `-b OUTPUT_BED_FN`: Output file in BED format. Optional.
 * `-t OUTPUT_TSV_FN`: Output file in tab-delimited format. Optional.
 * `-m MERGE_GAP`: minimum distance between CpG islands. CpG islands closer than this will be merged to one. Optional. Default: 0
 * `-w MIN_WIN_LEN`: Minimum CpG island length. Optional. Default: 0
 * `-c MIN_CG_FREQ`: Minimum GC frequency. Optional. Default: 0.5
 * `-r MIN_OBS_CG_RATIO`: Minimum normalized CpG ratio. Optional. Default: 0.6
 * `-v`: verbose mode
 * `-q`: quiet mode
 * `-p`: display progress bar

## Output files

### Tabulated TSV file

This tabulated file contains the following fields for each CpG island found:

* chromosome / start / end : Genomic coordinates
* length: Length of the interval
* num_CpG: Number of CpGs found
* CG_freq: G+C nucleotide frequency
* obs_exp_freq: Observed versus expected CpG frequency

### BED file

Minimal standard genomic [BED3](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format listing the coordinates of putative CpG islands.
