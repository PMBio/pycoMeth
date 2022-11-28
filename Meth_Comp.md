# Meth_Comp

Testing of segments for differential methylation between two or more samples.

Typically, you would a list of segments to be tested, either from an annotation file, CpG-Islands from the `CGI_Finder` 
subcommand or segments defined throuigh the `Meth_Seg`subcommand. Alternatively, you can also perform a windowed test by
providing an interval size.

## Usage
    pycometh Meth_Comp [-h] -i [H5_FILE_LIST [H5_FILE_LIST ...]] -f REF_FASTA_FN [-r READ_GROUPS_KEY]
                       [-a INTERVAL_BED_FN] [-b OUTPUT_BED_FN] [-t OUTPUT_TSV_FN] [-n INTERVAL_SIZE]
                       [-c MIN_NUM_READS_PER_INTERVAL] [-m MAX_MISSING] [-w WORKER_PROCESSES] [-l MIN_ABS_LLR]
                       [-s [SAMPLE_ID_LIST [SAMPLE_ID_LIST ...]]] [--pvalue_adj_method PVALUE_ADJ_METHOD]
                       [--pvalue_threshold PVALUE_THRESHOLD] [--only_tested_sites] [--hypothesis HYPOTHESIS]
                       [--do_independent_hypothesis_weighting] [-v] [-q] [-p]

## Arguments
 * `-i H5_FILE_LIST`: Input MetH5 files containing methylation calls
 * `-f REF_FASTA_FN`: Filename to genome reference in FASTA format
 * `-r READ_GROUPS_KEYS`: One or more read-group keys which store read-groups in the MetH5 file (space separated). Optional.
 * `-a INTERVAL_BED_FN`: Segments to be tested in tab-separated format or BED format. Must contain the 3 columns chromosome, 
start, and end choordinate. Optional. If no file is provided, a windowed test is performed instead.
 * `-b OUTPUT_BED_FN`: BED delimited output file. Optional.
 * `-t OUTPUT_TSV_FN`: Tab delimited output file. Optional.
 * `-n INTERVAL_SIZE`: The interval/window size in bases if a windowed test is to be performed. Optional. If not provided, `-a` must be provided instead. Default: 1000
 * `-c MIN_NUM_READS_PER_INTERVAL`: Minimum number of reads per interval for a sample to be considered in the test. Optional. Default: 10
 * `-m MAX_MISSING`: The maximum number of samples allowed missing per interval. If greater than 0, samples which do not 
have sufficient coverage in the interval will be excluded and the test may proceed between the remaining samples. Optional. Default: 0
 * `-w WORKER_PROCESSES`: Number of processes to spawn for testing. Optional. Default: 4
 * `-l MIN_ABS_LLR`: Mininum absolute log-likelihood ratio of methylation. Methylation calls with calls of higher uncertainty are ignored. Optional. Default: 2
 * `-s SAMPLE_ID_LIST`: A unique sample name for each MetH5 file. Optional. Will be numeric ids if none provided.
 * `--pvalue_adj_method PVALUE_ADJ_METHOD`: Method for adjustment of p-values for multiple testing. Optional. 
Default: `fdr_bh`. Choose from:
   - `bonferroni` : one-step correction
   - `sidak` : one-step correction
   - `holm-sidak` : step down method using Sidak adjustments
   - `holm` : step-down method using Bonferroni adjustments
   - `simes-hochberg` : step-up method  (independent)
   - `hommel` : closed method based on Simes tests (non-negative)
   - `fdr_bh` : Benjamini/Hochberg  (non-negative)
   - `fdr_by` : Benjamini/Yekutieli (negative)
   - `fdr_tsbh` : two stage fdr correction (non-negative)
   - `fdr_tsbky` : two stage fdr correction (non-negative)
 * `--pvalue_threshold PVALUE_THRESHOLD`: Significance threshold. Optional. Default: 0.05
 * `--only_tested_sites`: Exclude intervals which were not tested (e.g. due to insufficient coverage) from the output
 * `--hypothesis`: The testing hypothesis. Choose from:
    - `llr_diff`: Will perform a non-parametric test for shift in methylation log-likelihood ratio
    - `bs_diff`: Will compute a methylation rate per read in the interval and then perform a non-parametric test for a shift in read-methylation
    - `count_dependency`: Will binarize methylation calls and compute a contingency table for a fisher exact (2 samples) or chi-squared test (more than 2 samples)
 * `--do_independent_hypothesis_weighting`: Adjusts raw p-values by wheighing them using methylation rate variance in the segment as a weight. 
 * `-v`: verbose mode
 * `-q`: quiet mode
 * `-p`: display a progress bar

## Output format

The TSV output format contains one line per tested segment, with the following columns:

* chromosome
* start position
* end position
* number of samples with sufficient coverage
* raw p-value
* adjusted p-value
* number of unique cpg positions
* sample ids
* medium log-likelihood ratio for each sample, order matches the sample id column
* methylation rate difference between samples, order matches the sample id column with first entry being sample1 - sample2, the next being sample2 - sample3, etc
* average covereage per sample
* comment

