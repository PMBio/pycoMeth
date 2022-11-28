# Meth_Comp

**Generate a fully responsive interactive HTML report for significant differentially methylated intervals found with `pycoMeth Meth_Comp`.**



## Usage
    pycometh Comp_Report [-h] -i [H5_FILE_LIST [H5_FILE_LIST ...]] -f REF_FASTA_FN [-r READ_GROUPS_KEY] -c
                         METHCOMP_FN -g GFF3_FN [-o OUTDIR] [-s [SAMPLE_ID_LIST [SAMPLE_ID_LIST ...]]] 
                         [-d MAX_TSS_DISTANCE] [--pvalue_threshold PVALUE_THRESHOLD] [-n N_TOP]
                         [--min_diff_bs MIN_DIFF_BS] [--n_len_bin N_LEN_BIN] [--export_static_plots]
                         [--report_non_significant] [-v] [-q] [-p] [--min_diff_llr MIN_DIFF_LLR]

## Arguments
 * `-i H5_FILE_LIST`: Input MetH5 files containing methylation calls
 * `-f REF_FASTA_FN`: Filename to genome reference in FASTA format
 * `-r READ_GROUPS_KEYS`: One or more read-group keys which store read-groups in the MetH5 file (space separated). Optional.
 * `-c METHCOMP_FN`: Output TSV file from `pycometh Meth_Comp` 
 * `-g GFF3_FN`: Reference annotation as Ensembl GFF3 file
 * `-o OUTDIR`: Output directory. Optional. Default: `./`
 * `-s SAMPLE_ID_LIST`: Sample names as provided when running `pycometh Meth_Comp` 
 * `-d MAX_TSS_DISTANCE`: Maximum distance of an interval from a gene's TSS to be reported. Optional. Default: 100000
 * `--pvalue_threshold PVALUE_THRESHOLD`: Significance threshold. Optional. Default: 0.05
 * `-n N_TOP`: Number of intervals for which a detailed report should be produced. Optional. Default: 100
 * `--min_diff_bs MIN_DIFF_BS`: Minimum difference in methylation rate (minimum effect size) between samples to be reported. Optional. Default: 0.5
 * `--n_len_bin N_LEN_BIN`: Number of genomic intervals for the longest chromosome of the ideogram figure. Optional. Default: 500
 * `--min_diff_llr MIN_DIFF_LLR`: Minimal llr boundary for negative and positive median llr. 1 is recommanded for vizualization purposes.
 * `--export_static_plots`: Plot static images in addition to dynamic HTML reports. 
 * `-v`: verbose mode
 * `-q`: quiet mode
 * `-p`: display progress bar

## Input files

### pycoMeth CpG_Aggregate or Interval_Aggregate output file

A **TSV** output files generated with `pycoMeth Meth_comp`.

### Reference GFF3 file

An **Ensembl GFF3** file containing genomic annotations to extract transcripts with TSS close to the top most significant CpG Intervals.

* Main Ensembl: https://www.ensembl.org/index.html
* Ensembl Bacteria: http://bacteria.ensembl.org/index.html
* Ensembl Fungi: http://fungi.ensembl.org/index.html
* Ensembl Plants: http://plants.ensembl.org/index.html
* Ensembl Protists: http://protists.ensembl.org/index.html
* Ensembl Metazoa: http://metazoa.ensembl.org/index.html

### Reference FASTA file

FASTA reference file used for read alignment and Nanopolish. This file is required and used to sort the chromosomes in the ideogram plot.

## Output files

### Summary HTML report (pycoMeth_summary_report.html)

Entry point to the report containing information of methylation status at genomic level.

[Example Summary report](examples/pycometh_report1.html)

The report contains the following items:

* Overall summary
  Simple table containing counts of intervals and CpGs significantly differentially methylated.

* Methylation category counts by sample
  Interactive plotly stacked bar plot showing the number of methylated, unmethylated and ambiguous intervals for each samples.

* Methylation log-likelihood ratio by CpG interval
  Interactive plotly heatmap of the median llr for all significant intervals.
  Sites are ordered by genomic coordinates and samples are clustered by hierarchical clustering.

* Distribution of CpG methylation log-likelihood ratio by sample
  Interactive plotly ridgeplot of median llr distribution for all sample at interval resolution.
  Samples are ordered by descending overall median llr.

* Ideogram of the distribution of significant CpG intervals by chromosomic intervals

  Interactive chromosome ideogram plot showing the number of significant sites binned in large intervals.

* Top differentially methylated intervals
  Table containing details of the n top differentially methylated candidates ranked by adjusted pvalue.
  The table contains links to individual intervals HTML reports.

### Top intervals HTML reports

Individual reports for the top differentially methylated interval candidates.
The right side navigation bar allows to explore all the other intervals or go back to the summary report.

[Example Interval report](examples/pycometh_report_details.html)

The report contains the following items:

* Interval details
Simple table containing details about the interval

* Methylation log-likelihood ratio by CpG position
Interactive plotly heatmap of the median llr for all CpG positions in the interval.
Sites are ordered by genomic coordinates and samples are clustered by hierarchical clustering.

* Distribution of CpG methylation log-likelihood ratio by sample
Interactive plotly ridgeplot of median llr distribution for all sample at CpG resolution.
Samples are ordered by descending overall median llr.

* Closest transcripts TSS
Table containing the list of all transcripts having their TSS within 100kb of the interval edges.
Transcripts are ranked by TSS distance to the interval
