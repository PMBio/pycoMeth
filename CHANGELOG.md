# pycoMeth Changelog

### 2022/06/09 v2.2.1
* Meth_Seg now allows for filtering which read-groups to be used in the segmentation (parameter `--read_groups_to_include`)  
* Multiple bugfixes

### 2022/05/16 v2.2.0
* Implemented Fisher-Exact test for 2-sample testing and Chi-squared test for multi-sample testing (with Fisher exact as a post-hoc test)

### 2022/03/07 v2.1.1
* Implemented post-hoc testing and reporting for multi-sample tests
* Multiple bugfixes and performance improvements

### 2022/02/16 v2.0.0
* First release of PycoMeth 2.0 using the [MetH5](http://github.com/snajder-r/MetH5Format) as an input.
* Meth_Seg is now a firm component of PycoMeth, allowing for de novo methylaiton segmentation
* Additional testing modes: paired-CpG methylation rate testing as well as unpaired read-methylation rate testing
* CpG_Aggregate and Interval_Aggregate removed as these are no longer necessary since random access to methylation calls is possible in the MetH5 Format
* Meth_Comp can now be computed on multiple CPUs. 

### 2020/07/16 v-0.4.14

* Switch from orca to kaleido for static export
* Export summary of intervals in table file for all intervals and not just for top hits
* Add new plot for distance between CpG islands and closest tss

### 2020/07/09 v-0.4.8

* Speed improvement to Comp_report
* Add interactive API to Comp_report
* Add ability to export static plots using orca

### 2020/01/15 v-0.4.5

* Add tabular text reports to Comp_report
* Add (default) option to write out all the intervals in Meth_Comp with reasons why excluded or included in DM analysis
* Improve Comp_report summary table and include new fields from Meth_Comp
* Tidy output folder structure for reports generated by meth_report
* Add Chromosome ideogram plot to summary report
* A fasta reference is now required to run Comp_Report
 
### 2020/04/01 v-0.4.0

* General improvement of logging message output
* Implement fancy color logger
* Add position tracking within intervals
* Add Comp_Report module to generate HTML reports of significant candidates

### 09/10/2019 v-0.1.0

* Deep refactoring of Freq_meth_calculate to Aggregate
* Major doc update including doc folders reorganisation
* Implement autodoc from docstring
* Fix and test CLI
