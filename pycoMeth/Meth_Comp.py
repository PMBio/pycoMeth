# -*- coding: utf-8 -*-

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Standard library imports
import itertools
import random
from typing import IO, List, Dict, Any, Generator
from math import sqrt
import fileinput

# Third party imports
from tqdm import tqdm
import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu, wilcoxon, fisher_exact, chi2_contingency
from statsmodels.stats.multitest import multipletests
from multiprocessing import Pool
from meth5 import MetH5File

# Local imports
from pycoMeth.common import *
from pycoMeth.FileParser import FileParser
from pycoMeth.CoordGen import CoordGen, Coord
from pycoMeth.loader import MetH5Loader

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
        pvalue_threshold,
        hypothesis,
        do_independent_hypothesis_weighting,
    ):
        self.min_diff_llr = min_diff_llr
        self.min_samples = min_samples
        self.pvalue_method = pvalue_method
        self.min_pval = np.nextafter(float(0), float(1))
        self.min_num_reads_per_interval = min_num_reads_per_interval
        self.sample_hf_files: Dict[str, MetH5File] = {}
        self.llr_threshold = 2.0  # TODO expose parameter
        self.min_diff_bs = 0.25  # TODO expose parameter
        self.pvalue_threshold = pvalue_threshold
        self.hypothesis = hypothesis
        self.do_independent_hypothesis_weighting = do_independent_hypothesis_weighting
        
        self.loader = MetH5Loader(h5_read_groups_key, sample_id_list, h5_file_list)
    
    def compute_ihw_weight(self, test_values: List[List[float]]) -> float:
        flat_list = [v for vl in test_values for v in vl]
        mean = sum(flat_list) / len(flat_list)
        variance = sum((v - mean) ** 2 for v in flat_list) / (len(flat_list) - 1)
        std = sqrt(variance)
        return std
    
    def compute_site_betascores(
        self, raw_pos_list: List[List[Any]], raw_llr_list: List[List[float]]
    ) -> List[List[float]]:
        """
        Computes betascores (methylation frequency) based on sites. Two lists
        :param raw_pos_list: positions for each llr for site-frequency
        :param raw_llr_list: one list of log-likelihood ratios per sample to be compared
        :return:
        """
        unique_pos = list(set().union(*[set(pos) for pos in raw_pos_list]))
        met_list = [
            [sum(llr > self.llr_threshold for s_pos, llr in zip(poss, llrs) if s_pos == pos) for pos in unique_pos]
            for poss, llrs in zip(raw_pos_list, raw_llr_list)
        ]
        nonambig_list = [
            [
                sum(abs(llr) > self.llr_threshold for s_pos, llr in zip(s_poss, llrs) if s_pos == pos)
                for pos in unique_pos
            ]
            for s_poss, llrs in zip(raw_pos_list, raw_llr_list)
        ]
        is_valid_idx = [True for _ in nonambig_list[0]]
        for na in nonambig_list:
            is_valid_idx = [a and b for a, b in zip(is_valid_idx, [n > 0 for n in na])]
        met_list = [[b for b, i in zip(met, is_valid_idx) if i] for met in met_list]
        nonambig_list = [[b for b, i in zip(nonambig, is_valid_idx) if i] for nonambig in nonambig_list]
        bs_list = [[m / na for m, na in zip(met, nonambig)] for met, nonambig in zip(met_list, nonambig_list)]
        return bs_list
    
    def compute_read_betascores(
        self, raw_read_list: List[List[Any]], raw_llr_list: List[List[float]]
    ) -> List[List[float]]:
        """
        Computes betascores (methylation frequency) based on reads. Two lists
        :param raw_read_list: readname for each llr for site-frequency
        :param raw_llr_list: one list of log-likelihood ratios per sample to be compared
        :return:
        """
        num_samples = len(raw_llr_list)
        bs_list = []
        for sample in range(num_samples):
            unique_reads = list(set(raw_read_list[sample]))
            if len(unique_reads) > 1000:
                # excessive coverage indicates weird mapping - we will subsample to not cause a huge memory spike
                unique_reads = [r for r in unique_reads if random.random() < 100 / len(unique_reads)]
            met_list = [
                sum(
                    llr > self.llr_threshold
                    for read, llr in zip(raw_read_list[sample], raw_llr_list[sample])
                    if read == r
                )
                for r in unique_reads
            ]
            nonambig_list = [
                sum(
                    abs(llr) > self.llr_threshold
                    for read, llr in zip(raw_read_list[sample], raw_llr_list[sample])
                    if read == r
                )
                for r in unique_reads
            ]
            bs = [m / na for m, na in zip(met_list, nonambig_list) if na > 0]
            bs_list.append(bs)
        return bs_list
    
    def compute_contingency_table(self, raw_llr_list: List[List[float]]) -> List[List[float]]:
        """
        Computes contingency table as for fisher exact test by thresholding llrs and counting
        methylation/unmethylated calls for each sample
        :param raw_llr_list: one list of log-likelihood ratios per sample to be compared
        :return: contingency table of shape (n_samples, 2)
        """
        num_samples = len(raw_llr_list)
        contingency_table = []
        for sample in range(num_samples):
            n_met = sum(llr > self.llr_threshold for llr in raw_llr_list[sample])
            n_called = sum(abs(llr) > self.llr_threshold for llr in raw_llr_list[sample])
            n_unmet = n_called - n_met
            contingency_table.append([n_met, n_unmet])
        return contingency_table
    
    def compute_posthoc_test(self, test_values):
        try:
            posthoc_pvalue_list = []
            for sample_one in range(len(test_values)):
                values_one = test_values[sample_one]
                if self.pvalue_method == "chi_squared":
                    # If the original test was a chi_squared test, we use fisher exact to compute post-hoc test
                    # on the corresponding rows of the contingency table
                    values_others = [sum(v[i] for j, v in enumerate(test_values) if j != sample_one) for i in range(2)]
                    _, pvalue = fisher_exact([values_one, values_others])
                elif self.pvalue_method == "KW":
                    # Merge sample values before and after
                    values_others = [v for vl in test_values[:sample_one] for v in vl]
                    values_others += [v for vl in test_values[sample_one + 1 :] for v in vl]
                    if abs(np.mean(values_others) - np.mean(values_one)) > 0.25:  # FIXME expose parameter
                        _, pvalue = mannwhitneyu(values_one, values_others)
                    else:
                        pvalue = 1
                else:
                    raise ValueError("Internal error: Attempted to compute post-hoc when not appropriate")
                
                posthoc_pvalue_list.append(pvalue)
        except ValueError:
            posthoc_pvalue_list = [1.0] * len(test_values)
        return posthoc_pvalue_list
    
    def compute_pvalue(self, interval, label_list, raw_llr_list, raw_pos_list, raw_reads_list):
        counters_to_increase = []
        res = OrderedDict()
        res["chromosome"] = interval.chr_name
        res["start"] = interval.start
        res["end"] = interval.end
        avg_coverage = [len(pos) / len(set(pos)) for pos in raw_pos_list]
        non_ambig_llr_count = [sum(1 for l in llrs if abs(l) > self.llr_threshold) for llrs in raw_llr_list]
        pos_count = [sum(1 for l in llrs if l > self.llr_threshold) for llrs in raw_llr_list]
        n_samples = sum(1 for c in non_ambig_llr_count if c > 0)
        overall_bs_list = [pos / total for pos, total in zip(pos_count, non_ambig_llr_count) if total > 0]
        
        # Collect median llr
        med_llr_list = [np.median(llrs) for llrs in raw_llr_list]
        
        # "Lazy load" this variable as a slight performance boon
        read_beta_scores = None
        
        if self.hypothesis == "llr_diff":
            test_values = raw_llr_list
        elif self.hypothesis == "bs_diff":
            if self.pvalue_method == "paired":
                test_values = self.compute_site_betascores(raw_pos_list, raw_llr_list)
            else:
                read_beta_scores = self.compute_read_betascores(raw_reads_list, raw_llr_list)
                test_values = read_beta_scores
        elif self.hypothesis == "count_dependency":
            test_values = self.compute_contingency_table(raw_llr_list)
        
        if n_samples < self.min_samples:
            # Not enough samples
            comment = "Insufficient samples"
            pvalue = np.nan
        elif all(len(vals) < 3 for vals in test_values) and self.pvalue_method == "paired":
            comment = "Insufficient coverage"
            pvalue = np.nan
        # Sufficient samples and effect size
        else:
            comment = "Valid"
        
        # Update counters result table
        counters_to_increase.append(comment)
        
        if len(overall_bs_list) > 0:
            difference = np.diff(overall_bs_list).tolist()
        else:
            difference = []
        if comment == "Valid":
            post_hoc_pvalues = []
            try:
                # Run stat test
                if self.pvalue_method == "KW":
                    statistics, pvalue = kruskal(*test_values)
                    if pvalue < self.pvalue_threshold:
                        post_hoc_pvalues = self.compute_posthoc_test(test_values)
                elif self.pvalue_method == "MW":
                    statistics, pvalue = mannwhitneyu(test_values[0], test_values[1], alternative="two-sided")
                elif self.pvalue_method == "paired":
                    statistics, pvalue = wilcoxon(test_values[0], test_values[1])
                elif self.pvalue_method == "fisher_exact":
                    statistics, pvalue = fisher_exact(test_values)
                elif self.pvalue_method == "chi_squared":
                    statistics, pvalue, _, _ = chi2_contingency(test_values)
                    if pvalue < self.pvalue_threshold:
                        post_hoc_pvalues = self.compute_posthoc_test(test_values)
            except ValueError:
                # This happens for example if all values are equal in mannwhitneyu
                pvalue = 1
            
            # Fix and categorize p-values
            if pvalue is np.nan or pvalue is None or pvalue > 1 or pvalue < 0:
                counters_to_increase.append("Sites with invalid pvalue")
            # Correct very low pvalues to minimal float size
            elif pvalue == 0:
                pvalue = self.min_pval
            
            # Compute statistic used for independent hypothesis weighting
            if self.do_independent_hypothesis_weighting:
                if read_beta_scores is None:
                    read_beta_scores = self.compute_read_betascores(raw_reads_list, raw_llr_list)
                ihw_weight = self.compute_ihw_weight(read_beta_scores)
            
            res["pvalue"] = pvalue
            res["adj_pvalue"] = np.nan
            res["n_samples"] = n_samples
            if self.do_independent_hypothesis_weighting:
                res["ihw_weight"] = ihw_weight
            res["labels"] = list_to_str(label_list)
            res["med_llr_list"] = list_to_str(med_llr_list)
            res["raw_llr_list"] = list_to_str(raw_llr_list)
            res["difference"] = list_to_str(difference)
            if self.pvalue_method in {"KW", "chi_squared"}:
                res["post_hoc_pvalues"] = list_to_str(post_hoc_pvalues)
            res["comment"] = comment
            res["raw_pos_list"] = list_to_str(raw_pos_list)
            res["avg_coverage"] = list_to_str(avg_coverage)
            res["unique_cpg_pos"] = len(set(itertools.chain.from_iterable(raw_pos_list)))
        else:
            res["pvalue"] = np.nan
            res["adj_pvalue"] = np.nan
            res["n_samples"] = 0
            if self.do_independent_hypothesis_weighting:
                res["ihw_weight"] = 0.0
            res["labels"] = "[]"
            res["med_llr_list"] = "[]"
            res["raw_llr_list"] = "[]"
            res["comment"] = comment
            res["difference"] = "[]"
            if self.pvalue_method in {"KW", "chi_squared"}:
                res["post_hoc_pvalues"] = "[]"
            res["raw_pos_list"] = "[]"
            res["avg_coverage"] = "[]"
            res["unique_cpg_pos"] = "[]"
        
        return res, counters_to_increase
    
    def __call__(self, interval):
        try:
            label_list, raw_llr_list, raw_pos_list, raw_read_list = self.loader.read_raw_llrs(interval)
            return self.compute_pvalue(interval, label_list, raw_llr_list, raw_pos_list, raw_read_list)
        except:
            import traceback
            
            print(traceback.format_exc())
            raise


def initializer(args: Dict):
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
    read_groups_key: str = None,
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
    paired_test: bool = False,
    worker_processes: int = 4,
    hypothesis: str = "bs_diff",
    do_independent_hypothesis_weighting: bool = False,
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
     * read_groups_key
         Key in h5 file containing read groups to be used. (optional)
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
     * paired_test
         Test with a paired test on beta scores instead of unpaired on llrs
     * only_tested_sites
         Do not include sites that were not tested because of insufficient samples or effect size in the report
     * worker_processes
         Number of processes to be launched
     * hypothesis
        "llr_diff" if the hypotheis is a shift in llrs, "bs_diff" if the hypothesis is a shift in mean read methylation
        rate, or "count_dependency" if the hypothesis is dependencies between groups in the contingency table of
        methylated/unmethylated calls
     * do_independent_hypothesis_weighting
         Whether to include independent hypothesis weighting in the p-value adjustment
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
    
    sample_id_list = MetH5Loader.interpret_sample_ids_from_arguments(sample_id_list, read_groups_key, h5_file_list)
    
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
        if hypothesis == "count_dependency":
            log.debug("Multiple comparison mode for count depdencies (Chi-squared test)")
            pvalue_method = "chi_squared"
            min_samples = 3
        else:
            pvalue_method = "KW"
            log.debug("Multiple comparison mode (Kruskal_Wallis test)")
            if min_samples < 3:
                log.debug("Automatically raise number of minimal samples to 3")
                min_samples = 3
    # 2 values = Mann_Withney test
    elif all_samples == 2:
        if min_samples:
            log.debug("No missing samples allowed for 2 samples comparison")
            min_samples = 2
        if paired_test:
            pvalue_method = "paired"
            log.debug("Paired comparison mode (Wilcoxon)")
        else:
            if hypothesis == "count_dependency":
                log.debug("Pairwise comparison mode for count depdencies (Fisher exact test)")
                pvalue_method = "fisher_exact"
            else:
                pvalue_method = "MW"
                log.debug("Pairwise comparison mode (Mann_Withney test)")
    else:
        raise pycoMethError("Meth_Comp needs at least 2 input files")
    
    log.warning("Opening H5 files")
    try:
        # Define StatsResults to collect valid sites and perform stats
        stats_results = StatsResults(
            pvalue_adj_method=pvalue_adj_method,
            pvalue_threshold=pvalue_threshold,
            only_tested_sites=only_tested_sites,
            do_independent_hypothesis_weighting=do_independent_hypothesis_weighting,
        )
        
        log.info("Starting asynchronous file parsing")
        with tqdm(
            total=num_intervals,
            unit=" intervals",
            unit_scale=True,
            desc="\tProgress",
            disable=not progress,
        ) as pbar:
            
            log.info("Launching %d worker processes" % worker_processes)
            if worker_processes == 1:
                initializer(
                    dict(
                        h5_read_groups_key=read_groups_key,
                        sample_id_list=sample_id_list,
                        h5_file_list=h5_file_list,
                        min_diff_llr=min_diff_llr,
                        min_samples=min_samples,
                        pvalue_method=pvalue_method,
                        min_num_reads_per_interval=min_num_reads_per_interval,
                        pvalue_threshold=pvalue_threshold,
                        hypothesis=hypothesis,
                        do_independent_hypothesis_weighting=do_independent_hypothesis_weighting,
                    )
                )
            else:
                pool = Pool(
                    worker_processes,
                    initializer=initializer,
                    initargs=[
                        dict(
                            h5_read_groups_key=read_groups_key,
                            sample_id_list=sample_id_list,
                            h5_file_list=h5_file_list,
                            min_diff_llr=min_diff_llr,
                            min_samples=min_samples,
                            pvalue_method=pvalue_method,
                            min_num_reads_per_interval=min_num_reads_per_interval,
                            pvalue_threshold=pvalue_threshold,
                            hypothesis=hypothesis,
                            do_independent_hypothesis_weighting=do_independent_hypothesis_weighting,
                        )
                    ],
                )
            
            # Continue reading lines from all files
            log.debug("Starting deep parsing")
            fp_done = 0
            
            # Init file writer
            with Comp_Writer(
                bed_fn=output_bed_fn,
                tsv_fn=output_tsv_fn,
                verbose=verbose,
                output_raw_lists=False,
                with_ihw_weight=do_independent_hypothesis_weighting,
                with_posthoc_test=pvalue_method in {"KW", "chi_squared"},
            ) as writer:
                try:
                    
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
                        if worker_processes == 1:
                            callback(worker(interval))
                        else:
                            ar = pool.apply_async(
                                worker_function,
                                args=[interval],
                                callback=callback,
                                error_callback=error_callback,
                            )
                            async_results.append(ar)
                    if worker_processes > 1:
                        for i, ar in enumerate(async_results):
                            if abort:
                                break
                            ar.wait(timeout=3 * 24 * 3600)
                except:
                    writer.abort()
                    raise
            # Exit condition
            if not stats_results.res_list:
                log.info("No valid p-Value could be computed")
            else:
                # Convert results to dataframe and correct pvalues for multiple tests
                log.info("Adjust pvalues")
                stats_results.multitest_adjust()
                
                rewriter = Comp_ReWriter([f for f in (output_bed_fn, output_tsv_fn) if f is not None])
                rewriter.write_adjusted_pvalues(stats_results.res_list)
    finally:
        # Print counters
        log_dict(stats_results.counter, log.info, "Results summary")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~StatsResults HELPER CLASS~~~~~~~~~~~~~~~~~~~~~~~~~~~#


class DictList:
    def __init__(self, universal_key):
        self.dict = {}
        self.universal_key = universal_key
        pass
    
    def __len__(self):
        if len(self.dict) == 0:
            return 0
        
        return len(self.dict[self.universal_key])
    
    def __getitem__(self, i):
        class ListAccessor:
            def __init__(innerSelf):
                innerSelf.i = i
            
            def __getitem__(innerSelf, key):
                return self.dict[key][innerSelf.i]
            
            def __setitem__(innerSelf, key, value):
                self.dict[key][innerSelf.i] = value
        
        return ListAccessor()
    
    def __setitem__(self, i, d: Dict):
        for key, val in d.items():
            self.dict[key][i] = val
    
    def append(self, d):
        if len(self.dict) == 0:
            self.dict = {key: [val] for key, val in d.items()}
        else:
            if len(set(d.keys()).intersection(set(self.dict.keys()))) != len(d.keys()):
                raise ValueError("All keys must be present in all entries")
            for key, val in d.items():
                self.dict[key].append(val)
    
    def __iter__(self):
        for i in range(len(self)):
            yield self[i]


class StatsResults:
    def __init__(
        self,
        pvalue_adj_method="fdr_bh",
        pvalue_threshold=0.01,
        only_tested_sites=False,
        do_independent_hypothesis_weighting=True,
    ):
        """"""
        # Save self variables
        self.pvalue_adj_method = pvalue_adj_method
        self.pvalue_threshold = pvalue_threshold
        self.only_tested_sites = only_tested_sites
        self.do_independent_hypothesis_weighting = do_independent_hypothesis_weighting
        
        # Init self collections
        self.res_list =  DictList("pvalue")
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
        if "ihw_weight" in res:
            reduced_res["ihw_weight"] = res["ihw_weight"]
        self.res_list.append(reduced_res)
        return res
    
    def multitest_adjust(self):
        """"""
        # Collect non-nan pvalues
        pvalue_idx = []
        pvalue_list = []
        ihw_weight_list = []
        
        for i, res in enumerate(self.res_list):
            if not np.isnan(res["pvalue"]):
                pvalue_idx.append(i)
                pvalue_list.append(res["pvalue"])
                if self.do_independent_hypothesis_weighting:
                    ihw_weight_list.append(res["ihw_weight"])
        print(pvalue_list)
        # Adjust values
        if len(pvalue_list) == 0:
            return
        
        if self.do_independent_hypothesis_weighting:
            # Re-center weights (must average to 1) and is most faithful to FDR if scaled between 0 and 2
            mean_weight = sum(ihw_weight_list) / len(ihw_weight_list)
            ihw_weight_list = [w - mean_weight for w in ihw_weight_list]
            min_weight = max(abs(min(ihw_weight_list)), 0.01)
            ihw_weight_list = [max(w / min_weight + 1, 0.01) for w in ihw_weight_list]
            # Weight p-values
            pvalue_list = [p / w for p, w in zip(pvalue_list, ihw_weight_list)]
        
        adj_pvalue_list = multipletests(
            pvals=pvalue_list,
            alpha=self.pvalue_threshold,
            method=self.pvalue_adj_method,
        )[1]
        
        # add adjusted values to appropriate category
        for i, adj_pvalue in zip(pvalue_idx, adj_pvalue_list):
            
            # Fix and categorize p-values
            if adj_pvalue is np.nan or adj_pvalue is None or adj_pvalue > 1 or adj_pvalue < 0:
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
            if self.res_list[i]["comment"] == "Valid":
                # Overwriting comment, but only if it was a site that was tested
                # (not if it is a site that was excluded for coverage or other reasons)
                self.counter[comment] += 1
                self.res_list[i]["comment"] = comment
            self.res_list[i]["adj_pvalue"] = adj_pvalue


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~Comp_Writer HELPER CLASS~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Comp_Writer:
    """Extract data for valid sites and write to BED and/or TSV file."""
    
    def __init__(
        self,
        bed_fn=None,
        tsv_fn=None,
        verbose=True,
        output_raw_lists=False,
        with_ihw_weight=False,
        with_posthoc_test=False,
    ):
        """"""
        self.log = get_logger(name="Comp_Writer", verbose=verbose)
        self.bed_fn = bed_fn
        self.tsv_fn = tsv_fn
        self.output_raw_lists = output_raw_lists
        self.with_ihw_weight = with_ihw_weight
        self.with_posthoc_test = with_posthoc_test
        
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
        fp = gzip.open(self.bed_fn, "wt") if self.bed_fn.endswith(".gz") else open(self.bed_fn, "w")
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
        fp = gzip.open(self.tsv_fn, "wt") if self.tsv_fn.endswith(".gz") else open(self.tsv_fn, "w")
        # Write header line
        
        self.header = [
            "chromosome",
            "start",
            "end",
            "n_samples",
            "pvalue",
            "adj_pvalue",
            "unique_cpg_pos",
            "labels",
            "med_llr_list",
            "difference",
        ]
        if self.with_posthoc_test:
            self.header = self.header + ["post_hoc_pvalues"]
        if self.output_raw_lists:
            self.header = self.header + ["raw_llr_list", "raw_pos_list"]
        if self.with_ihw_weight:
            self.header = self.header + ["ihw_weight"]
        self.header = self.header + [
            "avg_coverage",
            "comment",
        ]
        fp.write(str_join(self.header, sep="\t", line_end="\n"))
        return fp
    
    def _write_tsv(self, res):
        """Write line to TSV file."""
        res_line = [res[k] for k in self.header]
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
                is_bed = filename.endswith("bed")
                with fileinput.input(filename, inplace=True) as fi_fp:
                    for res, line in zip((None, *res_list), fi_fp):
                        line = line.strip()
                        if res is None:
                            header = {k: i for i, k in enumerate(line.split(sep))}
                            print(line)
                        else:
                            updated_line = line.split(sep)
                            if is_bed:
                                if np.isnan(res["adj_pvalue"]) or res["adj_pvalue"] <= 0:
                                    score = 0
                                else:
                                    score = int(-np.log10(res["adj_pvalue"]))
                                
                                updated_line[4] = str(score)
                            else:
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
        raise pycoMethError("Invalid number of columns in read groups file (should be 2 or 3)")
    
    if not all([col in read_groups.columns for col in should_colnames]):
        raise pycoMethError("Invalid column names in read groups file (should be %s)" % should_colnames.join(", "))
    
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
        for line in bed:
            ct = coordgen(line.chrom, line.start, line.end)
            yield (ct)
