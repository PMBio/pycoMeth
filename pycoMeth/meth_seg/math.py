import numpy as np
from scipy.stats import rankdata
import scipy

from typing import Tuple


def llr_to_p(llr, prior=0.5):
    """
    Convert log-likelihood ratios log(p(x|a)/p(x|~a)) to posterior
    probabilty p(a|x) given a prior p(a). For unbiased prediction,
    leave prior at 0.5
    """
    return 1 / (1 + np.exp(-llr) * (1 / prior - 1))


def p_to_llr(p, prior=0.5):
    """
    Converts the posterior probability p(a|x) into a log-likelihood ratio
    log(p(x|a)/p(x|~a)) given a prior pa(a)
    """
    return -np.log(prior * (1 - p) / (p * (1 - prior)))


def llr_to_uncertainty(llr, method="linear"):
    if method == "linear":
        p = llr_to_p(llr)
        return 0.5 - np.abs(0.5 - p)


def fdr_from_pvals(p_vals: np.ndarray) -> np.ndarray:
    """
    Computes FDR from p-values using the Benjamini-Hochberg method.
    :param p_vals: numpy array of p-values
    :return: numpy array of adjusted p-values
    """
    ranked_p_values = rankdata(p_vals)
    fdr = p_vals * len(p_vals) / ranked_p_values
    fdr[fdr > 1] = 1
    
    return fdr


def bs_from_llrs(llrs: np.ndarray, thres: float = 1, min_reads: int = 1) -> float:
    """
    Computes methylation beta score from a list of log-likelihood ratios
    :param llrs: Log-likelihood ratio array
    :param thres: threshold for absolute llr - excluding all llrs with an absolute llr lower than this threshold
                  (default: 1.0)
    :param min_reads: return np.nan if length of llrs after threshold filtering is less than min_reads (default: 1)
    :return: methylation beta score
    """
    llrs_used = llrs[np.abs(llrs) > thres]
    if len(llrs_used) < min_reads:
        return np.nan
    return (llrs_used > 0).sum() / len(llrs_used)

def __ensure_numpy(x) -> np.ndarray:
    if not isinstance(x, np.ndarray):
        x = np.array(x)
    return x

def nangmean(x: np.ndarray) -> float:
    """ Computes geometric mean while ignoring NaNs """
    x = __ensure_numpy(x)
    x = x[~np.isnan(x)]
    return scipy.stats.gmean(x)

def maxabs(x: np.ndarray) -> float:
    x = __ensure_numpy(x)
    """ Returns the value with the maximum magnitude """
    return x[np.unravel_index(np.argmax(np.abs(x)), x.shape)]

def compute_differential_methylation(
    llrs_a: np.ndarray, llrs_b: np.ndarray
) -> Tuple:

    # Paired test
    a_nan = llrs_a.copy()
    a_nan[a_nan == 0] = np.nan
    b_nan = llrs_b.copy()
    b_nan[b_nan == 0] = np.nan

    # Filtering sites for which both haplotypes have at least one read,
    # in order to avoid warnings
    good_sites = ((~np.isnan(a_nan)).sum(axis=0) > 0) & (
        (~np.isnan(b_nan)).sum(axis=0) > 0
    )
    a_nan = a_nan[:, good_sites]
    b_nan = b_nan[:, good_sites]

    pp = scipy.stats.ttest_rel(np.nanmean(a_nan, axis=0), np.nanmean(b_nan, axis=0))

    # Unpaired test
    up = scipy.stats.mannwhitneyu(llrs_a[llrs_a != 0], llrs_b[llrs_b != 0])

    if np.isnan(up[0]):
        # Workaround because ttest_ind returns (nan, nan) if it fails
        return None, None

    return up[1], pp[1]