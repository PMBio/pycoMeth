import numpy as np
from meth5.sparse_matrix import SparseMethylationMatrixContainer

from pycoMeth.meth_seg.math import llr_to_p
from pycoMeth.meth_seg.emissions import BernoulliPosterior
from pycoMeth.meth_seg.hmm import SegmentationHMM
from pycoMeth.meth_seg.postprocessing import cleanup_segmentation


def segment(sparse_matrix: SparseMethylationMatrixContainer, max_segments_per_window: int) -> np.ndarray:
    llrs = np.array(sparse_matrix.met_matrix.todense())
    obs = llr_to_p(llrs)
    samples = sparse_matrix.read_samples
    
    unique_samples = list(set(samples))
    
    id_sample_dict = {i: s for i, s in enumerate(unique_samples)}
    sample_id_dict = {v: k for k, v in id_sample_dict.items()}
    
    sample_ids = np.array([sample_id_dict[s] for s in samples])
    
    emission_lik = BernoulliPosterior(len(unique_samples), max_segments_per_window, prior_a=None)
    hmm = SegmentationHMM(
        max_segments=max_segments_per_window, t_stay=0.1, t_move=0.8, e_fn=emission_lik, eps=np.exp(-512)
    )
    segment_p, posterior = hmm.baum_welch(obs, tol=np.exp(-8), samples=sample_ids)
    
    segmentation, _ = hmm.MAP(posterior)
    
    segment_p_array = np.concatenate([v[np.newaxis, :] for v in segment_p.values()], axis=0)
    segmentation = cleanup_segmentation(segment_p_array, segmentation, min_parameter_diff=0.2)
    
    return segmentation
