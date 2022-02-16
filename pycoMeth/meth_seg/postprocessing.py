import numpy as np


def cleanup_segmentation(
    segment_p: np.ndarray, segments: np.ndarray, min_length: int = 5, min_parameter_diff=0.1
) -> np.ndarray:
    """
    Cleans up a segmentation by merging segments that are too close or too
    similar in their parameter space.
    :param segment_p: segment parameters. Shape (C,M) where C is the number
    of clusters and M is the number of segments
    :param segments: segmentation as returned from the viterbi or MAP algorithm.
    Shape is (N) where N is the number of genomic indices (i.e. CpG sites).
    :param min_length: minimum number of genomic loci per segment
    :param min_parameter_diff: minimum difference in parameter values between
    neighboring segments
    :return: a new segmentation of the same shape as segments
    """
    new_segments = segments.copy()
    lastx = 0
    for segment in sorted(list(set(segments))):
        if len(set(new_segments)) <= 1:
            # No need to go on if it's all just one segment
            break
        length = (new_segments == segment).sum()
        if segment == new_segments[-1]:
            candidate_replace = new_segments[new_segments != segment][-1]
        else:
            candidate_replace = segment + 1
        absdif = np.abs(segment_p[:, segment] - segment_p[:, candidate_replace]).max()
        if length < min_length or absdif < min_parameter_diff:
            new_segments[new_segments == segment] = candidate_replace
    
    return np.array(new_segments)
