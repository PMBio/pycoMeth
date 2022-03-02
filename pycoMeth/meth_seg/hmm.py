from typing import Tuple, Dict

import numpy as np
import scipy.optimize

from pycoMeth.meth_seg.emissions import EmissionLikelihoodFunction


def arraylogexpsum(x):
    ret = x[0]
    for i in range(1, len(x)):
        ret = logaddexp(ret, x[i])
    return ret  # if ret > -256 else -512


def logaddexp(a, b):
    ret = np.logaddexp(a, b)
    return ret


class SegmentationHMM:
    def __init__(
        self,
        max_segments: int,
        t_stay: float,
        t_move: float,
        e_fn: EmissionLikelihoodFunction,
        seg_penalty: float = 0,
        eps: np.float64 = np.exp(-512),
    ):
        self.eps = eps
        self.e_fn = e_fn
        
        self.num_segments = max_segments
        
        if t_stay + t_move + seg_penalty > 1:
            raise ValueError("t_stay + t_move + seg_penalty may not exceed 1")
        
        self.seg_penalty = np.array([(seg_penalty * i / max_segments) for i in range(max_segments)], dtype=np.float64)
        self.t_move = np.array([t_move - self.seg_penalty[i] for i in range(max_segments)], dtype=np.float64)
        self.t_stay = np.array([t_stay for i in range(max_segments)], dtype=np.float64)
        self.t_end = np.array([1 - self.t_move[i] - self.t_stay[i] for i in range(max_segments)], dtype=np.float64)
        self.t_move = np.log(self.t_move + eps)
        self.t_stay = np.log(self.t_stay + eps)
        self.t_end = np.log(self.t_end + eps)
    
    def t_fn(self, i, j):
        if i == j:
            return self.t_stay[i]
        if i == (j - 1):
            # Probability to move to the next state
            return self.t_move[i]  # + sim_penalty
        if j == (self.num_segments - 1):
            # Probability to go the last segment
            return self.t_end[i]
        
        raise RuntimeError("Transition %d to %d is not a valid transition in segmentation " "HMM " % (i, j))
    
    def forward(self, observations, obs_c):
        e_fn = self.e_fn.likelihood
        M = self.num_segments
        R = observations.shape[0]
        N = observations.shape[1]
        F = np.zeros((N, M), dtype=np.float) + self.eps
        F[0, 0] = 1 - F[0, :].sum() - self.eps
        F = np.log(F)
        start_prob = np.zeros(M) + self.eps
        start_prob[0] = 1
        start_prob = np.log(start_prob)
        
        for k in range(N):
            o = observations[:, k]
            for i in range(M):
                e = e_fn(i, o, obs_c)
                
                if k == 0:
                    F[k, i] = e + start_prob[i]
                    continue
                    
                    # Stay probability
                F[k, i] = e + F[k - 1, i] + self.t_fn(i, i)
                
                # Move probabilty
                if i > 0:
                    F[k, i] = logaddexp(F[k, i], e + F[k - 1, i - 1] + self.t_fn(i - 1, i))
                
                # End probability
                if i == M - 1:
                    # if end state we could have come from anywhere to the
                    # end state:
                    for j in range(M - 2):  # exclude last 2 because those were already
                        # handled above
                        F[k, i] = logaddexp(F[k, i], e + F[k - 1, j] + self.t_fn(j, i))
        evidence = F[-1, -1]
        return F, evidence
    
    def backward(self, observations, obs_c):
        e_fn = self.e_fn.likelihood
        R = observations.shape[0]
        M = self.num_segments
        N = observations.shape[1]
        B = np.zeros((N, M), dtype=np.float64) + self.eps
        B[-1, -1] = 1
        B = np.log(B)
        
        for k in range(N - 1, 0, -1):
            o = observations[:, k]
            k = k - 1
            for i in range(M):
                e_stay = e_fn(i, o, obs_c)
                
                if i == M - 1:
                    # If i is end state, we can only stay
                    B[k, i] = e_stay + B[k + 1, i] + self.t_fn(i, i)
                else:
                    e_move = e_fn(i + 1, o, obs_c)
                    # Move and stay probability
                    B[k, i] = logaddexp(
                        B[k + 1, i] + self.t_fn(i, i) + e_stay, B[k + 1, i + 1] + self.t_fn(i, i + 1) + e_move
                    )
                    if i < M - 2:
                        # End probability only if i<M-2 because otherwise it
                        # was covered by move or stay
                        e_end = e_fn(M - 1, o, obs_c)
                        B[k, i] = logaddexp(B[k, i], B[k + 1, M - 1] + self.t_fn(i, M - 1) + e_end)
        
        o = observations[:, 0]
        evidence = B[0, 0] + e_fn(0, o, obs_c)
        return B, evidence
    
    def viterbi(self, observations, obs_c):
        e_fn = self.e_fn.likelihood
        M = self.num_segments
        N = observations.shape[1]
        
        V = np.zeros((N, M), dtype=np.float64) + self.eps
        V[0, 0] = 1
        V = np.log(V)
        P = np.zeros((N, M), dtype=np.int32)
        
        start_prob = np.zeros(M) + self.eps
        start_prob[0] = 1
        start_prob = np.log(start_prob)
        
        for k in range(0, N - 1):
            o = observations[:, k]
            for i in range(M):
                e = e_fn(i, o, obs_c)
                
                if k == 0:
                    V[k, i] = np.max(e + start_prob[i])
                    continue
                
                p = np.zeros(M) - np.inf
                
                p[i] = V[k - 1, i] + self.t_fn(i, i)
                
                if i > 0:
                    p[i - 1] = V[k - 1, i - 1] + self.t_fn(i - 1, i)
                
                if i == M - 1:
                    # last two have been covered by stay and move
                    for j in range(M - 2):
                        p[j] = V[k - 1, j] + self.t_fn(j, i)
                p = e + p
                
                V[k, i] = np.max(p)
                P[k, i] = np.argmax(p)
            # Rescaling prevents underflow
            V[k, :] = V[k, :] - arraylogexpsum(V[k, :])
        V[-1, :] = np.log(self.eps)
        V[-1, -1] = np.max(V[-2, :])
        P[-1, -1] = np.argmax(V[-2, :])
        X = np.zeros(N, dtype=np.int32)
        Z = np.zeros(N, dtype=np.float32)
        X[N - 1] = M - 1
        Z[N - 1] = 0
        
        for k in range(N - 2, -1, -1):
            X[k] = P[k + 1, X[k + 1]]
            Z[k] = V[k + 1, X[k + 1]]
        
        return X, Z
    
    def MAP(self, posterior):
        M = self.num_segments
        N = posterior.shape[0]
        
        V = np.zeros((N, M), dtype=np.float) + self.eps
        V[0, 0] = 1
        V = np.log(V)
        P = np.zeros((N, M), dtype=np.int32)
        
        start_prob = np.zeros(M) + self.eps
        start_prob[0] = 1
        start_prob = np.log(start_prob)
        
        for k in range(0, N - 1):
            for i in range(M):
                e = posterior[k, i]
                
                if k == 0:
                    V[k, i] = np.max(e + start_prob[i])
                    continue
                
                p = np.zeros(M) - np.inf
                
                p[i] = V[k - 1, i] + self.t_fn(i, i)
                
                if i > 0:
                    p[i - 1] = V[k - 1, i - 1] + self.t_fn(i - 1, i)
                
                if i == M - 1:
                    for j in range(M - 2):  # last two have been covered by stay and
                        # move
                        p[j] = V[k - 1, j] + self.t_fn(j, i)
                p = e + p
                
                V[k, i] = np.max(p)
                P[k, i] = np.argmax(p)
            # Rescaling prevents underflow
            V[k, :] = V[k, :] - arraylogexpsum(V[k, :])
        V[-1, :] = np.log(self.eps)
        V[-1, -1] = np.max(V[-2, :])
        P[-1, -1] = np.argmax(V[-2, :])
        X = np.zeros(N, dtype=np.int32)
        Z = np.zeros(N, dtype=np.float32)
        X[N - 1] = M - 1
        Z[N - 1] = 0
        
        for k in range(N - 2, -1, -1):
            X[k] = P[k + 1, X[k + 1]]
            Z[k] = V[k + 1, X[k + 1]]
        
        return X, Z
    
    def baum_welch(
        self,
        observations: np.ndarray,
        tol: float = np.exp(-4),
        it_hook=None,
        samples: np.ndarray = None,
        verbose: bool = False,
    ) -> Tuple[Dict[int, np.ndarray], np.ndarray]:
        """
        Run the baum_welch algorithm, an expectation maximization algorithm,
        to find a segmentation of the methylation signal.

        Note that this algorithm is rather memory expensive. It will take
        O(CM) memory where C is the number of samples and M the maximum number
        of segments. If no samples are provided, C is equal to the number of
        reads, meaning the memory requirement grows with the read coverage.

        :param observations: a numpy array of shape RxN, where R is the
        number of reads and N is the number of gennomic positions (or CpG
        sites). The values need to be in the range (0,1) and are methylation
        predictions for the individual CpG sites. In order to speed up
        computation, missing predictions can be labeled with the value -1.
        This should lead to the same result as setting it to the value 0.5,
        but reduces the number of computations required significantly.
        :param tol: The absolute maximum difference in a parameter value that
        determines convergence. If the difference is below tol, the algorithm
        aborts
        :param it_hook: A function hook that will be called after each
        iteration. Takes the same parameters as the return value of this
        function
        :param samples: A 1-dimensional numpy array of length R, assigns each
        read to a sample id. Sample ids must be integer, and must start from
        0 and have no gaps.
        :return: tuple with estimated parameters and posteriors. Estimated
        paramater type depends on the given emission probability class.
        Posterior is of shape NxM and gives the posterior probability of each
        genomic site n being in each segment m
        """
        # Initial guess of parameters
        R = observations.shape[0]
        N = observations.shape[1]
        
        if N < 2:
            raise ValueError("Observations must contain at least 2 CpG-sites")
        if any((observations != -1).sum(axis=0) == 0):
            raise ValueError("Observations must not include reads with zero observations")
        if any((observations != -1).sum(axis=1) == 0):
            raise ValueError("Observations must not include sites with zero observations")
        
        if samples is None:
            # No samples, then we use the identity
            self.obs_c = np.arange(R)
        else:
            self.obs_c = samples
        
        C = len(set(self.obs_c))
        
        for it in range(100):
            F, f_evidence = self.forward(observations, self.obs_c)
            B, b_evidence = self.backward(observations, self.obs_c)
            # Sanity check: fwd and bwd algorithm should return same evidence
            if np.abs(f_evidence - b_evidence) > 10e-6:
                print("WARNING: forward evidence %f does not equal backward " "evidence %f." % (f_evidence, b_evidence))
            
            posterior = F + B - b_evidence
            
            # Maximize
            segment_p_new = {}
            
            for c in range(C):
                old_params = self.e_fn.get_cluster_params(c)
                to_minimize = self.e_fn.minimization_objective(observations[self.obs_c == c], np.exp(posterior))
                bounds = self.e_fn.get_param_bounds()
                
                estimated_p = scipy.optimize.minimize(to_minimize, old_params, method="SLSQP", bounds=bounds).x
                segment_p_new[c] = np.log(estimated_p)
            
            diff = self.e_fn.update_params(segment_p_new)
            
            segment_p = segment_p_new
            
            if it_hook is not None:
                it_hook(segment_p, posterior)
            
            if verbose:
                print("Iteration %d, parameter difference: %f" % (it, diff))
            if diff < tol:
                break
        
        return segment_p, posterior
