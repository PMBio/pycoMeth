import math
from abc import ABC, abstractmethod
from typing import Dict

import numpy as np


class EmissionLikelihoodFunction(ABC):
    """
    Abstract class defining the interface for emission likelihood functions
    The types for parameters is not defined here (only that it's stored as Dict) so that
    the concrete implementation of the likelihood function may define it's own
    domain of parameters, as long as the update_params function takes the same structure
    of parameters as the get_params and get_cluster_params method return
    """
    
    @abstractmethod
    def update_params(self, params: Dict):
        """
        :param params: Dictionary where the key is the cluster id and the value contains
        the parameter(s)
        :type params: Dict
        :return:
        """
        pass
    
    @abstractmethod
    def likelihood(self, segment_index: int, observations: np.array):
        pass
    
    @abstractmethod
    def get_cluster_params(self, cluster: int):
        pass
    
    @abstractmethod
    def minimization_objective(self, observations: np.array, posterior: np.array):
        pass
    
    @abstractmethod
    def get_param_bounds(self):
        pass


class BernoulliPosterior(EmissionLikelihoodFunction):
    """
    This class models a Bernoulli likelihood with uncertain observations.
    Observations are provided as probabilities p(a|S), and the likelihood
    is parameterized with bernoulli-likelihoods p(a|mu) for each segment.
    Optionally, a gamma prior p(mu) can be defined as well.
    The likelihood is modeled as:
    L = (1-p(a|S)) * (1-mu) / (1-p(a)) + p(a|S) * mu / p(a)
    """
    
    def __init__(
        self,
        number_of_clusters,
        number_of_segments,
        prior_a: float = None,
        eps=np.exp(-512),
        initial_segment_p: np.array = None,
    ):
        # Just for type hinting and to get rid of warnings
        self.segment_p: np.array = None
        self.segment_prior: np.array = None
        self.prior_lognormfactor: float = 0
        
        self.eps: float = eps
        self.prior_a: float = prior_a
        
        if initial_segment_p is None:
            self.segment_p = np.zeros((number_of_clusters, number_of_segments))
            self.segment_p[:, ::2] = np.log(1 / 5)
            self.segment_p[:, 1::2] = np.log(4 / 5)
        else:
            if not (
                initial_segment_p.shape[0] == number_of_clusters and initial_segment_p.shape[0] == number_of_segments
            ):
                raise ValueError("Initial parameters must be of shape (" "number_of_clusters, number_of_segments)")
            self.segment_p = initial_segment_p.copy()
        
        if self.prior_a is not None:
            # Precompute the normfactor (in log space) of the prior gamma
            # distribution
            self.prior_lognormfactor = np.log(math.gamma(2 * prior_a) / (math.gamma(prior_a) ** 2))*2
            self.prior_pdf = scipy.stats.beta(self.prior_a, self.prior_a).pdf
        else:
            self.prior_pdf = None
    
    def update_prior(self):
        """
        Updates the prior beta distribution (in log space) based on
        self.segment_p and self.prior_a
        """
        self.segment_prior = self.segment_p * (self.prior_a - 1)
        
        self.segment_prior += np.log(1 - np.exp(self.segment_p) + self.eps) * (self.prior_a - 1)
        self.segment_prior += self.prior_lognormfactor
    
    def update_params(self, segment_p_list: Dict[int, np.ndarray]):
        """
        Updates bernoulli parameters
        :param segment_p_list: bernoulli parameters. A list with one numpy array
        per cluster
        :return: difference between old and new parameters to assess convergence
        """
        maxdiff = 0
        for sample in segment_p_list.keys():
            diff = np.max(np.abs(np.exp(self.segment_p[sample, :]) - np.exp(segment_p_list[sample])))
            maxdiff = max(diff, maxdiff)
            self.segment_p[sample, :] = segment_p_list[sample]
        
        if self.prior_a is not None:
            self.update_prior()
        return maxdiff
    
    def likelihood(self, segment_index: int, observations: np.array, observations_cluster_assignment: np.array):
        p = self.segment_p[observations_cluster_assignment, :]
        idx = observations != -1
        
        ret_a = np.log(-np.expm1(p[idx, segment_index]) + self.eps)
        ret_a += np.log1p(-observations[idx]) + np.log(0.5)
        
        ret_b = p[idx, segment_index]
        ret_b += np.log(observations[idx] + self.eps) + np.log(0.5)
        
        ret = np.logaddexp(ret_a, ret_b)
        if self.segment_prior is not None:
            ret += self.segment_prior[observations_cluster_assignment, :][idx, segment_index]
        return ret.sum()
    
    def minimization_objective(self, observations: np.array, posterior_exp: np.array):
        """
        Returns a curried function that only takes the candidate parameters
        mu and returns a minimization object (in this case the total
        likelihood p(S|mu',psi) given S and posteriors p(psi|mu))
        :param observations: observations numpy array
        :param posterior_exp: posterior of segmentation as estimated by the hmm
        - must be in linear space (not log space)
        :return: a function that takes parameters mu and returns likelihood
        """
        
        def curried_objective(x):
            # x is in log-space but these computations are easier in
            # lin-space
            # since there are a lot of additions going on
            # x = np.exp(x)
            ls = np.zeros(self.segment_p.shape[1])
            ps = np.zeros(self.segment_p.shape[1])
            for r in range(observations.shape[0]):
                o = observations[r, :]
                idx = o != -1
                o = o[idx]
                pki = posterior_exp[idx, :]
                
                l_a = np.outer((1 - o), (1 - x)) / 2
                l_b = np.outer(o, x) / 2
                l = np.log(l_a + l_b + self.eps)
                ls += (l * pki).sum(axis=0)
                ps += pki.sum(axis=0)
            
            ret = ls / (ps + self.eps)
            if self.prior_pdf is not None:
                prior = np.array([self.prior_pdf(xi) for xi in x])
                ret = ret * prior
            
            ret = ret.sum()
            
            return -ret
        
        return curried_objective
    
    def get_cluster_params(self, cluster):
        """
        Returns parameters for optimization for one cluster
        :param cluster: cluster index
        :return: parameters in linear space
        """
        # Clipping is needed due to a scipy bug
        return np.clip(np.exp(self.segment_p[cluster, :]), 0.02, 0.98)
    
    def get_params(self):
        # Clipping is needed due to a scipy bug
        return np.clip(np.exp(self.segment_p), 0.02, 0.98)
    
    def get_param_bounds(self):
        """
        :return: parameters bounds for optimization (interval [0,1])
        """
        return [(0.01, 0.99)] * self.segment_p.shape[1]
