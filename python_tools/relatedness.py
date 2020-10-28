"""
This module contains methods to compute metrics of genomic relatedness.
"""

import numpy as np
from scipy import spatial
import scipy.stats as sci_stats

def dist(vect_a, vect_b):
    """
    This method calculates the n-dimensional of two vectors.

    Parameters:
    -----------
    vect_a: np.ndarray
        Vector with high-dimensional position
    vect_b: np.ndarray
        Vector with high-dimensional position with the same size as vect_a

    Returns:
    --------
    d: float
        The distance modulus between the two vectors
    """
    d = np.linalg.norm(vect_a-vect_b)
    return d

def L1(vect_a, vect_b):
    """
    This method calculates the L1 norm of the difference between two vectors.

    Parameters:
    -----------
    vect_a: np.ndarray
        Vector with high-dimensional position
    vect_b: np.ndarray
        Vector with high-dimensional position with the same size as vect_a

    Returns:
    --------
    d: float
        The distance norm L1 between the two vectors
    """
    d = np.mean(np.abs(vect_a-vect_b))
    return d

def jaccard_dist(vect_a, vect_b):
    """
    This method calculates the Jaccard distance between two vectors.

    Parameters:
    -----------
    vect_a: np.ndarray
        Vector with high-dimensional position
    vect_b: np.ndarray
        Vector with high-dimensional position with the same size as vect_a

    Returns:
    --------
    d: float
        The Jaccard distance between the two vectors
    """
    PHS = np.sum((vect_a > 0)*(vect_b > 0))/float(np.sum((vect_a > 0)+(vect_b > 0)))
    return 1 - PHS



def relatedness_mat(data_1, data_2, method = 'L2'):
    """
    This method builds a matrix with the relatedness between two
    data samples for each pair comparison.

    Parameters:
    -----------
    data_1: np.ndarray
        Data with shape (samples, properties)
    data_2: np.ndarray
        Data with shape (samples, properties)
    method: str {'L2', 'L1', 'jaccard', 'bin_sharing', 'pcorr'}
        Metrics used to measure the relatedness (default is 'L2')

            'L2':
                L2-norm distance, Euclidean

            'L1':
                L1-norm distance

            'jaccard':
                Jaccard distance

            'bin_sharing':
                Binary sharing

            'pcorr':
                Pearson correlation coefficient

    Returns:
    --------
    rel_mat: np.ndarray
        Matrix of shape (data_1 samples, data_2 samples) with their
        relatedness
    """
    if method == 'L2': #TODO test
        rel_mat = spatial.distance_matrix(data_1, data_2)
    elif method == 'L1': #TODO test
        rel_mat = spatial.distance_matrix(data_1, data_2, p = 1)
    else:
        rel_mat = np.zeros((data_1.shape[0],data_2.shape[0]))
        for i in range(rel_mat.shape[0]):
            for j in range(rel_mat.shape[1]):
                if method == 'jaccard':
                    rel_mat[i,j] = jaccard_dist(data_1[i], data_2[j])
                if method == 'bin_sharing':
                    rel_mat[i,j] = jaccard_dist(data_1[i], data_2[j]) < 1
                if method == 'pcorr':
                    rel_mat[i,j] = sci_stats.pearsonr(data_1[i], data_2[j])[0]
    return rel_mat
