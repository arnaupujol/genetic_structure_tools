#This module defines statistical operations for data analysis

import numpy as np
import scipy.special as sci_esp
from genomic_tools.relatedness import relatedness_mat

def mat_mean(mat):
    """
    This method computes the mean value of a diagonally
    symmetric matrix, excluding the diagonal
    matrix and one of the symmetric halfs.

    Parameters:
    -----------
    mat: np.ndarray
        Matrix with diagonaly symmetry.

    Returns:
    --------
    mean: float
        Mean value of matrix.
    """
    mat_len = len(mat)
    total = .5*(mat_len**2. - mat_len)
    mean = .0
    for i in range(mat_len):
        for j in range(mat_len):
            if j > i:
                mean += mat[i,j]/total
    return mean

def mat_vals(mat, mask = None, diag = True):
    """
    This method returns the values in a
    diagonally symmetric matrix excluding
    the diagonal terms and the elements masked.

    Parameters:
    -----------
    mat: np.ndarray
        Matrix with diagonaly symmetry.
    mask: np.array
        Mask selecting the elements in mat used.
    diag: bool
        It specifies whether the matrix is diagonally
        symmetric and repeated cases are excluded

    Returns:
    --------
    vals: np.array
        array of all the non-diagonal values
        of the matrix.
    """
    vals = []
    mat_len = mat.shape
    if mat_len[0] != mat_len[1]:
        diag = False
    if mask is None:
        if diag:
            mask = np.ones(mat_len[0], dtype = bool)
        else:
            mask = np.ones(mat_len, dtype = bool)
    for ii in range(mat_len[0]):
        for jj in range(mat_len[1]):
            if diag:
                if jj > ii and mask[ii] and mask[jj]:
                    vals.append(mat[ii,jj])
            else:
                if mask[ii,jj]:
                    vals.append(mat[ii,jj])
    return np.array(vals)

def mat_mean_err(mat, jk_num = 20, jk_num_2 = 20, rand_order = True, mask = None, diag = True):
    """
    This method measures the mean and Jack-Knife (JK) error of
    the values in a diagonally symmetric matrix excluding the
    diagonal terms.

    Parameters:
    -----------
    mat: np.ndarray
        Matrix with diagonaly symmetry.
    jk_num: int
        Number of JK subsamples to identify
    jk_num_2: int
        Number of JK subsamples for 2nd to identify when diag is False
    rand_order: bool
        If True, the indeces are assigned in a random order
    mask: np.array
        Mask selecting the elements in mat used
    diag: bool
        It specifies whether the matrix is diagonally
        symmetric and repeated cases are excluded

    Returns:
    --------
    mean_val: float
        Mean value of the array
    err: float
        JK error of the mean
    """
    mean_val = np.mean(mat_vals(mat, diag = diag))
    jk_ids = get_jk_indeces_1d(mat[0], jk_num, rand_order)
    if diag:
        mean_arr_jk = np.array([np.mean(mat_vals(mat, jk_ids != i, diag = diag)) for i in range(jk_num)])
    else:
        mean_arr_jk = []
        jk_ids_2 = get_jk_indeces_1d(mat, jk_num_2, rand_order)
        for i in range(jk_num):
            for j in range(jk_num_2):
                mask = np.ones_like(mat, dtype=bool)
                mask[:,jk_ids == i] = False
                mask[jk_ids_2 == j] = False
                mean_arr_jk.append(np.mean(mat[mask]))
        mean_arr_jk = np.array(mean_arr_jk)
    err = jack_knife(np.array([mean_val]), mean_arr_jk)
    if not diag:
        err = jk_2d_corr(err, jk_num, jk_num_2)
    return mean_val, err

def get_jk_indeces_1d(array, jk_num, rand_order = True):
    """
    This method assigns equally distributed indeces to the elements of an array.

    Parameters:
    -----------
    array: np.array
        Data array
    jk_num: int
        Number of JK subsamples to identify
    rand_order: bool
        If True, the indeces are assigned in a random order

    Returns:
    --------
    jk_indeces: np.array
        Array assigning an index (from 0 to jk_num - 1) to
        each of the data elements
    """
    ratio = int(len(array)/jk_num)
    res = int(len(array)%jk_num > 0)
    jk_indeces = (np.arange(len(array), dtype = int)/ratio).astype(int)
    jk_indeces[-res:] = np.random.randint(jk_num, size = res)
    np.random.shuffle(jk_indeces)
    return jk_indeces

def mean_err(array, jk_num = 50, rand_order = True):
    """
    This method measures the mean and Jack-Knife (JK) error of
    the values in an array

    Parameters:
    -----------
    array: np.array
        Data array
    jk_num: int
        Number of JK subsamples to identify
    rand_order: bool
        If True, the indeces are assigned in a random order

    Returns:
    --------
    mean_val: float
        Mean value of the array
    err: float
        JK error of the mean
    """
    mean_val = np.mean(array)
    jk_ids = get_jk_indeces_1d(array, jk_num, rand_order)
    mean_arr_jk = np.array([np.mean(array[jk_ids != i]) for i in range(jk_num)])
    err = jack_knife(np.array([mean_val]), mean_arr_jk)
    return mean_val, err

def jack_knife(var, jk_var):
    """
    This method gives the Jack-Knife error of var from the jk_var subsamples.

    Parameters:
    -----------
    var: float
        The mean value of the variable
    jk_var: np.ndarray
        The variable from the subsamples. The shape of the jk_var must be (jk subsamples, bins)

    Returns:
    --------
    jk_err: float
        The JK error of var.
    """
    jk_dim = np.prod(jk_var.shape)
    err = (jk_dim - 1.)/jk_dim * (jk_var - var)**2.
    jk_err = np.sqrt(np.sum(err, axis = 0))
    return jk_err


def sig2pow(sig):
    """
    This method returns the confidence interval
    of a distribution from a sigma factor.

    Parameters:
    -----------
    sig: float
        Value indicating how many sigmas away the
        signal is.

    Returns:
    --------
    p: float
        Power indicating interval of confidence.
    """
    return sci_esp.erf(sig/np.sqrt(2.))

def count_cases(variable, ignore_null = True):
    """This method counts the appearences of each case
    or label in a list of values.

    Parameters:
    -----------
    variable: (pd.Series, np.array)
        List of values or labels in the variable
    ignore_null: bool, default is True
        If True, is does not include the null cases

    Returns:
    --------
    cases: list
        List of appearing cases
    counts: list
        Number of appearences per case
    """
    cases = []
    counts = []
    #If we include the null cases, we add them here
    if not ignore_null:
        cases.append('None')
        counts.append(np.sum(variable.isnull()))
    #We take only the no null values and count their cases
    for case in variable[variable.notnull()].unique():
        cases.append(case)
        counts.append(np.sum(variable.notnull()*(variable == case)))
    return cases, counts

def get_percent(mat, p, mask = None, diag = True):
    """This method returns the percentile value for all
    the values within a matrix.

    Parameters:
    -----------
    mat: np.ndarray
        Matrix of values.
    p: float
        It defines the percentile to obtain
    mask: np.array
        Mask selecting the elements in mat used.
    diag: bool
        It specifies whether the matrix is diagonally
        symmetric and repeated cases are excluded

    Returns:
    --------
    val: float
        value corresponding to the percentile
    """
    return np.percentile(mat_vals(mat, mask = mask, diag = diag), p)

def percentile_change(mat_pre, mat_post, p, mask = None, diag = True, inverse = False):
    """
    This method measures the change in proportion of populaiton
    within a precentile defined on the first population.

    Parameters:
    -----------
    mat_pre: np.ndarray
        Matrix of values of the first population
    mat_post: np.ndarray
        Matrix of values of the second population
    p: float
        It defines the percentile to obtain from the first population.
    mask: np.array
        Mask selecting the elements in mat used
    diag: bool
        It specifies whether the matrix is diagonally
        symmetric and repeated cases are excluded
    inverse: bool
        If True, the fractions are obtained from values below the percentile.

    Returns:
    --------
    pval: float
        Value corresponding to the percentile in mat_pre
    frac_pre: float
        Fraction of values > pval on mat_pre
    frac_post: float
        Fraction of values > pval on mat_post
    """
    pval = get_percent(mat_pre, p)
    print("Value corresponding to the "+str(p)+" percentile on the first population: " +str(pval))
    #fraction of values higher than pval for first population
    if inverse:
        frac_pre = np.sum(mat_vals(mat_pre) < pval)/float(mat_vals(mat_pre).shape[0])
    else:
        frac_pre = np.sum(mat_vals(mat_pre) > pval)/float(mat_vals(mat_pre).shape[0])
    print("Fraction of values inside this range for the first population: " + str(frac_pre))
    #fraction of values higher than pval for second population
    if inverse:
        frac_post = np.sum(mat_vals(mat_post) < pval)/float(mat_vals(mat_post).shape[0])
    else:
        frac_post = np.sum(mat_vals(mat_post) > pval)/float(mat_vals(mat_post).shape[0])
    print("Fraction of values inside this range for the second population: " + str(frac_post))
    return pval, frac_pre, frac_post

def get_frac(data, pval, lower_than = False):
    """This method estimates the fraction of cases above a certain value.

    Parameters:
    -----------
    data: np.array
        The data values
    pval: float
        Value defining the threshold
    lower_than: bool
        If true, the fraction is calculated for cases lower
        than pval instead of higher

    Returns:
    --------
    frac: float
        Fraction of population inside the threshold
    """
    if lower_than:
        frac = np.sum(data < pval)/float(data.shape[0])
    else:
        frac = np.sum(data > pval)/float(data.shape[0])
    return frac

def get_frac_err(mat, pval, lower_than = False, jk_num = 20, rand_order = True, mask = None, diag = True):
    """This method estimates the fraction of cases above a certain value and its error.

    Parameters:
    -----------
    mat: np.ndarray
        Matrix with pairwise values
    pval: float
        Value defining the threshold
    lower_than: bool
        If true, the fraction is calculated for cases lower
        than pval instead of higher
    jk_num: int
        Number of Jack-Knife subsamples used to obtain the error
    rand_order: bool
        It specifies whether the Jack-knife subsamples are randomly
        ordered from the data
    mask: np.ndarray
        Mask specifying the elements used from mat
    diag: bool
        It specifies whether the matrix is diagonally
        symmetric and repeated cases are excluded

    Returns:
    --------
    frac: float
        Fraction of population inside the threshold
    err: float
        Error on the estimation of the fraction
    """
    frac = get_frac(mat_vals(mat, mask = mask, diag = diag), pval, lower_than = lower_than)
    jk_ids = get_jk_indeces_1d(mat[0], jk_num, rand_order)
    frac_jk = np.array([get_frac(mat_vals(mat, jk_ids != i, diag = diag), pval, lower_than = lower_than) for i in range(jk_num)])
    err = jack_knife(np.array([frac]), frac_jk)
    return frac, err

def get_diff_err_pow(mean_1, err_1, mean_2, err_2, verbose = True):
    """This method calculates the difference between two
    measurements, the errors and the power of the change.

    Parameters:
    -----------
    mean_1: float
        First measurement
    err_1: float
        Error on the first measurement
    mean_2: float
        Second measurement
    err_2: float
        Error on the second measurement
    verbose: bool
        If True, it prints the results (default True)

    Returns:
    diff_means: float
        Difference between the measurements
    diff_err: float
        Error on the difference from error propagation
    power: float
        Power of the difference being different than zero.
    """
    diff_means = mean_2 - mean_1
    diff_err = np.sqrt(err_1**2 + err_2**2)
    power = 1 - sig2pow(np.abs(diff_means)/diff_err)
    if verbose:
        print ("Difference in measurements: " + str(round(diff_means, 4)) + "+/-" + str(round(diff_err, 4)) + ', p = ' + str(round(power, 4)))
    return diff_means, diff_err, power
