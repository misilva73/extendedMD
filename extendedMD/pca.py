import numpy as np


def extract_pca_ts(multi_dim_ts):
    """
    This function reduces a multi-dimensional time-series to 1d using PCA

    :param multi_dim_ts: multi-dimensional time-series as a pandas dataframe
    :type multi_dim_ts: DataFrame
    :return: 1-dimensional time-series
    :rtype: 1d array
    """
    # compute vector with the mean of each time series
    means_vec = multi_dim_ts.agg('mean').values
    # compute the covariance matrix of the multi-dim time-series data
    cov_mat = multi_dim_ts.cov().values
    # extract eigenvalues and eigenvectors (the PCs) of the covariance matrix
    e_val, e_vec = np.linalg.eigh(cov_mat)
    # get the eigenvector with the highest eigenvalue (i.e. the 1st PC)
    pc1_vec = e_vec[np.argmax(e_val)]
    # compute the 1-dim time series as the data's projection on the 1st PC
    ts_1d = np.dot((multi_dim_ts.values - means_vec), pc1_vec)
    return ts_1d
