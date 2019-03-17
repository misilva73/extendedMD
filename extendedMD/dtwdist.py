from dtaidistance import dtw
import numpy as np


def add_distance_vec_to_pattern_dic(pattern_dic_list, r=None):
    """
    This function returns a new pattern_dic_list where each dic has an extra entry with the vector of distances between
    the BS sequence related to that dic and all the other BS subsequences in the pattern

    :param pattern_dic_list: list of all the dictionaries related to BS sequences in the pattern
    :type pattern_dic_list: list of dic
    :param r: distance upper bound  (for DTW distance computation)
    :type r: float
    :return: new_pattern_dic_list
    :rtype: list of dic
    """
    pattern_ts_list = [dic['ts'] for dic in pattern_dic_list]
    dist_mat = compute_dtw_dist_mat(pattern_ts_list, r)
    new_pattern_dic_list = []
    for i, dic in enumerate(pattern_dic_list):
        new_dic = dic
        new_dic['dist_vec'] = dist_mat[i]
        new_pattern_dic_list.append(new_dic)
    return new_pattern_dic_list


def compute_dtw_dist_mat(ts_list, r=None):
    """
    This function computes the pairwise distance matrix of a list of time-series with Dynamic Time Warping distance.
    It is based on dtaidistance package

    :param ts_list: list of time-series to compare pairwise
    :type ts_list: list of 1D array
    :param r: distance upper bound - if distance is higher than R, then computation stops and the distance is set as inf
              this parameter serves merely for speeding up computation
    :type r: float
    :return: distance matrix
    :rtype: 2D array
    """
    dist_matrix_vec = dtw.distance_matrix(ts_list, parallel=True, max_dist=r)
    dist_matrix = np.triu(dist_matrix_vec) + np.triu(dist_matrix_vec).T
    np.fill_diagonal(dist_matrix, 0)
    return dist_matrix


def compute_dwt_dist_between_ts_and_list(single_ts, ts_list, r):
    """
    This function computes the list of DTW distances between a single time-series and a list of time-series.

    :param single_ts: single time-series
    :type single_ts: 1d array
    :param ts_list: list of time-series
    :type ts_list: list of 1d array
    :param r: distance upper bound  (for DTW distance computation)
    :type r: float
    :return: list of DTW distances
    :rtype: list of float
    """
    dist_list = []
    for ts in ts_list:
        dist = dtw.distance(single_ts, ts, max_dist=r)
        dist_list.append(dist)
    return dist_list
