from dtaidistance import dtw
import numpy as np


def add_distance_vec_to_pattern_dic(pattern_dic_list, R=None):
    """
    This function returns a new pattern_dic_list where each dic has an extra entry with the vector of distances between
    the BS sequence related to that dic and all the other BS subsequences in the pattern
    :param pattern_dic_list: list of all the dictionaries related to BS sequences in the pattern
    :param R: distance upper bound  (for DTW distance computation)
    :return: new_pattern_dic_list
    """
    pattern_ts_list = [dic['ts'] for dic in pattern_dic_list]
    dist_mat = compute_dtw_dist_mat(pattern_ts_list, R)
    new_pattern_dic_list = []
    for i, dic in enumerate(pattern_dic_list):
        new_dic = dic
        new_dic['dist_vec'] = dist_mat[i]
        new_pattern_dic_list.append(new_dic)
    return new_pattern_dic_list


def compute_dtw_dist_mat(ts_list, R=None):
    """
    This function computes the pairwise distance matrix of a list of time-series with Dynamic Time Warping distance.
    It is based on dtaidistance package
    :param ts_list: list of time-series to compare pairwise
    :param R: distance upper bound - if distance is higher than R, then computation stops and the distance is set as inf
              this parameter serves merely for speeding up computation
    :return: dist_matrix
    """
    dist_matrix_vec = dtw.distance_matrix(ts_list, parallel=True, max_dist=R)
    dist_matrix = np.triu(dist_matrix_vec) + np.triu(dist_matrix_vec).T
    np.fill_diagonal(dist_matrix, 0)
    return dist_matrix