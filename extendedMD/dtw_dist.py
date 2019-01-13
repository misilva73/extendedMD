"""Compute DTW distance matrix"""
from dtaidistance import dtw
import numpy as np

def compute_dtw_dist_mat(ts_list, R=None):
    dist_matrix_vec = dtw.distance_matrix(ts_list, parallel=True, max_dist=R)
    dist_matrix = np.triu(dist_matrix_vec) + np.triu(dist_matrix_vec).T
    np.fill_diagonal(dist_matrix, 0)
    return dist_matrix