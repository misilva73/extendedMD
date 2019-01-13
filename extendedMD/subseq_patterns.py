"""Functions to extract BS subsequences and BS patterns"""
import numpy as np


def get_bs_subsequences_list(bs_seq, bs_point, subseq_size):
    subseq_bs_list = []
    subseq_point_list = []
    for i in range(len(bs_seq) - subseq_size + 1):
        subseq_bs = bs_seq[i:(i+subseq_size)]
        subseq_bs_list.append(subseq_bs)
        subseq_point = list(set(np.concatenate(bs_point[i:(i+subseq_size)])))
        subseq_point_list.append(subseq_point)
    return subseq_bs_list, subseq_point_list


def get_all_subsequences_in_pattern(ts, subseq_bs_list, subseq_point_list, pattern):
    pattern_ts_list = []
    pattern_pos_list = []
    for i in range(len(subseq_bs_list)):
        if subseq_bs_list[i] == pattern:
            pattern_pos_list.append(i)
            pattern_ts = ts[subseq_point_list[i]]
            pattern_ts_list.append(pattern_ts)
    return pattern_ts_list, pattern_pos_list


def find_index_of_pattern_center_and_members(dist_mat, pattern_pos_list, R):
    """Finds the index of the motif's center and members"""
    count_subseq_list = []
    sum_dist_list = []
    for row in dist_mat:
        count_subseq = sum(row < R)
        count_subseq_list.append(count_subseq)
        sum_dist = sum(row[row < R])
        sum_dist_list.append(sum_dist)
    center_index, mean_dist = find_index_of_pattern_center(count_subseq_list, sum_dist_list)
    center_row = dist_mat[center_index]
    members_index = [index for index,value in enumerate(center_row) if value < R]
    # extract the center and members position in the subseq list
    center_pos = pattern_pos_list[center_index]
    members_pos = [pattern_pos_list[i] for i in members_index]
    return center_pos, members_pos, mean_dist


def find_index_of_pattern_center(count_subseq_list, sum_dist_list):
    """Finds the index of the motif's center, i.e., subseq with max count and min sum of dist"""
    count_subseq_list = np.array(count_subseq_list)
    sum_dist_list = np.array(sum_dist_list)
    max_count = max(count_subseq_list)
    center_index = np.argmin(sum_dist_list[count_subseq_list == max_count])
    mean_dist = np.min(sum_dist_list[count_subseq_list == max_count])/float(max_count)
    return center_index, mean_dist