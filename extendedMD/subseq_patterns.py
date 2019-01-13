"""Functions to extract BS subsequences and BS patterns"""
import numpy as np


def get_bs_subsequences_list(bs_seq, bs_point, subseq_size):
    # initialize lists to save the subsequences
    subseq_bs_list = []
    subseq_point_list = []
    for i in range(len(bs_seq) - subseq_size + 1):
        # extract the bs subsequence and append it to the list
        subseq_bs = bs_seq[i:(i+subseq_size)]
        subseq_bs_list.append(subseq_bs)
        # extract the pointers of the bs subsequence and append it to the list
        subseq_point = list(set(np.concatenate(bs_point[i:(i+subseq_size)])))
        subseq_point_list.append(subseq_point)
    return subseq_bs_list, subseq_point_list


def get_all_subsequences_in_pattern(ts, subseq_bs_list, subseq_point_list, pattern):
    # initialize list to save the subsequences that bellong to the pattern
    pattern_ts_list = []
    pattern_pos_list = []
    # for each bs subsequence:
    for i in range(len(subseq_bs_list)):
        # if the subsequence bellongs to the pattern:
        if subseq_bs_list[i] == pattern:
            # save the position of the subsequence in the list
            pattern_pos_list.append(i)
            # extract the initial ts data of that subsequence and append it to the list
            pattern_ts = ts[subseq_point_list[i]]
            pattern_ts_list.append(pattern_ts)
    return pattern_ts_list, pattern_pos_list


def find_index_of_pattern_center_and_members(dist_mat, pattern_pos_list, R):
    # initialize lists to save row counters
    count_subseq_list = []
    sum_dist_list = []
    # for each row in the dist mat:
    for row in dist_mat:
        # compute the count of subseq with a distance < R and append it to the list
        count_subseq = sum(row < R)
        count_subseq_list.append(count_subseq)
        # sum the distances of those subseq and append it to the list
        sum_dist = sum(row[row < R])
        sum_dist_list.append(sum_dist)
    # convert lists to arrays (needed for the argmin)
    count_subseq_list = np.array(count_subseq_list)
    sum_dist_list = np.array(sum_dist_list)
    # find the maximum count of subseq
    max_num = max(count_subseq_list)
    # find the center, i.e., subseq with max count and min sum of dist
    center_index = np.argmin(sum_dist_list[count_subseq_list == max_num])
    # extract the dist row of the center
    center_row = dist_mat[center_index]
    # extract the members
    members_index = [index for index,value in enumerate(center_row) if value < R]
    #extract the center and members position in the subseq list
    center_pos = pattern_pos_list[center_index]
    members_pos = [pattern_pos_list[i] for i in members_index]
    return center_pos, members_pos