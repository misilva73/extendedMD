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


def find_index_of_pattern_center_and_members(dist_mat, subseq_point_list, pattern_pos_list, R):
    """Finds the index of the motif's center and members"""
    pattern_point_list = [subseq_point_list[i] for i in pattern_pos_list]
    count_members_list = count_members_in_pattern(dist_mat, pattern_point_list, R)
    max_count = max(count_members_list)
    sum_dist_list = []
    for i in range(len(dist_mat)):
        if count_members_list[i]==max_count:
            row = dist_mat[i]
            pruned_members_index = prune_pattern_members(row, pattern_point_list, R)
            sum_dist = sum([row[i] for i in pruned_members_index])
        else:
            sum_dist = np.inf
        sum_dist_list.append(sum_dist)
    center_index = np.argmin(sum_dist_list)
    mean_dist = sum_dist_list[center_index]/float(count_members_list[center_index])
    center_row = dist_mat[center_index]
    members_index = prune_pattern_members(center_row, pattern_point_list, R)
    # extract the center and members position in the subseq list
    center_pos = pattern_pos_list[center_index]
    members_pos = [pattern_pos_list[i] for i in members_index]
    return center_pos, members_pos, mean_dist


def prune_pattern_members(dist_row, pattern_point_list, R):
    members_in_radius_index = [index for index,value in enumerate(dist_row) if value < R]
    dist_in_radius = [dist_row[i] for i in members_in_radius_index]
    members_with_no_overlap_index = [members_in_radius_index[0]]
    for i in range(1, len(members_in_radius_index)):
        if lists_overlap(pattern_point_list[members_in_radius_index[i]],
                         pattern_point_list[members_with_no_overlap_index[-1]]):
            dist_in = dist_in_radius[members_with_no_overlap_index[-1]]
            dist_out = dist_in_radius[members_in_radius_index[i]]
            if dist_in > dist_out:
                members_with_no_overlap_index = members_with_no_overlap_index[:-1]
                members_with_no_overlap_index.append(members_in_radius_index[i])
        else:
            members_with_no_overlap_index.append(members_in_radius_index[i])
    return members_with_no_overlap_index


def lists_overlap(l1, l2):
    intersection_set = set(l1).intersection(set(l2))
    overlaping_test = (len(intersection_set) > 0)
    return overlaping_test


def count_overlaps(point_list):
    overlap_count = 0
    for i in range(len(point_list)-1):
        if lists_overlap(point_list[i], point_list[i+1]):
            overlap_count = overlap_count + 1
    return overlap_count


def count_members_in_pattern(dist_mat, pattern_point_list, R):
    count_members_list = []
    for dist_row in dist_mat:
        members_point_list = [pattern_point_list[index] for index,value in enumerate(dist_row) if value < R]
        count_members = len(members_point_list) - count_overlaps(members_point_list)
        count_members_list.append(count_members)
    return count_members_list