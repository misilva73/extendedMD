import numpy as np


def find_pattern_center_and_members(pattern_dic_list, r):
    """
    This function finds the all the members of a motif (all non-overlapping subsequences with a distance to the motif
    center less than R) and and the motif's center (the subsequence with the maximal count of members and minimal mean
    distance)

    :param pattern_dic_list: list of dictionaries related to all the subsequences in the pattern
    :type pattern_dic_list: list of dic
    :param r: maximum distance to the center of the motif
    :type r: float
    :return:
        - center_dic (:py:class:`dic`) - dictionary of the motif's center subsequence
        - members_dic_list (:py:class:`list of dic`) - list of dictionaries related to all the motif's members
        - min_mean_dist (:py:class:`float`) - mean distance between the motif's center and all the motif's members
    """
    center_dic = pattern_dic_list[0]
    members_dic_list = find_all_pruned_members(0, pattern_dic_list, r)
    max_members_count = len(members_dic_list)
    min_mean_dist = np.mean([dist for dist in center_dic['dist_vec'] if dist < r])
    for candidate_index, candidate_dic in enumerate(pattern_dic_list[1:], 1):
        candidate_members_dic_list = find_all_pruned_members(candidate_index, pattern_dic_list, r)
        members_count = len(candidate_members_dic_list)
        mean_dist = np.mean([dist for dist in candidate_dic['dist_vec'] if dist < r])
        if members_count > max_members_count:
            center_dic = candidate_dic
            members_dic_list = candidate_members_dic_list
            max_members_count = members_count
            min_mean_dist = mean_dist
        elif (members_count == max_members_count) & (mean_dist < min_mean_dist):
            center_dic = candidate_dic
            members_dic_list = candidate_members_dic_list
            max_members_count = members_count
            min_mean_dist = mean_dist
        else:
            continue
    return center_dic, members_dic_list, min_mean_dist


def find_all_pruned_members(center_index, pattern_dic_list, r):
    """
    This function finds all the motif members assuming that the center of the motif is the BS sequence with index
    center_index. It returns the dictionaries related to the all the non-overlapping BS subsequences that have a
    distance to the center subsequence lower than R

    :param center_index: index of the center motif in the pattern_dic_list
    :type center_index: int
    :param pattern_dic_list: list of dictionaries related to all the subsequences in the pattern
    :type pattern_dic_list: list of dic
    :param r: maximum distance to the center
    :type r: float
    :return: list of dictionaries of the motif members
    :rtype: list of dic
    """
    center_dic = pattern_dic_list[center_index]
    unpruned_members_list = [pattern_dic_list[i] for i, dist in enumerate(center_dic['dist_vec']) if dist < r]
    if len(unpruned_members_list) < 2:
        return unpruned_members_list
    else:
        pruned_member_list = [unpruned_members_list[0]]
        for member_dic in unpruned_members_list[1:]:
            last_prunned_member = pruned_member_list[-1]
            if lists_overlap(last_prunned_member['pointers'], member_dic['pointers']):
                dist_in = last_prunned_member['dist_vec'][center_index]
                dist_out = member_dic['dist_vec'][center_index]
                if dist_in > dist_out:
                    pruned_member_list = pruned_member_list[:-1]
                    pruned_member_list.append(member_dic)
                else:
                    continue
            else:
                pruned_member_list.append(member_dic)
    return pruned_member_list


def lists_overlap(l1, l2):
    """
    This functions return whether two lists overlap

    :param l1: list 1
    :type l1: list
    :param l2: list 2
    :type l2: list
    :return: whether the list overlap. If TRUE, then the list overlap
    :rtype: bool
    """
    intersection_set = set(l1).intersection(set(l2))
    overlaping_test = (len(intersection_set) > 0)
    return overlaping_test
