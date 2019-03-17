import math


def compute_motif_mdl_cost(members_dic_list, bs_len):
    """
    This function compute the MDL cost associated to the motif that generated the members in members_dic_list,
    as proposed by Tanaka et all

    :param members_dic_list: list of dictionaries related to the BS subsequences that belong to a motif
    :type members_dic_list: list of dic
    :param bs_len: list of the lengths of each BS sequence
    :type bs_len: list of int
    :return: mdl cost
    :rtype: float
    """
    split_bs_len_list = split_bs_len(members_dic_list, bs_len)
    mdl_cost = compute_segmentation_mdl_cost(split_bs_len_list)
    return mdl_cost


def split_bs_len(members_dic_list, bs_len):
    """
    This function splits the BS sequence based on the position of the motif members
    (i.e. subsequences in members_dic_list). This is a middle step in the MDL cost computation proposed by Tanaka et all

    :param members_dic_list: list of dictionaries related to the BS subsequences that belong to a motif
    :type members_dic_list: list of dic
    :param bs_len: list of the lengths of each BS sequence
    :type bs_len: list of int
    :return: list of split BS lengths
    :rtype: list of list of int
    """
    bs_position_list = [dic['bs_position'] for dic in members_dic_list]
    break_points = [pos[0] for pos in bs_position_list]
    subseq_size = len(bs_position_list[0])
    split_bs_len_list = []
    if break_points[0] > 0:
        first_seq = bs_len[0:break_points[0]]
        split_bs_len_list.append(first_seq)
    for i in range(len(break_points) - 1):
        pattern_seq = bs_len[break_points[i]:break_points[i] + subseq_size]
        split_bs_len_list.append(pattern_seq)
        next_seq = bs_len[break_points[i] + subseq_size:break_points[i + 1]]
        split_bs_len_list.append(next_seq)
    final_pattern_seq = bs_len[break_points[-1]:break_points[-1] + subseq_size]
    split_bs_len_list.append(final_pattern_seq)
    if break_points[-1] + subseq_size < len(bs_len):
        final_seq = bs_len[break_points[-1] + subseq_size:]
        split_bs_len_list.append(final_seq)
    return split_bs_len_list


def compute_segmentation_mdl_cost(seg_len_list):
    """
    This function computes the MDL cost associated with the lengths split BS sequence.
    This is a middle step in the MDL cost computation proposed by Tanaka et all
    as proposed by Tanaka et all

    :param seg_len_list: list of split BS lengths
    :type seg_len_list: list of list of int
    :return: mdl cost
    :rtype: float
    """
    par_cost_list = []
    data_cost_list = []
    for seg in seg_len_list:
        seg_len_sum = float(sum(seg))
        seg_par_cost = math.log2(seg_len_sum)
        par_cost_list.append(seg_par_cost)
        seg_data_cost_list = [-l*math.log2(l/seg_len_sum) for l in seg]
        seg_data_cost = sum(seg_data_cost_list)
        data_cost_list.append(seg_data_cost)
    par_cost = sum(par_cost_list)
    data_cost = sum(data_cost_list)
    seg_cost = len(seg_len_list)*math.log2(sum(sum(seg_len_list, [])))
    mdl_cost = round(par_cost + data_cost + seg_cost, 2)
    return mdl_cost
