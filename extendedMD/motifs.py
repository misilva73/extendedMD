from extendedMD.dtwdist import add_distance_vec_to_pattern_dic
from extendedMD.bs import get_bs_subsequences_dic_list
from extendedMD.patterns import find_pattern_center_and_members
from extendedMD.mdl import compute_motif_mdl_cost


def find_all_motif_candidates(ts, bs_seq, bs_len, bs_pointers, r):
    """
    This function finds all the motif candidates of all sizes and saves each motif instance in a dictionary. The result
    is a list of dictionaries.

    :param ts: original 1-d time-series
    :type ts: 1d array
    :param bs_seq: list of modified sax words (i.e. BS sequence)
    :type bs_seq: list of str
    :param bs_len: list of the lengths of each BS sequence
    :type bs_len: list of int
    :param bs_pointers: list of pointers to the original time-series
    :type bs_pointers: list of list of int
    :param r: maximum distance to the center of the motif
    :type r: float
    :return: motif_dic_list
    :rtype: list of dic
    """
    motif_dic_list = []
    subseq_size = 1
    while True:
        bs_subseq_dic_list = get_bs_subsequences_dic_list(ts, bs_seq, bs_pointers, subseq_size)
        subseq_bs_list = [dic['pattern'] for dic in bs_subseq_dic_list]
        subseq_bs_set = [list(item) for item in set(tuple(row) for row in subseq_bs_list)]
        # if there are no repeating patterns of size subseq_size, then there are no more motifs and the loop ends
        if len(subseq_bs_list) == len(subseq_bs_set):
            break
        else:
            for pattern in subseq_bs_set:
                # if the pattern appears only once, then it is not a motif
                if subseq_bs_list.count(pattern) < 2:
                    continue
                else:
                    motif_dic = build_motif_dic_from_pattern(bs_subseq_dic_list, bs_len, pattern, r)
                    # if the pattern only has one member, then it is not a motif
                    if motif_dic == {}:
                        continue
                    else:
                        motif_dic_list.append(motif_dic)
            print('Motif candidates of size {} successfully extracted'.format(subseq_size))
            subseq_size = subseq_size + 1
    return motif_dic_list


def build_motif_dic_from_pattern(bs_subseq_dic_list, bs_len, pattern, r):
    """
    This function builds the motif dictionary related to the pattern given as input

    :param bs_subseq_dic_list: list of dictionaries where each dic represents a single BS subsequence
    :type bs_subseq_dic_list: list of dic
    :param bs_len: list of the lengths of each BS sequence
    :type bs_len: list of int
    :param pattern: list of sax words that made up the pattern
    :type pattern: list of str
    :param r: maximum distance to the center of the motif
    :type r: float
    :return: motif_dic
    :rtype: dic
    """
    initial_pattern_dic_list = [dic for dic in bs_subseq_dic_list if dic['pattern'] == pattern]
    pattern_dic_list = add_distance_vec_to_pattern_dic(initial_pattern_dic_list, r)
    center_dic, members_dic_list, mean_dist = find_pattern_center_and_members(pattern_dic_list, r)
    # if the pattern only has one member, then it is not a motif
    if len(members_dic_list) < 2:
        return {}
    else:
        mdl_cost = compute_motif_mdl_cost(members_dic_list, bs_len)
        motif_dic = transform_members_list_into_motif_dic(center_dic, members_dic_list, mdl_cost, round(mean_dist, 2))
        return motif_dic


def transform_members_list_into_motif_dic(center_dic, members_dic_list, mdl_cost, mean_dist):
    """
    This functions creates the dictionary related to the motif composed by the subsequences in members_dic_list

    :param center_dic: dictionary related to the center of motif
    :type center_dic: dic
    :param members_dic_list: list of dictionaries related to all the subsequences that belong to the motif
    :type members_dic_list: list of dic
    :param mdl_cost: MDL cost of the motif (as defined by Tanaka et all)
    :type mdl_cost: float
    :param mean_dist: mean distance between all the pair of motif's members
    :type mean_dist: float
    :return: motif_dic
    :rtype: dic
    """
    members_ts_pointers = [dic['pointers'] for dic in members_dic_list]
    motif_dic = {
        'pattern': center_dic['pattern'],
        'mdl_cost': mdl_cost,
        'mean_dist': mean_dist,
        'members_ts_pointers': members_ts_pointers,
        'center_ts_pointers': center_dic['pointers']
    }
    return motif_dic
