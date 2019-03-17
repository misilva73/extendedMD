import numpy as np
import math
import extendedMD.dtwdist as dtwdist


def prune_motifs_with_mdl(ts, motif_dic_list, r):
    """
    This function returns the most relevant motifs from the original list of motif extracted from the emd algorithm,
    based on the computed MDL cost and avoiding overlapping motifs

    :param ts: 1-dimensional time-series either resulting from the PCA method or the original 1-dimensional time-series
    :type ts: 1d array
    :param motif_dic_list: list of motif dictionaries returned from the emd algorithm
    :type motif_dic_list: list of dic
    :param r: maximum distance to the center of the motif
    :type r: float
    :return: list of dictionaries with the most relevant motifs. The list is ordered based on the MDL cost
    :rtype: list of dic
    """
    sorted_dic_list = sorted(motif_dic_list, key=lambda dic: dic['mdl_cost'])
    pruned_motif_dic_list = prune_motifs(ts, sorted_dic_list, r)
    return pruned_motif_dic_list


def prune_motifs_with_dist(ts, motif_dic_list, r, mdl_bins):
    """

    :param ts: 1-dimensional time-series either resulting from the PCA method or the original 1-dimensional time-series
    :type ts: 1d array
    :param motif_dic_list: list of motif dictionaries returned from the emd algorithm
    :type motif_dic_list: list of dic
    :param r: maximum distance to the center of the motif
    :type r: float
    :param mdl_bins: number of bins to break the MDL cost range
    :type mdl_bins: int
    :return: list of dictionaries with the most relevant motifs. The list is ordered based on MDL cost and motif's compactness
    :rtype: list of dic
    """
    mdl_sorted_dic_list = sorted(motif_dic_list, key=lambda dic: dic['mdl_cost'])
    step = math.floor(len(mdl_sorted_dic_list) / mdl_bins)
    dist_sorted_dic_list = []
    for i in range(mdl_bins):
        temp_dic_list = mdl_sorted_dic_list[i * step:(i + 1) * step]
        temp_dist_sorted_dic_list = sorted(temp_dic_list, key=lambda dic: dic['mean_dist'])
        dist_sorted_dic_list += temp_dist_sorted_dic_list
    if mdl_bins * step < len(mdl_sorted_dic_list):
        temp_dic_list = mdl_sorted_dic_list[mdl_bins * step:]
        temp_dist_sorted_dic_list = sorted(temp_dic_list, key=lambda dic: dic['mean_dist'])
        dist_sorted_dic_list += temp_dist_sorted_dic_list
    pruned_motif_dic_list = prune_motifs(ts, dist_sorted_dic_list, r)
    return pruned_motif_dic_list


def prune_motifs(ts, sorted_dic_list, r):
    """

    :param ts: 1-dimensional time-series either resulting from the PCA method or the original 1-dimensional time-series
    :type ts: 1d array
    :param sorted_dic_list: list of motif dictionaries returned from the emd algorithm, ordered by relevance
    :type sorted_dic_list: list of dic
    :param r: maximum distance to the center of the motif
    :type r: float
    :return: list of dictionaries with the most relevant motifs
    :rtype: list of dic
    """
    pruned_motif_dic_list = [sorted_dic_list[0]]
    first_center_ts = extract_ts_from_pointers(ts, sorted_dic_list[0]['center_ts_pointers'])
    pruned_center_ts_list = [first_center_ts]
    for motif_dic in sorted_dic_list[1:]:
        cur_center_ts = extract_ts_from_pointers(ts, motif_dic['center_ts_pointers'])
        dist_list = dtwdist.compute_dwt_dist_between_ts_and_list(cur_center_ts, pruned_center_ts_list, 2 * r)
        dist_test_list = [dist <= 2 * r for dist in dist_list]
        if sum(dist_test_list) == 0:
            pruned_motif_dic_list.append(motif_dic)
            pruned_center_ts_list.append(cur_center_ts)
        else:
            continue
    return pruned_motif_dic_list


def extract_ts_from_pointers(ts, pointers):
    """

    :param ts: 1-dimensional time-series either resulting from the PCA method or the original 1-dimensional time-series
    :type ts: 1d array
    :param pointers: list of indexes related to the subsequence one wishes to extract from ts
    :type pointers: list of int
    :return: time-series subsequence
    :rtype: 1d array
    """
    ts_from_pointers = np.array([ts[i] for i in pointers])
    return ts_from_pointers
