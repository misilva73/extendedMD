import numpy as np
import math
import extendedMD.dtwdist as dtwdist


def prune_motifs_with_mdl(ts, motif_dic_list, r):
    sorted_dic_list = sorted(motif_dic_list, key=lambda dic: dic['mdl_cost'])
    pruned_motif_dic_list = prune_motifs(ts, sorted_dic_list, r)
    return pruned_motif_dic_list


def prune_motifs_with_dist(ts, motif_dic_list, r, mdl_splits):
    mdl_sorted_dic_list = sorted(motif_dic_list, key=lambda dic: dic['mdl_cost'])
    step = math.floor(len(mdl_sorted_dic_list)/mdl_splits)
    dist_sorted_dic_list = []
    for i in range(mdl_splits):
        temp_dic_list = mdl_sorted_dic_list[i * step:(i + 1) * step]
        temp_dist_sorted_dic_list = sorted(temp_dic_list, key=lambda dic: dic['mean_dist'])
        dist_sorted_dic_list += temp_dist_sorted_dic_list
    if mdl_splits * step < len(mdl_sorted_dic_list):
        temp_dic_list = mdl_sorted_dic_list[mdl_splits * step:]
        temp_dist_sorted_dic_list = sorted(temp_dic_list, key=lambda dic: dic['mean_dist'])
        dist_sorted_dic_list += temp_dist_sorted_dic_list
    pruned_motif_dic_list = prune_motifs(ts, dist_sorted_dic_list, r)
    return pruned_motif_dic_list


def prune_motifs(ts, sorted_dic_list, r):
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
    ts_from_pointers = np.array([ts[i] for i in pointers])
    return ts_from_pointers
