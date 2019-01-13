import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def create_motif_table(pattern_list, motif_point_list, mdl_cost_list, mean_dist):
    motif_dic = {'pattern' : pattern_list,
                 'pattern_len' : [len(pattern) for pattern in pattern_list],
                 'n_members' : [len(temp_motif) for temp_motif in motif_point_list],
                 'overlap_ratio' : compute_motif_overlap_ratio(motif_point_list),
                 'mdl_cost' : mdl_cost_list,
                 'mean_dist' : mean_dist}
    motif_df = pd.DataFrame(motif_dic) \
                        .assign(mdl_cost = lambda x: x['mdl_cost'].round(1)) \
                        .sort_values('mdl_cost')
    return motif_df


def compute_motif_overlap_ratio(motif_point_list):
    overlap_ratio_list = []
    for motif_point in motif_point_list:
        values, counts = np.unique(np.concatenate(motif_point), return_counts = True)
        overlap_ratio = np.mean(counts>1)
        overlap_ratio_list.append(overlap_ratio)
    return overlap_ratio_list


def plot_single_motif(ts_1d, motif_index, motif_point_list, pattern_list):
    motif_pointers = motif_point_list[motif_index]
    motif_pattern = pattern_list[motif_index]
    # Plots:
    fig = plt.figure(figsize=(12,5))
    plt.suptitle('Pattern {} : {}% overlap'.format(motif_pattern))
    #subplot 1
    plt.subplot(2,1,1)
    plt.plot(ts_1d)
    for temp_point in motif_pointers:
        plt.plot(temp_point, ts_1d[temp_point], 'r')   
    #subplot 2
    plt.subplot(2,1,2)
    for temp_point in motif_pointers:
        plt.plot(ts_1d[temp_point])
    return fig


def plot_k_motifs(k, ts_1d, mdl_cost_list, motif_point_list, pattern_list):
    _mdl_cost_list = mdl_cost_list
    for i in range(k):
        motif_index = np.argmin(_mdl_cost_list)
        fig = plot_single_motif(ts_1d, motif_index, motif_point_list, pattern_list)
        plt.show()
        del _mdl_cost_list[motif_index]
    return 'Plot completed'