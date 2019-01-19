import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def create_motif_table(pattern_list, motif_point_list, mdl_cost_list, mean_dist):
    motif_dic = {'pattern' : pattern_list,
                 'pattern_len' : [len(pattern) for pattern in pattern_list],
                 'n_members' : [len(temp_motif) for temp_motif in motif_point_list],
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


def plot_single_motif(ts_1d, motif_index, motif_point_list, center_point_list, pattern_list, yaxis_label):
    motif_pointers = motif_point_list[motif_index]
    center_pointers = center_point_list[motif_index]
    motif_pattern = pattern_list[motif_index]
    # Plots:
    fig = plt.figure(figsize=(12,6))
    plt.suptitle(motif_pattern)
    #subplot 1
    plt.subplot(3,1,1)
    plt.plot(ts_1d)
    for temp_point in motif_pointers:
        plt.plot(temp_point, ts_1d[temp_point], 'r')
    plt.ylabel(yaxis_label)
    #subplot 2
    plt.subplot(3,1,2)
    for temp_point in motif_pointers:
        plt.plot(ts_1d[temp_point])
    plt.ylabel(yaxis_label)
    #subplot 3
    plt.subplot(3,1,3)
    plt.plot(ts_1d[center_pointers])
    plt.ylabel(yaxis_label)
    return fig


def plot_k_motifs(k, ts_1d, mdl_cost_list, motif_point_list, center_point_list, pattern_list):
    sorted_index = np.argsort(mdl_cost_list)
    for motif_index in sorted_index[0:k]:
        fig = plot_single_motif(ts_1d, motif_index, motif_point_list, center_point_list, pattern_list, 'pca time-series')
        plt.show()
    return 'Plot completed'


def plot_k_multdim_motifs(k, multidim_ts, mdl_cost_list, motif_point_list, center_point_list, pattern_list):
    sorted_index = np.argsort(mdl_cost_list)
    for motif_index in sorted_index[0:k]:
        fig1 = plot_single_motif(multidim_ts.iloc[:,0].values, motif_index, motif_point_list,
                                 center_point_list, pattern_list, multidim_ts.columns[0])
        plt.show()
        fig2 = plot_single_motif(multidim_ts.iloc[:,1].values, motif_index, motif_point_list,
                                 center_point_list, pattern_list, multidim_ts.columns[1])
        plt.show()