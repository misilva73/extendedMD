import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


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


def plot_single_motif(ts_1d, events_ts, motif_index, motif_point_list, center_point_list, pattern_list, yaxis_label):
    motif_pointers = motif_point_list[motif_index]
    center_pointers = center_point_list[motif_index]
    motif_pattern = pattern_list[motif_index]
    event_df = pd.DataFrame([ts_1d, events_ts.values]).T.reset_index()
    event_df.columns = ['index', 'var', 'event']
    # Plots:
    fig = plt.figure(figsize=(12,6))
    plt.suptitle(motif_pattern)
    #subplot 1
    #plt.subplot(3,1,1)
    plt.subplot2grid((2, 2), (0, 0), colspan=2)
    plt.plot(ts_1d, 'tab:gray')
    for temp_point in motif_pointers:
        plt.plot(temp_point, ts_1d[temp_point], 'tab:blue')
    sns.scatterplot(x="index", y="var", hue="event", data=event_df[event_df['event']>0], legend=False,
                    palette=sns.xkcd_palette(['red', 'tangerine', 'greenish yellow']))
    plt.ylabel(yaxis_label)
    plt.xlabel('')
    plt.ylim(min(ts_1d), max(ts_1d))
    #subplot 2
    plt.subplot2grid((2, 2), (1, 0))
    for temp_point in motif_pointers:
        plt.plot(ts_1d[temp_point], 'tab:blue')
    plt.ylabel(yaxis_label)
    plt.ylim(min(ts_1d), max(ts_1d))
    #subplot 3
    plt.subplot2grid((2, 2), (1, 1))
    plt.plot(ts_1d[center_pointers], 'tab:blue')
    plt.ylabel(yaxis_label)
    plt.ylim(min(ts_1d), max(ts_1d))
    return fig

def plot_k_motifs(k, ts_1d, events_ts, mdl_cost_list, motif_point_list, center_point_list, pattern_list):
    sorted_index = np.argsort(mdl_cost_list)
    for motif_index in sorted_index[0:k]:
        fig = plot_single_motif(ts_1d, events_ts, motif_index, motif_point_list, center_point_list, pattern_list, 'pca time-series')
        plt.show()
    return 'Plot completed'


def plot_k_multdim_motifs(k, multidim_ts, events_ts, mdl_cost_list, motif_point_list, center_point_list, pattern_list):
    sorted_index = np.argsort(mdl_cost_list)
    for motif_index in sorted_index[0:k]:
        for column in multidim_ts.columns:
            fig1 = plot_single_motif(multidim_ts[column].values, events_ts, motif_index, motif_point_list, center_point_list, pattern_list, column)
            plt.show()