import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def create_motif_table(motif_dic_list):
    """
    This function creates a pandas dataframe with a summary of each motif in motif_dic_list

    :param motif_dic_list: list of dictionaries, where a dic is related to a single motif
    :type motif_dic_list: list of dic
    :return: motif_table_df
    :rtype: pandas.DataFrame
    """
    pattern_list = [dic['pattern'] for dic in motif_dic_list]
    mdl_cost_list = [dic['mdl_cost'] for dic in motif_dic_list]
    mean_dist_list = [dic['mean_dist'] for dic in motif_dic_list]
    n_members_list = [len(dic['members_ts_pointers']) for dic in motif_dic_list]
    motif_table_dic = {
        'pattern': pattern_list,
        'pattern_len': [len(pattern) for pattern in pattern_list],
        'n_members': n_members_list,
        'mdl_cost': mdl_cost_list,
        'mean_dist': mean_dist_list
    }
    motif_table_df = pd.DataFrame(motif_table_dic)
    return motif_table_df


def plot_single_motif(ts, events_ts, motif_dic, y_label):
    """
    This function creates the base visualization for a single motif:
     1) plot with the whole time-series highlighting the labels and the position of each motif's member
     2) plot with all the motif's members
     3) plot with the motif's center

    :param ts: original 1-dimensional time-series
    :type ts: 1d array
    :param events_ts: list of labels for each entry in ts (max of 10 labels for avoiding color cycling in the plots)
                      label == 0 is assumed to represent no events and thus it is not plotted
    :type events_ts: list
    :param motif_dic: dictionary related to the motif
    :type motif_dic: dic
    :param y_label: label to add in the y-axis of the first plot
    :type y_label: str
    :return: fig - figure with the motif plot
    :rtype: matplotlib.figure.Figure
    """
    member_pointers = motif_dic['members_ts_pointers']
    center_pointers = motif_dic['center_ts_pointers']
    raw_event_df = pd.DataFrame([ts, events_ts]).T.reset_index()
    raw_event_df.columns = ['index', 'var', 'event']
    event_df = raw_event_df[raw_event_df['event'] > 0]
    n_events = len(set(event_df['event'].values))
    # Plots:
    fig = plt.figure(figsize=(12, 6))
    plt.suptitle(motif_dic['pattern'])
    # subplot 1
    plt.subplot2grid((2, 2), (0, 0), colspan=2)
    plt.plot(ts, 'xkcd:grey', alpha=0.5)
    for temp_point in member_pointers:
        plt.plot(temp_point, ts[temp_point], 'xkcd:dark grey')
    sns.scatterplot(x="index", y="var", hue="event", data=event_df, legend=False,
                    palette=sns.color_palette(n_colors=n_events))
    plt.ylabel(y_label)
    plt.xlabel('')
    plt.ylim(min(ts), max(ts))
    # subplot 2
    plt.subplot2grid((2, 2), (1, 0))
    for temp_point in member_pointers:
        plt.plot(ts[temp_point], 'xkcd:dark grey')
    plt.ylabel("Motif's members")
    plt.ylim(min(ts), max(ts))
    # subplot 3
    plt.subplot2grid((2, 2), (1, 1))
    plt.plot(ts[center_pointers], 'xkcd:dark grey')
    plt.ylabel("Motif's center")
    plt.ylim(min(ts), max(ts))
    return fig


def plot_k_motifs(k, ts, events_ts, motif_dic_list, y_label='Full 1-d time-series'):
    """
    This function prints the base visualisation for the first k motifs in motif_dic_list for the original 1-d time-series

    :param k: number of motifs to plot
    :type k: int
    :param ts: original 1-dimensional time-series
    :type ts: 1d array
    :param events_ts: list of labels for each entry in ts
    :type events_ts: list
    :param motif_dic_list: list of dictionaries, where a dic is related to a single motif
    :type motif_dic_list: list of dic
    :param y_label: label to add in the y-axis of the first plot (optional and defaults to 'Full 1-d time-series')
    :type y_label: str
    :return: None
    """
    for motif_dic in motif_dic_list[0:k]:
        plot_single_motif(ts, events_ts, motif_dic, y_label)
        plt.show()


def plot_k_multdim_motifs(k, multidim_ts, events_ts, motif_dic_list):
    """
    This function prints the base visualisation for the first k motifs in motif_dic_list for the original multidimensional
    time-series. It shows one 1-d plot for each dimension in multidim_ts

    :param k: number of motifs to plot
    :type k: int
    :param multidim_ts: original multidimensional time-series
    :type multidim_ts: pandas.DataFrame
    :param events_ts: list of labels for each entry in multidim_ts
    :type events_ts: list
    :param motif_dic_list: list of dictionaries, where a dic is related to a single motif
    :type motif_dic_list: list of dic
    :return: None
    """
    for motif_dic in motif_dic_list[0:k]:
        for column in multidim_ts.columns:
            plot_single_motif(multidim_ts[column].values, events_ts, motif_dic, column)
            plt.show()
