from extendedMD.pca import extract_pca_ts
from extendedMD.sax import extract_sax_sequence
from extendedMD.bs import extract_modified_bs_sequence, generate_bs_pointers
from extendedMD.motifs import find_all_motif_candidates
import pandas as pd


def find_motifs_from_emd(ts, r, win_size, paa_size, alphabet_size, adaptive_break_points=True, z_threshold=0.01):
    """
    Returns the full list of motifs from either a multi-dimensional or a 1-dimensional time-series by
    running the extendedMD algorithm. If the time-series received is multi-dimensional (i.e. a pandas dataframe),
    then the algorithm starts by applying PCA to reduce it to a single dimension.

    :param ts: original time-series
    :type ts: Union[1d array, DataFrame]
    :param r: maximum distance to the center of the motif
    :type r: float
    :param win_size: size fo the sliding window that generated each sax word
    :type win_size: int
    :param paa_size: number of characters in a single sax word
    :type paa_size: int
    :param alphabet_size: number of unique characters to use in the sax representation
    :type alphabet_size: int
    :param adaptive_break_points: Whether to use a representation with adaptive break-points
    :type adaptive_break_points: bool
    :param z_threshold: z_threshold for the znorm method from saxpy
    :type z_threshold: float
    :return:
        - motif_candidates_dic_list (:py:class:`list of dic`) - list of motif dictionaries
        - ts_1d (:py:class:`1d array`) - 1-dimensional time-series either resulting from the PCA method or the original 1-dimensional time-series
    """
    if isinstance(ts, pd.DataFrame):
        ts_1d = extract_pca_ts(ts)
    else:
        ts_1d = ts
    sax_sequence = extract_sax_sequence(ts_1d, win_size, paa_size, alphabet_size, adaptive_break_points, z_threshold)
    bs_sequence, bs_lengths = extract_modified_bs_sequence(sax_sequence)
    bs_pointers = generate_bs_pointers(bs_lengths, win_size)
    motif_candidates_dic_list = find_all_motif_candidates(ts_1d, bs_sequence, bs_lengths, bs_pointers, r)
    return motif_candidates_dic_list, ts_1d
