"""Runs the full ExtendedMD algorithm"""
from extendedMD.pca import extract_pca_ts
from extendedMD.sax import extract_sax_sequence
from extendedMD.modified_bs import extract_modified_bs_sequence
from extendedMD.modified_bs import generate_bs_pointers
from extendedMD.motif_candidates import find_all_motif_candidates
import pandas as pd


def find_motifs_from_emd(ts, R, win_size, paa_size, alphabet_size, adaptive_break_points=True, z_threshold=0.01):
    """Returns the full list of motifs form a multi-dimensional time-series by running the extendedMD algorithm"""
    if isinstance(ts, pd.DataFrame):
        ts_1d = extract_pca_ts(ts)
    else:
        ts_1d = ts
    sax_sequence = extract_sax_sequence(ts_1d, win_size, paa_size, alphabet_size, adaptive_break_points, z_threshold)
    bs_seq, bs_len = extract_modified_bs_sequence(sax_sequence)
    bs_point = generate_bs_pointers(bs_len, win_size)
    mdl_cost_list, motif_point_list, motif_center_list, pattern_list, mean_dist_list = find_all_motif_candidates(ts_1d, bs_seq, bs_len, bs_point, R)
    return ts_1d, mdl_cost_list, motif_point_list, motif_center_list, pattern_list, mean_dist_list