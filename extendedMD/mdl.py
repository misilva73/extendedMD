"""Compute segmentation MDL"""
import math
import numpy as np


def break_bs_len_seq(bs_len, break_points, subseq_size):
    bs_len_break = []
    if break_points[0] > 0:
        first_seq = bs_len[0:break_points[0]]
        bs_len_break.append(first_seq)
    for i in range(len(break_points)-1):
        pattern_seq = bs_len[break_points[i]:break_points[i]+subseq_size]
        bs_len_break.append(pattern_seq)
        next_seq = bs_len[break_points[i]+subseq_size:break_points[i+1]]
        bs_len_break.append(next_seq)
    final_pattern_seq = bs_len[break_points[-1]:break_points[-1]+subseq_size]
    bs_len_break.append(final_pattern_seq)
    if break_points[-1]+subseq_size < len(bs_len)-1:
        final_seq = bs_len[break_points[-1]+subseq_size:]
        bs_len_break.append(final_seq)
    return bs_len_break


def compute_pattern_mdl(seg_len_list):
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
    seg_cost = len(seg_len_list)*math.log2(sum(np.concatenate(seg_len_list)))
    mdl_cost = par_cost + data_cost + seg_cost
    return mdl_cost