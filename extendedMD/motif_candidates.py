"""Find all the motif candidates from pre-computed a BS sequence"""
from extendedMD.dtw_dist import compute_dtw_dist_mat
from extendedMD.subseq_patterns import get_bs_subsequences_list
from extendedMD.subseq_patterns import get_all_subsequences_in_pattern
from extendedMD.subseq_patterns import find_index_of_pattern_center_and_members
from extendedMD.mdl import break_bs_len_seq, compute_pattern_mdl


def find_all_motif_candidates(ts, bs_seq, bs_len, bs_point, R):
    # initialize the output lists
    mdl_cost_list = []
    motif_point_list = []
    pattern_list = []
    # initialize the size of the BS subsequences
    subseq_size = 1
    while True:
        subseq_bs_list, subseq_point_list = get_bs_subsequences_list(bs_seq, bs_point, subseq_size)
        subseq_bs_set = [list(item) for item in set(tuple(row) for row in subseq_bs_list)]
        if len(subseq_bs_list)==len(subseq_bs_set):
            break
        for pattern in subseq_bs_set:
            # if the pattern only has one subsequence, then it is not a motif
            if subseq_bs_list.count(pattern) < 2:
                continue
            pattern_ts_list, pattern_pos_list = get_all_subsequences_in_pattern(ts, subseq_bs_list, subseq_point_list, pattern)
            dist_mat = compute_dtw_dist_mat(pattern_ts_list, R)
            # compute the center the the members of the motif candidate and
            # return their position in the BS subsequence
            center_pos, members_pos = find_index_of_pattern_center_and_members(dist_mat, pattern_pos_list, R)
            # if the pattern only has one subsequence, then it is not a motif
            if len(members_pos) < 2:
                continue
            # get the segmentation of the BS subseq fo lenghts, based on the motif (the motif segmentation)
            bs_segmentation_len_list = break_bs_len_seq(bs_len, members_pos, subseq_size)
            # Exclude empty lists (happens when there's overlapping patterns) ------------------------ NEED TO CORRECT THIS!!!!
            bs_segmentation_len_list = [item for item in bs_segmentation_len_list if len(item)>0]
            # compute the MDL of the motif segmentation and append it to the list
            mdl_cost = compute_pattern_mdl(bs_segmentation_len_list)
            mdl_cost_list.append(mdl_cost)
            # extract the motif pointers (i.e. indices of the ts where all the motif members are located)
            # and append it to the list
            motif_point = [subseq_point_list[i] for i in members_pos]
            motif_point_list.append(motif_point)
            # append the pattern to the pattern list
            pattern_list.append(pattern)
        print('motif candidates of size {} successfully extracted'.format(subseq_size))
        # go to the next subsequence size
        subseq_size = subseq_size + 1
    return mdl_cost_list, motif_point_list, pattern_list