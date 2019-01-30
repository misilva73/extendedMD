# ---------------------------------------------------------------------------------------------------------------------
# Testing functions
# ---------------------------------------------------------------------------------------------------------------------
import numpy as np
import pandas as pd

def compare_dic_lists(dic_list1, dic_list2):
    if len(dic_list1) == len(dic_list2):
        list_compare = [compare_dics(dic_list1[i], dic_list2[i]) for i in range(len(dic_list1))]
        return sum(list_compare) == 0
    else:
        return False

def compare_dics(dic1, dic2):
    if dic1.keys() == dic2.keys():
        key_compare = [compare_items(dic1[key], dic2[key]) for key in dic1.keys()]
        return sum(key_compare) == 0
    else:
        return False

def compare_items(item1, item2):
    if type(item1) == type(item2):
        if isinstance(item1, np.ndarray):
            return np.array_equal(item1, item2)
        else:
            return item1 == item2

# ---------------------------------------------------------------------------------------------------------------------
# Global variables
# ---------------------------------------------------------------------------------------------------------------------
ts_1d = np.array([1, 1, 2, 2, 4, 4, 7, 6, 4, 4, 2, 2, 1, 1] +
                 [1, 1, 2, 2, 4, 4, 6, 6, 4, 4, 2, 1, 1, 1])
ts_multi_dim = pd.DataFrame({
    'var1': [-1, 0, 1],
    'var2': [-1, 0, 1]
})
pca_ts = np.array([-1.41, 0., 1.41])
sax_seq_adapt = ['ac', 'ac', 'ac', 'ac', 'ac', 'bb', 'ca', 'ca', 'ca', 'ca', 'ca', 'ca', 'cc',
                 'ac', 'ac', 'ac', 'ac', 'ac', 'ac', 'bb', 'ca', 'ca', 'ca', 'ca', 'ca']
sax_seq_fixed = ['aa', 'ab', 'ac', 'bc', 'cc', 'cc', 'cc', 'cb', 'ca', 'ba', 'aa', 'aa', 'aa',
                 'aa', 'aa', 'ab', 'ac', 'bc', 'cc', 'cc', 'cc', 'cb', 'ca', 'ba', 'aa']
bs_seq = ['ac', 'bb', 'ca', 'cc', 'ac', 'bb', 'ca']
bs_len = [5, 1, 6, 1, 6, 1, 5]
bs_pointers = [
    [0, 1, 2, 3, 4, 5, 6, 7],
    [5, 6, 7, 8],
    [6, 7, 8, 9, 10, 11, 12, 13, 14],
    [12, 13, 14, 15],
    [13, 14, 15, 16, 17, 18, 19, 20, 21],
    [19, 20, 21, 22],
    [20, 21, 22, 23, 24, 25, 26, 27]
]
bs_subseq_dic_list = [
    {
        'pattern': ['ac', 'bb'],
        'pointers': [0, 1, 2, 3, 4, 5, 6, 7, 8],
        'ts': np.array([1, 1, 2, 2, 4, 4, 7, 6, 4]),
        'bs_position': [0, 1]
    },
    {
        'pattern': ['bb', 'ca'],
        'pointers': [5, 6, 7, 8, 9, 10, 11, 12, 13, 14],
        'ts': np.array([4, 7, 6, 4, 4, 2, 2, 1, 1, 1]),
        'bs_position': [1, 2]
    },
    {
        'pattern': ['ca', 'cc'],
        'pointers': [6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
        'ts': np.array([7, 6, 4, 4, 2, 2, 1, 1, 1, 1]),
        'bs_position': [2, 3]
    },
    {
        'pattern': ['cc', 'ac'],
        'pointers': [12, 13, 14, 15, 16, 17, 18, 19, 20, 21],
        'ts': np.array([1, 1, 1, 1, 2, 2, 4, 4, 6, 6]),
        'bs_position': [3, 4]
    },
    {
        'pattern': ['ac', 'bb'],
        'pointers': [13, 14, 15, 16, 17, 18, 19, 20, 21, 22],
        'ts': np.array([1, 1, 1, 2, 2, 4, 4, 6, 6, 4]),
        'bs_position': [4, 5]
    },
    {
        'pattern': ['bb', 'ca'],
        'pointers': [19, 20, 21, 22, 23, 24, 25, 26, 27],
        'ts': np.array([4, 6, 6, 4, 4, 2, 1, 1, 1]),
        'bs_position': [5, 6]
    }
]
ts_list = [
    np.array([1, 2, 3, 3]),
    np.array([1, 2, 2, 3, 3, 3]),
    np.array([3, 3, 3, 3])
]
dist_mat = np.array(
    [[0., 0., 2.24],
     [0., 0., 2.45],
     [2.24, 2.45, 0.]]
)
members_dic_list_no_dist = [
    {
        'pattern': ['ca'],
        'pointers': [6, 7, 8, 9, 10, 11, 12, 13, 14],
        'ts': np.array([7, 6, 4, 4, 2, 2, 1, 1, 1]),
        'bs_position': [2]
    }, {
        'pattern': ['ca'],
        'pointers': [20, 21, 22, 23, 24, 25, 26, 27],
        'ts': np.array([6, 6, 4, 4, 2, 2, 1, 1]),
        'bs_position': [6]
    }
]
members_dic_list_with_dist = [
    {
        'pattern': ['ca'],
        'pointers': [6, 7, 8, 9, 10, 11, 12, 13, 14],
        'ts': np.array([7, 6, 4, 4, 2, 2, 1, 1, 1]),
        'bs_position': [2],
        'dist_vec': np.array([0, 1])
    }, {
        'pattern': ['ca'],
        'pointers': [20, 21, 22, 23, 24, 25, 26, 27],
        'ts': np.array([6, 6, 4, 4, 2, 2, 1, 1]),
        'bs_position': [6],
        'dist_vec': np.array([1, 0])
    }
]
split_bs_len_list = [[5, 1], [6], [1, 6, 1], [5]]


# ---------------------------------------------------------------------------------------------------------------------
# pca module tests
# ---------------------------------------------------------------------------------------------------------------------
import extendedMD.pca as pca

test = pca.extract_pca_ts(ts_multi_dim)
if not compare_items(test.round(2), pca_ts):
    print('Error in pca.extract_pca_ts')


# ---------------------------------------------------------------------------------------------------------------------
# sax module tests
# ---------------------------------------------------------------------------------------------------------------------
import extendedMD.sax as sax

test = sax.extract_sax_sequence(ts_1d, 4, 2, 3, True)
if not compare_items(test, sax_seq_adapt):
    print('Error in sax.extract_sax_sequence with adaptive break-points')

test = sax.extract_sax_sequence(ts_1d, 4, 2, 3, False)
if not compare_items(test, sax_seq_fixed):
    print('Error in sax.extract_sax_sequence with fixed break-points')


# ---------------------------------------------------------------------------------------------------------------------
# bs module tests
# ---------------------------------------------------------------------------------------------------------------------
import extendedMD.bs as bs

test_seg, test_len = bs.extract_modified_bs_sequence(sax_seq_adapt)
if not (test_seg == bs_seq) & (test_len == bs_len):
    print('Error in bs.extract_modified_bs_sequence')

test = bs.generate_bs_pointers(bs_len, 4)
if not compare_items(test, bs_pointers):
    print('Error in bs.generate_bs_pointers')

test = bs.get_bs_subsequences_dic_list(ts_1d, bs_seq, bs_pointers, 2)
if not compare_dic_lists(test, bs_subseq_dic_list):
    print('Error in bs.get_bs_subsequences_dic_list')


# ---------------------------------------------------------------------------------------------------------------------
# dtwdist module tests
# ---------------------------------------------------------------------------------------------------------------------
import extendedMD.dtwdist as dtwdist

test = dtwdist.compute_dtw_dist_mat(ts_list)
if not compare_items(test.round(2), dist_mat):
    print('Error in dtwdist.compute_dtw_dist_mat')

test = dtwdist.add_distance_vec_to_pattern_dic(members_dic_list_no_dist, 2)
if not compare_dic_lists(test, members_dic_list_with_dist):
    print('Error in dtwdist.add_distance_vec_to_pattern_dic')


# ---------------------------------------------------------------------------------------------------------------------
# mdl module tests
# ---------------------------------------------------------------------------------------------------------------------
import extendedMD.mdl as mdl

test = mdl.split_bs_len(members_dic_list_with_dist, bs_len)
if not compare_items(test, [[5, 1], [6], [1, 6, 1], [5]]):
    print('Error in mdl.split_bs_len')

test = mdl.compute_segmentation_mdl_cost(split_bs_len_list)
if not compare_items(test, 41.46):
    print('Error in mdl.compute_segmentation_mdl_cost')

test = mdl.compute_motif_mdl_cost(members_dic_list_with_dist, bs_len)
if not compare_items(test, 41.46):
    print('Error in mdl.compute_motif_mdl_cost')


# ---------------------------------------------------------------------------------------------------------------------
# patterns module tests
# ---------------------------------------------------------------------------------------------------------------------
import extendedMD.patterns as patterns

center_test, members_test, mean_dist_test = patterns.find_pattern_center_and_members(members_dic_list_with_dist, 2)
print(center_test)
print(members_test)
print(mean_dist_test)
print('---')


test = patterns.find_all_pruned_members(0, members_dic_list_with_dist, 2)
print(test)


test1 = patterns.lists_overlap([1, 2, 3], [3, 4, 5])
test2 = patterns.lists_overlap([1, 2, 3], [4, 5, 6])
if not compare_items(test1, True) & compare_items(test2, False):
    print('Error in patterns.lists_overlap')


# ---------------------------------------------------------------------------------------------------------------------
# motifs module tests
# ---------------------------------------------------------------------------------------------------------------------
