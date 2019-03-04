import numpy as np


def extract_modified_bs_sequence(sax_sequence):
    """
    This functions extracts the modified Behaviour Subsequence (BS) sequence list, which is the original sax word
    sequence where every consecutive pairs of sax words are equal are fused into the same sax word.

    :param sax_sequence: list of original sax words
    :type sax_sequence: list of str
    :return:
        - bs_sequence (:py:class:`list of str`) - list of modified sax words
        - bs_lengths (:py:class:`list of int`) - list of lengths of each modified sax word
    """
    bs_sequence = []
    bs_lengths = []
    curr_len = 1
    for i in range(len(sax_sequence)):
        curr_bs = sax_sequence[i]
        if i < len(sax_sequence)-1:
            next_bs = sax_sequence[i+1]
        else:
            next_bs = ''
        if curr_bs == next_bs:
            curr_len = curr_len + 1
        else:
            bs_sequence.append(curr_bs)
            bs_lengths.append(curr_len)
            curr_len = 1
    return bs_sequence, bs_lengths


def generate_bs_pointers(bs_lengths, bs_size):
    """
    It generates the pointers (i.e. time indexes) of each modified sax word into the original time-series data

    :param bs_lengths: list of modified sax words
    :type bs_lengths: list of str
    :param bs_size: window size (in the original time-series) of a single sax word
    :type bs_size: int
    :return: list of pointers to the original time-series
    :rtype: list of list of int
    """
    bs_pointers = []
    start_pointer = 0
    for bs_len_item in bs_lengths:
        end_pointer = start_pointer + bs_size + bs_len_item - 1
        pointer_list = list(range(start_pointer, end_pointer))
        bs_pointers.append(pointer_list)
        start_pointer = start_pointer + bs_len_item
    return bs_pointers


def get_bs_subsequences_dic_list(ts, bs_seq, bs_pointers, subseq_size):
    """
    This function extracts a list with all the BS subsequences with fixed size from a BS sequence

    :param ts: original 1-d
    :type ts: list of float
    :param bs_seq: list of modified sax words (i.e. BS sequence)
    :type bs_seq: list of str
    :param bs_pointers: list of pointers to the original time-series
    :type bs_pointers: list of list of int
    :param subseq_size: number of sax words in a single BS subsequence
    :type subseq_size: int
    :return: list of dictionaries where each dic represents a single BS subsequence. The dic has 4 entries:

        - pattern - the list of sax words related to that BS subsequence
        - pointers - the pointers list to the original time-series related to that subsequence
        - ts - time-series subsequence related to that BS subsequence (numpy array!)
        - bs_position - list of indexes of the the bs_seq list related to that BS subsequence
    :rtype: list of dic
    """
    bs_subseq_dic_list = []
    for i in range(len(bs_seq) - subseq_size + 1):
        subseq_pointers = sorted(list(set(sum(bs_pointers[i:(i + subseq_size)], []))))
        subseq_dic = {
            'pattern': bs_seq[i:(i+subseq_size)],
            'pointers': subseq_pointers,
            'ts': np.array([ts[i] for i in subseq_pointers]),
            'bs_position': list(range(i, (i + subseq_size)))
        }
        bs_subseq_dic_list.append(subseq_dic)
    return bs_subseq_dic_list
