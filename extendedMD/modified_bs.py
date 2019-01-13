"""Generate the modified BS sequence"""

def extract_modified_bs_sequence(sax_sequence):
    # initialize the lists to save the bs and their lenghts
    bs_seq = []
    bs_len = []
    # initialize the bs lenght counter
    curr_len = 1
    for i in range(len(sax_sequence)):
        # set the current bs element
        curr_bs = sax_sequence[i]
        # set the next bs element
        if i<len(sax_sequence)-1:
            next_bs = sax_sequence[i+1]
        else: # if the current element is the last, then thre's no "next_bs"
            next_bs = ''
        # test if the current bs is equal to the next bs
        if curr_bs==next_bs:
            # if yes, add 1 to the current lenght counter
            curr_len = curr_len + 1
        else:
            # if no, save the bs and its lenght in the corresponding lists
            bs_seq.append(curr_bs)
            bs_len.append(curr_len)
            # and initialize the lenght counter
            curr_len = 1
    return bs_seq, bs_len

def generate_bs_pointers(bs_len, bs_size):
    bs_pointers = []
    start_pointer = 0
    for bs_len_item in bs_len:
        end_pointer = start_pointer + bs_size + bs_len_item - 1
        pointer_list = list(range(start_pointer, end_pointer))
        bs_pointers.append(pointer_list)
        start_pointer = start_pointer + bs_len_item
    return bs_pointers