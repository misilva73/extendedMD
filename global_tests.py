import numpy as np
import pandas as pd
import extendedMD.emd as emd

# Global variables
r = 2
win_size = 4
paa_size = 2
alphabet_size = 3


# Test with a 1-dimensional time-series
ts_1d = np.array([1, 1, 2, 2, 4, 4, 7, 6, 4, 4, 2, 2, 1, 1] +
                 [1, 1, 2, 2, 4, 4, 6, 6, 4, 4, 2, 1, 1, 1])
#motif_list =
test_motif_list, ts_1d = emd.find_motifs_from_emd(ts_1d, r, win_size, paa_size, alphabet_size)
for dic in test_motif_list:
    print(dic)
# assert test_motif_list == motif_list

# Test with a multi-dimensional time-series
ts_multi_dim = pd.DataFrame({
    'var1': ts_1d,
    'var2': ts_1d+2
})
#motif_list =
test_motif_list, ts_1d = emd.find_motifs_from_emd(ts_multi_dim, r, win_size, paa_size, alphabet_size)
# assert test_motif_list == motif_list
