from saxpy.znorm import znorm
from saxpy.paa import paa
from saxpy.alphabet import cuts_for_asize
from saxpy.sax import ts_to_string

def extract_sax_sequence(ts, win_size, paa_size, alphabet_size=3, z_threshold=0.01):
    """Applies the sax transformation to a 1-dim time series"""
    sax_sequence = []
    cuts = cuts_for_asize(alphabet_size)
    ts_znorm = znorm(ts, z_threshold)
    for t in range(0, len(ts) - win_size + 1):
        ts_win = ts_znorm[t:(t+win_size)]
        paa_rep = paa(ts_win, paa_size)
        sax_word = ts_to_string(paa_rep, cuts)
        sax_sequence.append(sax_word)
    return sax_sequence