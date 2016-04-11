from __future__ import absolute_import
from six import string_types

def seq_to_string(seq):
    """e.g. convert [1, 3, 5] to "@1,3,5
    """
    if isinstance(seq, string_types):
        return seq
    else:
        # assume 1D array
        return "@" + ",".join(str(s) for s in seq)
