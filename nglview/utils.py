from __future__ import absolute_import
import sys

PY2 = sys.version_info[0] == 2
PY3 = sys.version_info[0] == 3

if PY3:
    string_types = str,
else:
    string_types = basestring,


def seq_to_string(seq):
    """e.g. convert [1, 3, 5] to "@1,3,5
    """
    if isinstance(seq, string_types):
        return seq
    else:
        # assume 1D array
        return "@" + ",".join(str(s) for s in seq)

def _camelize(snake):
    """
    
    Examples
    --------
    >>> _camelize('remote_call')
    remoteCall
    >>> _camelize('remoteCall')
    remoteCall
    """
    words = snake.split('_')
    return words[0] + "".join(x.title() for x in words[1:])
