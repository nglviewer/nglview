import sys
import unittest
import doctest
import nglview
from nglview import widget, show, viewer_control

import warnings
warnings.filterwarnings('ignore')

PY3 = sys.version_info[0] == 3

try:
    from make_dummy_comm import *
except ImportError:
    has_nglview = False
    nglview = Comm = DummyComm = _widget_attrs = displayed = undefined = Widget = None


def get_total_errors(modules):
    return sum([doctest.testmod(mod).failed for mod in modules])


@unittest.skipUnless(PY3, 'doctest with py3 only')
def test_nglview_show_module():
    """
    """
    assert not get_total_errors([show])


@unittest.skipUnless(PY3, 'doctest with py3 only')
def test_nglview_viewer_control():
    """
    """
    assert not get_total_errors([viewer_control])


@unittest.skipUnless(PY3, 'doctest with py3 only')
def test_nglview_widget():
    """
    """
    assert not get_total_errors([widget])
