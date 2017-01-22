import sys
import unittest
import doctest
import nglview
from nglview import widget, show

import warnings
warnings.filterwarnings('ignore')

PY3 = sys.version_info[0] == 3


try:
    from ipywidgets import Widget
    from ipykernel.comm import Comm
    import nglview
    has_nglview = True

    #------------------------------------------------------
    # Utility stuff from ipywidgets tests: create DummyComm
    # we dont need Jupyter notebook for testing
    #------------------------------------------------------
    class DummyComm(Comm):
        comm_id = 'a-b-c-d'
    
        def open(self, *args, **kwargs):
            pass
    
        def send(self, *args, **kwargs):
            pass
    
        def close(self, *args, **kwargs):
            pass
    
    _widget_attrs = {}
    displayed = []
    undefined = object()

    _widget_attrs['_comm_default'] = getattr(Widget, '_comm_default', undefined)
    Widget._comm_default = lambda self: DummyComm()
    _widget_attrs['_ipython_display_'] = Widget._ipython_display_
    def raise_not_implemented(*args, **kwargs):
        raise NotImplementedError()
    Widget._ipython_display_ = raise_not_implemented
except ImportError:
    has_nglview = False
    nglview = Comm = DummyComm = _widget_attrs = displayed = undefined = Widget  = None

def get_total_errors(modules):
    return sum([doctest.testmod(mod).failed for mod in modules])

@unittest.skipUnless(PY3, 'doctest with py3 only')
def test_nglview_show_module():
    """
    """
    assert not get_total_errors([show,])
