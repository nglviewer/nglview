import doctest
import warnings

from nglview import show, viewer_control, widget

warnings.filterwarnings('ignore')


def get_total_errors(modules):
    return sum([doctest.testmod(mod).failed for mod in modules])


def test_nglview_show_module():
    """
    """
    assert not get_total_errors([show])


def test_nglview_viewer_control():
    """
    """
    assert not get_total_errors([viewer_control])


def test_nglview_widget():
    """
    """
    assert not get_total_errors([widget])
