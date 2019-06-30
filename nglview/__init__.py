import warnings

# for doc
from . import adaptor, datafiles, show, widget, widget_box
from ._version import get_versions
from .adaptor import *
from .base_adaptor import *
# TODO: do not use import *
# interface
from .config import BACKENDS
from .data_source import DatasourceRegistry
from .show import *
# utils
from .utils import js_utils, widget_utils
from .widget import NGLWidget, write_html

__version__ = get_versions()['version']
del get_versions

with warnings.catch_warnings():
    warnings.simplefilter("ignore")

# Register nbextension


def _jupyter_nbextension_paths():
    return [{
        'section': 'notebook',
        'src': 'static',
        'dest': 'nglview-js-widgets',
        'require': 'nglview-js-widgets/extension'
    }]


__all__ = ['NGLWidget', 'write_html'
           ] + widget.__all__ + adaptor.__all__ + show.__all__
