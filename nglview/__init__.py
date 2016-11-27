from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

from .install import install, enable_nglview_js
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")

# we already install from setup
# but it's better to install again. haizz
install()

# Register nbextension


def _jupyter_nbextension_paths():
    return [{
        'section': 'notebook',
        'src': 'static',
        'dest': 'nglview',
        'require': 'nglview/extension'
    }]

def _jupyter_labextension_paths():
    return [{
        'name': 'nglview',
        'src': 'staticlab',
    }]

enable_nglview_js()

# TODO: do not use import *
# interface
from .config import BACKENDS
from .widget import NGLWidget
from .base_adaptor import *
from .adaptor import *
from .show import *
from . import datafiles

# utils
from .utils import widget_utils, js_utils

# for doc
from . import widget_box, widget, adaptor, show

__all__ = ['NGLWidget'] + widget.__all__ + adaptor.__all__ + show.__all__
