import warnings
import sys

from . import adaptor, datafiles, show, widget
from .adaptor import *
from .base_adaptor import *
from .config import BACKENDS
from .data_source import DatasourceRegistry
from .show import *
from .utils import js_utils, widget_utils
from .widget import NGLWidget, write_html

if sys.version_info >= (3, 8):
    import importlib.metadata
    try:
        __version__ = importlib.metadata.version("nglview")
    except importlib.metadata.PackageNotFoundError:
        __version__ = "unknown"
    del importlib.metadata
else:
    # pkg_resources is deprecated, only use it if importlib.metadata is not available (Python < 3.8)
    import pkg_resources
    try:
        __version__ = pkg_resources.get_distribution("nglview").version
    except pkg_resources.DistributionNotFound:
        __version__ = "unknown"
    del pkg_resources

with warnings.catch_warnings():
    warnings.simplefilter("ignore")


# Register nbextension
# FIXME: do we still need this?
def _jupyter_nbextension_paths():
    return [{
        'section': 'notebook',
        'src': 'static',
        'dest': 'nglview-js-widgets',
        'require': 'nglview-js-widgets/extension'
    }]


__all__ = ['NGLWidget', 'write_html'
          ] + widget.__all__ + adaptor.__all__ + show.__all__
