import os
from IPython.display import display
from ipywidgets import HTML, DOMWidget
from traitlets import Unicode, List
from pathlib import Path

from ..base import _singleton, BaseWidget
from .._frontend import __frontend_version__


def _get_css_content(css_file):
    p = Path(__file__).resolve().parent / css_file
    return p.read_text()


@_singleton
class ThemeManager(BaseWidget):
    """EXPERIMENT
    """
    _view_name = Unicode("ThemeManagerView").tag(sync=True)
    _view_module = Unicode("nglview-js-widgets").tag(sync=True)
    _view_module_version = Unicode(__frontend_version__).tag(sync=True)
    _model_name = Unicode("ThemeManagerModel").tag(sync=True)
    _model_module = Unicode("nglview-js-widgets").tag(sync=True)
    _model_module_version = Unicode(__frontend_version__).tag(sync=True)

    _msg_q = List().tag(sync=True) # overwrite BaseWidget's trait to avoid caling base method in frontend

    def __init__(self):
        super().__init__()
        display(self)

    def _ipython_display_(self, **kwargs):
        if self._ready:
            return
        else:
            super()._ipython_display_(**kwargs)

    def remove(self):
        self._js("""
            var ele = document.getElementById('nglview_style')
            document.head.removeChild(ele)
        """)

    def dark(self):
        self._call(
            "setTheme",
            [_get_css_content('dark.css') + _get_css_content('main.css')])

    def light(self):
        self._call(
            "setTheme",
            [_get_css_content('light.css') + _get_css_content('main.css')])
