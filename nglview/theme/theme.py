import os
from IPython.display import display
from ipywidgets import HTML, DOMWidget
from traitlets import Unicode, List, observe
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
    _theme_css = Unicode(allow_none=True).tag(sync=True)
    # NOTE: Use `_theme_css` trait to sync backend and frontend to trigger them changed.
    # Use "sync" to make sure theme is created before displaying the NGLWidget.

    # overwrite BaseWidget's trait to avoid caling base method in frontend
    _msg_q = List().tag(sync=True)

    def __init__(self):
        super().__init__()
        self._theme = None
        display(self)

    @observe("_theme_css")
    def _on_theme_changed(self, _):
        """
        Without this method, we can't change the theme with below code

        # cell 1
        >>> import nglview as nv
        >>> 
        >>> view = nv.demo()
        >>> view.gui_style = 'ngl'
        >>> view._widget_theme.light()
        >>> # view._widget_theme.dark()
        >>> view

        # cell 2
        >>> import nglview as nv
        >>> 
        >>> view = nv.demo()
        >>> view.gui_style = 'ngl'
        >>> # view._widget_theme.light()
        >>> view._widget_theme.dark()
        >>> view
        """
        self._call('handleThemeChanged')

    def _ipython_display_(self, **kwargs):
        if self._ready:
            return
        else:
            super()._ipython_display_(**kwargs)

    def remove(self):
        self._js("""
            var ele = document.getElementById('nglview_style')
            if (ele){
                document.head.removeChild(ele)
            }
        """)
        self._theme = None
        self._theme_css = ''

    def dark(self):
        self._theme_css = _get_css_content('dark.css') + _get_css_content('main.css')
        self._theme = 'dark'

    def light(self):
        self._theme_css = _get_css_content('light.css') + _get_css_content('main.css')
        self._theme = 'light'
