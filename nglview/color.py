from ipywidgets import DOMWidget
from traitlets import Unicode, Bool, observe
from .base import BaseWidget
from ._frontend import __frontend_version__
from IPython.display import display
import time


COLOR_SCHEMES = [
    " ", "picking", "random", "uniform", "atomindex", "residueindex",
    "chainindex", "modelindex", "sstruc", "element", "resname", "bfactor",
    "hydrophobicity", "value", "volume", "occupancy"
]

_USER_COLOR_DICT = {}



class _ColorScheme:
    _color_dict = {}

    def __init__(self, args, label):
        # FIXME: validate `args`
        self._color_scheme = args
        self._label = f'user_{label}'
        _USER_COLOR_DICT[self._label] = self._color_scheme

    @property
    def data(self):
        return {'data': self._color_scheme, 'label': self._label}


class ColormakerRegistry(BaseWidget):
    _view_name = Unicode("ColormakerRegistryView").tag(sync=True)
    _view_module = Unicode("nglview-js-widgets").tag(sync=True)
    _view_module_version = Unicode(__frontend_version__).tag(sync=True)
    _model_name = Unicode("ColormakerRegistryModel").tag(sync=True)
    _model_module = Unicode("nglview-js-widgets").tag(sync=True)
    _model_module_version = Unicode(__frontend_version__).tag(sync=True)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        display(self)

    @observe("_ready")
    def _on_ready(self, change):
        if change.new:
            while self._msg_q:
                msg = self._msg_q.pop(0)
                self.send(msg)

    def _ipython_display_(self, **kwargs):
        super()._ipython_display_(**kwargs)

    def add_selection_scheme(self, scheme_id, arg):
        self._call("addSelectionScheme", scheme_id, arg)

    def add_scheme(self, func_str):
        self._call("addSchemeByFunc", func_str)
