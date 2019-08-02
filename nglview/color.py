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


def _singleton(cls):
    # https://www.python.org/dev/peps/pep-0318/#examples
    instances = {}
    def getinstance():
        if cls not in instances:
            instances[cls] = cls()
        print(instances)
        return instances[cls]
    return getinstance


@_singleton
class _ColormakerRegistry(BaseWidget):
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
        if self._ready:
            return str(self)
        super()._ipython_display_(**kwargs)

    def add_selection_scheme(self, scheme_id, arg):
        """
        Examples
        --------
        >>> ColormakerRegistry.add_selection_scheme("my_custom_scheme",
                [['blue', '1-10']])
        >>> view.add_cartoon(color="my_custom_scheme")
        """
        self._call("addSelectionScheme", scheme_id, arg)

    def add_scheme_func(self, scheme_id, func_str):
        """

        Examples
        --------
        >>> func_str = '''
             this.atomColor = function (atom) {
             if (atom.serial < 1000) {
               return 0x0000FF // blue
             } else if (atom.serial > 2000) {
               return 0xFF0000 // red
             } else {
               return 0x00FF00 // green
             }
             }
         '''
        >>> ColormakerRegistry.add_scheme_func('awesome', func_str)
        >>> view.add_cartoon(color='awesome')
        """

        code = """
            var schemeId = NGL.ColormakerRegistry.addScheme(function (params) {
                %s
            })
            
            this._updateId(schemeId, '%s')
        """ % (func_str, scheme_id)
        self._js(code)

    def add_scheme(self, scheme_id, obj):
        """
        Parameters
        ----------
        obj: List of List or str (of JS function)
        """
        if isinstance(obj, list):
            self.add_selection_scheme(scheme_id, obj)
        elif isinstance(obj, str):
            self.add_scheme_func(scheme_id, obj)
        else:
            raise ValueError(f"{obj} must be either list of list or string")

    def _remove_scheme(self, scheme_id):
        self._call("removeScheme", scheme_id)


ColormakerRegistry = _ColormakerRegistry()
