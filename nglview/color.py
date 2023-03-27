from ipywidgets import DOMWidget
from traitlets import Unicode, Bool, observe, List
from .base import BaseWidget, _singleton
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


@_singleton
class _ColormakerRegistry(BaseWidget):
    _view_name = Unicode("ColormakerRegistryView").tag(sync=True)
    _view_module = Unicode("nglview-js-widgets").tag(sync=True)
    _view_module_version = Unicode(__frontend_version__).tag(sync=True)
    _model_name = Unicode("ColormakerRegistryModel").tag(sync=True)
    _model_module = Unicode("nglview-js-widgets").tag(sync=True)
    _model_module_version = Unicode(__frontend_version__).tag(sync=True)
    _msg_q = List().tag(sync=True) # overwrite BaseWidget's trait to avoid caling base method in frontend

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        try:
            get_ipython() # only display in notebook
            self._ipython_display_()
        except NameError:
            pass

    def _ipython_display_(self, **kwargs):
        if self._ready:
            return
        try:
            # ipywidgets < 8
            super()._ipython_display_(**kwargs)
        except AttributeError:
            display(super()._repr_mimebundle_(), raw=True)

    def __repr__(self):
        # Prevent ipywidgets to print _ColormakerRegistry() in non-notebook
        # context
        return ""

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


ColormakerRegistry = _ColormakerRegistry()
