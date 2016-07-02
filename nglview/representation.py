from __future__ import absolute_import
from traitlets import (Unicode, Bool, Dict, List, Int, Float, Any, Bytes, observe,
                       CaselessStrEnum)
from ipywidgets import DOMWidget, interactive

# local
from .colors import color_scheme as COLOR_SCHEMES

class Representation(DOMWidget):
    parameters = Dict().tag(sync=False)

    def __init__(self, view, component_index, repr_index, *args, **kwargs):
        super(Representation, self).__init__(*args, **kwargs)
        self.component_index = component_index   
        self.repr_index = repr_index
        self._view = view

    @observe('parameters')
    def _on_parameters_changed(self, change):
        parameters = change['new']

        kwargs = dict(component_index=self.component_index,
                      repr_index=self.repr_index)
        kwargs.update(parameters)

        self._view._remote_call('setParameters',
                 target='Representation',
                 kwargs=kwargs)

    def _add_button(self):
        def func(opacity=1.,
                assembly='default',
                color_scheme=""):
            self.parameters = dict(opacity=opacity, assembly=assembly,
                    colorScheme=color_scheme)

        assembly_list = ['default', 'AU', 'BU1', 'UNITCELL', 'SUPERCELL']
        return interactive(func, opacity=(0., 1., 0.1),
                                 color_scheme=COLOR_SCHEMES,
                                 assembly=assembly_list)
