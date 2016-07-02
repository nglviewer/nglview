from traitlets import (Unicode, Bool, Dict, List, Int, Float, Any, Bytes, observe,
                       CaselessStrEnum)
from ipywidgets import DOMWidget, interactive

class Representation(DOMWidget):
    parameters = Dict().tag(sync=False)

    def __init__(self, view, cindex, repr_index, *args, **kwargs):
        super(Representation, self).__init__(*args, **kwargs)
        self.cindex = cindex   
        self.repr_index = repr_index
        self._view = view

    @observe('parameters')
    def _on_parameters_changed(self, change):
        parameters = change['new']

        kwargs = dict(component_index=self.cindex,
                      repr_index=self.repr_index)
        kwargs.update(parameters)

        self._view._remote_call('setParameters',
                 target='Representation',
                 kwargs=kwargs)

    def _add_button(self):
        def func(opacity=1.,
                assembly='defaul'):
            self.parameters = dict(opacity=opacity, assembly=assembly)

        assembly_list = ['default', 'AU', 'BU1', 'UNITCELL', 'SUPERCELL']
        return interactive(func, opacity=(0., 1., 0.1),
                                 assembly=assembly_list)
