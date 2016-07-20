from __future__ import absolute_import
from traitlets import Dict, observe
from ipywidgets import DOMWidget, interactive
from ipywidgets import VBox

# local
from .colors import color_schemes as COLOR_SCHEMES
from .widget_utils import make_default_slider_width

class Representation(DOMWidget):
    parameters = Dict().tag(sync=False)

    def __init__(self, view, component_index, repr_index, name=None, *args, **kwargs):
        super(Representation, self).__init__(*args, **kwargs)
        self.component_index = component_index   
        self.repr_index = repr_index
        self._view = view
        self.name = name

    @observe('parameters')
    def _on_parameters_changed(self, change):
        parameters = change['new']

        self._view.update_representation(component=self.component_index,
                repr_index=self.repr_index,
                **parameters)

    def _display(self):
        c_string = 'c' + str(self.component_index)
        r_string = str(self.repr_index)
        try:
            _repr_dict = self._view._repr_dict[c_string][r_string]['parameters']
        except KeyError:
            _repr_dict = dict()

        def func(opacity=_repr_dict.get('opacity', 1.),
                 assembly=_repr_dict.get('assembly', 'default'),
                 color_scheme=_repr_dict.get('colorScheme', ""),
                 wireframe=_repr_dict.get('wireframe', False)):
            parameters = dict(opacity=opacity,
                    assembly=assembly,
                    colorScheme=color_scheme,
                    wireframe=wireframe)
            if not color_scheme:
                parameters.pop('colorScheme')

            self.parameters = parameters 

        assembly_list = ['default', 'AU', 'BU1', 'UNITCELL', 'SUPERCELL']
        iwidget = interactive(func, opacity=(0., 1., 0.1),
                                 color_scheme=COLOR_SCHEMES,
                                 assembly=assembly_list)
        make_default_slider_width(iwidget)
        wbox = VBox([iwidget,])
        if self.name == 'surface':
            def func_extra(probe_radius=1.4,
                    isolevel=2.,
                    smooth=2.,
                    surface_type='ms',
                    box_size=10,
                    cutoff=0.):
                self.parameters = dict(probeRadius=probe_radius,
                        isolevel=isolevel,
                        smooth=smooth,
                        surfaceType=surface_type,
                        boxSize=box_size,
                        cutoff=cutoff)
            surface_types = ['vws', 'sas', 'ms', 'ses']
            # use continuous_update=False to avoid expensive surface calculation and update
            widget_extra = interactive(func_extra,
                    probe_radius=(0., 5., 0.1),
                    isolevel=(0., 10., 0.1),
                    smooth=(0, 10, 1),
                    surface_type=surface_types,
                    box_size=(0, 100, 2),
                    cutoff=(0., 100, 0.1),
                    continuous_update=False)

            make_default_slider_width(widget_extra)
            wbox.children = [iwidget, widget_extra]
        return wbox
