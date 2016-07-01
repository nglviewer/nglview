from ipywidgets import DOMWidget

class Representation(DOMWidget):
    params = Dict().tag(sync=False)

    def __init__(self, view, *args, **kwargs):
        super(Representation, self).__init__(*args, **kwargs)
        self._coponent_index = self.type = self.params = None
        self._view = view

    @property
    def component_index(self):
        pass

    @property
    def type(self):
        pass

    @observe('params'):
    def params(self, change):
        params = change['new']
        self._view._update_representation(component=self.component_index,
                params=params)

    def _add_button(self):
        def func(opacity=1.):
            params = dict(opacity=opacity)
            self._view._update_representation(component=self.component_index,
                    params=params)
        return interative(func, opacity=(0., 1., 0.1))
