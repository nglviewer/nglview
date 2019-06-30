from ipywidgets import Box
from traitlets import Bool, CaselessStrEnum, Unicode, observe

from .layout import make_form_item_layout
from .utils import js_utils
from .widget import NGLWidget


class BoxNGL(Box):
    _gui_style = CaselessStrEnum(['row', 'column'],
                                 default_value='row').tag(sync=True)
    _is_beautified = Bool(False)

    def __init__(self, *args, **kwargs):
        self.layout = make_form_item_layout()
        super().__init__(*args, **kwargs)

    @observe('_gui_style')
    def _update_gui_style(self, change):
        """row or column style
        """
        what = change['new']
        self.layout.flex_flow = what.lower()

    def _ipython_display_(self, *args, **kwargs):
        super()._ipython_display_(*args, **kwargs)
        self._beautify()

    def _update_size(self):
        for widget in self.children:
            if isinstance(widget, NGLWidget):
                widget._remote_call('setSize',
                                    target='Widget',
                                    args=['60%', '60%'])

    def _beautify(self):
        if not self._is_beautified:
            js_utils._set_notebook_width('60%', left_padding=None)
            self._update_size()
            self._is_beautified = True
