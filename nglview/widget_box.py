from __future__ import absolute_import
from ipywidgets import Box
from .widget import NGLWidget
from .layout import form_item_layout
from . import default

class HBoxNGL(Box):
    def __init__(self, *args, **kwargs):
        super(HBoxNGL, self).__init__(*args, **kwargs)

        self.layout = form_item_layout

    def _ipython_display_(self, *args, **kwargs):
        for widget in self.children:
            if isinstance(widget, NGLWidget):
                widget.displayed = True
        super(HBoxNGL, self)._ipython_display_(*args, **kwargs)

    def _update_padding(self, padding=default.DEFAULT_PADDING):
        for widget in self.children:
            if isinstance(widget, NGLWidget):
                widget.player._create_all_tabs()
                widget.player._update_padding(padding=padding)
