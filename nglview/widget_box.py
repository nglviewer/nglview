from __future__ import absolute_import
from ipywidgets import Box
from .widget import NGLWidget
from .layout import form_item_layout

class NGLHBox(Box):
    def __init__(self, *args, **kwargs):
        super(NGLHBox, self).__init__(*args, **kwargs)

        self.layout = form_item_layout

    def _ipython_display_(self, *args, **kwargs):
        for widget in self.children:
            if isinstance(widget, NGLWidget):
                widget.displayed = True
        super(NGLHBox, self)._ipython_display_(*args, **kwargs)
