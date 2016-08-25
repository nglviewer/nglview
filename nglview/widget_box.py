from __future__ import absolute_import
from ipywidgets import HBox
from .widget import NGLWidget

class NGLHBox(HBox):
    def __init__(self, *args, **kwargs):
        super(NGLHBox, self).__init__(*args, **kwargs)

    def _ipython_display_(self, *args, **kwargs):
        for widget in self.chilren:
            if isinstance(widget, NGLWidget):
                widget.displayed = True
        super(NGLHBox, self)._ipython_display_(*args, **kwargs)
