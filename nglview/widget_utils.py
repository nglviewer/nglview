
from __future__ import absolute_import
from ipywidgets import IntSlider, FloatSlider

from .default import DEFAULT_SLIDER_WIDTH

def get_widget_by_name(box, widget_name):
    for widget in box.children:
        if hasattr(widget, '_ngl_name') and widget._ngl_name == widget_name:
            return widget
    raise ValueError('can not find widget with name = {}'.format(widget_name))

def make_default_slider_width(box):
    for kid in box.children:
        if isinstance(kid, (IntSlider, FloatSlider)):
            kid.width = DEFAULT_SLIDER_WIDTH
