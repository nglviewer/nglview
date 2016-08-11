
from __future__ import absolute_import
from ipywidgets import IntSlider, FloatSlider

from nglview import default

def get_widget_by_name(box, widget_name):

    if hasattr(box, '_ngl_children'):
        children = box._ngl_children
    elif hasattr(box, 'children'):
        children = box.children
    else:
        children = None

    if children is not None:
        for widget in children:
            if hasattr(widget, '_ngl_name') and widget._ngl_name == widget_name:
                return widget
    return None

def make_default_slider_width(box):
    for kid in box.children:
        if isinstance(kid, (IntSlider, FloatSlider)):
            kid.width = default.DEFAULT_SLIDER_WIDTH
