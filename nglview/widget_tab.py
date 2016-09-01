from __future__ import absolute_import
from ipywidgets import Tab
from traitlets import Unicode

class NGLTab(Tab):
    _view_name = Unicode("NGLTab").tag(sync=True)
    _view_module = Unicode("nglview-js").tag(sync=True)
