# Code is copied/adapted from nglview
import base64
import ipywidgets as widgets
from traitlets import (Bool, Dict, Integer,
                       Unicode, observe)
from ._frontend import __frontend_version__
from .widget_base import WidgetBase



@widgets.register
class MolstarView(WidgetBase):
    # Name of the widget view class in front-end
    _view_name = Unicode('MolstarView').tag(sync=True)

    # Name of the widget model class in front-end
    _model_name = Unicode('MolstarModel').tag(sync=True)

    # Name of the front-end module containing widget view
    _view_module = Unicode('nglview-js-widgets').tag(sync=True)

    # Name of the front-end module containing widget model
    _model_module = Unicode('nglview-js-widgets').tag(sync=True)

    # Version of the front-end module containing widget view
    _view_module_version = Unicode(__frontend_version__).tag(sync=True)
    # Version of the front-end module containing widget model
    _model_module_version = Unicode(__frontend_version__).tag(sync=True)
    frame = Integer().tag(sync=True)
    loaded = Bool(False).tag(sync=False)
    molstate = Dict().tag(sync=True)
    _molstar_version = Unicode().tag(sync=True)

    def __init__(self):
        super().__init__()
        self._molstar_component_ids = []
        self._state = None

    def _handle_nglview_custom_message(self, widget, msg, buffers):
        msg_type = msg.get("type")
        data = msg.get("data")
        if msg_type == "exportImage":
            image = widgets.Widget.widgets[msg.get("model_id")]
            image.value = base64.b64decode(data)
        elif msg_type == "state":
            self._state = data
        elif msg_type == 'request_loaded':
            if not self.loaded:
                # FIXME: doublecheck this
                # trick to trigger observe loaded
                # so two viewers can have the same representations
                self.loaded = False
            self.loaded = msg.get('data')
        elif msg_type == 'getCamera':
            self._molcamera = data

    @observe('loaded')
    def on_loaded(self, change):
        if change['new']:
            self._fire_callbacks(self._callbacks_before_loaded)

    def _load_structure_data(self, data: str, format: str = 'pdb', preset="default"):
        self._remote_call("loadStructureFromData",
                          target="Widget",
                          args=[data, format, preset])

    def add_trajectory(self, trajectory):
        self._load_structure_data(trajectory.get_structure_string(),
                                  'pdb')  # FIXME
        self._trajlist.append(trajectory)
        self._update_max_frame()
        self._molstar_component_ids.append(trajectory.id)

    def add_structure(self, struc):
        self._load_structure_data(struc.get_structure_string(),
                                  'pdb')
        self._molstar_component_ids.append(struc.id)

    @observe('frame')
    def _on_frame_changed(self, change):
        """set and send coordinates at current frame
        """
        self._set_coordinates(self.frame)
