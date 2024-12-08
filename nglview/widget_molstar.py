# Code is copied/adapted from nglview
import base64
import ipywidgets as widgets
from .utils import widget_utils
from traitlets import (Bool, Dict, Integer,
                       Unicode)
from ._frontend import __frontend_version__
from .widget_base import WidgetBase
from .utils.py_utils import (_camelize_dict)



@widgets.register
class MolstarView(WidgetBase):
    _view_name = Unicode('MolstarView').tag(sync=True)
    _model_name = Unicode('MolstarModel').tag(sync=True)

    frame = Integer().tag(sync=True)
    loaded = Bool(False).tag(sync=False)
    molstate = Dict().tag(sync=True)
    _molstar_version = Unicode().tag(sync=True)

    def __init__(self):
        super().__init__()
        self._molstar_component_ids = []
        self._state = None
        widget_utils._add_repr_method_shortcut(self, self)

    def _handle_nglview_custom_message(self, widget, msg, buffers):
        msg_type = msg.get("type")
        data = msg.get("data")

        if msg_type == "exportImage":
            self._handle_export_image(msg)
        elif msg_type == "state":
            self._handle_state(data)
        elif msg_type == 'request_loaded':
            self._handle_request_loaded(msg)
        elif msg_type == 'getCamera':
            self._handle_get_camera(data)

    def _handle_export_image(self, msg):
        image = widgets.Widget.widgets[msg.get("model_id")]
        image.value = base64.b64decode(msg.get("data"))

    def _handle_state(self, data):
        self._state = data

    def _handle_request_loaded(self, msg):
        if not self.loaded:
            # FIXME: doublecheck this
            # trick to trigger observe loaded
            # so two viewers can have the same representations
            self.loaded = False
        self.loaded = msg.get('data')

    def _handle_get_camera(self, data):
        self._molcamera = data

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

    def add_component(self, component):
        raise NotImplementedError()

    def add_representation(self, **params):
        params = _camelize_dict(params)

        if 'component' in params:
            model_index = params.pop('component')
        else:
            model_index = 0

        for k, v in params.items():
            try:
                params[k] = v.strip()
            except AttributeError:
                # e.g.: opacity=0.4
                params[k] = v

        self._remote_call('addRepresentation',
                          args=[
                              params, model_index
                          ])

    def load_spec(self, state, **options):
        self._remote_call('loadMolstarSpec',
                  args=[state, options])