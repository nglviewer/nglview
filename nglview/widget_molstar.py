# Code is copied/adapted from nglview
import threading
import base64
import ipywidgets as widgets
from traitlets import (Bool, Dict, Integer,
                       Unicode, observe)
from ._frontend import __frontend_version__
from .widget import BaseWidget

from .remote_thread import RemoteCallThread


@widgets.register
class MolstarView(BaseWidget):
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

    def __init__(self):
        super().__init__()
        self._molstar_component_ids = []
        self._trajlist = []
        self._callbacks_before_loaded = []
        self._event = threading.Event()

    def _handle_custom_widget_msg(self, widget, msg, buffers):
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

    def _update_max_frame(self):
        self.max_frame = max(
            int(traj.n_frames) for traj in self._trajlist
            if hasattr(traj, 'n_frames')) - 1 # index starts from 0

    def _set_coordinates(self, index):
        '''update coordinates for all trajectories at index-th frame
        '''
        if self._trajlist:
            coordinates_dict = {}
            for trajectory in self._trajlist:
                traj_index = self._molstar_component_ids.index(trajectory.id)
                try:
                    coordinates_dict[traj_index] = trajectory.get_coordinates(
                        index)
                except (IndexError, ValueError):
                    coordinates_dict[traj_index] = np.empty((0), dtype='f4')
            self._send_coordinates(coordinates_dict)

    def _send_coordinates(self, arr_dict):
        self._coordinates_dict = arr_dict

        buffers = []
        coords_indices = dict()
        for index, arr in self._coordinates_dict.items():
            buffers.append(arr.astype('f4').tobytes())
            coords_indices[index] = index
        msg = {
            'type': 'binary_single',
            'data': coords_indices,
        }
        self.send(msg, buffers=buffers)

    @observe('frame')
    def _on_frame_changed(self, change):
        """set and send coordinates at current frame
        """
        self._set_coordinates(self.frame)
