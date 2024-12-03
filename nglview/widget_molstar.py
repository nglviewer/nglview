# Code is copied/adapted from nglview
import threading
import base64
import ipywidgets as widgets
from traitlets import (Bool, Dict, Integer,
                       Unicode, observe)
from ._frontend import __frontend_version__

from .remote_thread import RemoteCallThread


@widgets.register
class MolstarView(widgets.DOMWidget):
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
        self._remote_call_thread = RemoteCallThread(
            self,
            registered_funcs=[])
        self._remote_call_thread.daemon = True
        self._remote_call_thread.start()
        self._handle_msg_thread = threading.Thread(
            target=self.on_msg, args=(self._molstar_handle_message, ))
        # register to get data from JS side
        self._handle_msg_thread.daemon = True
        self._handle_msg_thread.start()
        self._state = None

    def render_image(self):
        image = widgets.Image()
        self._js(f"this.exportImage('{image.model_id}')")
        # image.value will be updated in _molstar_handle_message
        return image

    def handle_resize(self):
        self._js("this.plugin.handleResize()")

    @observe('loaded')
    def on_loaded(self, change):
        if change['new']:
            self._fire_callbacks(self._callbacks_before_loaded)

    def _thread_run(self, func, *args):
        thread = threading.Thread(
            target=func,
            args=args,
        )
        thread.daemon = True
        thread.start()
        return thread

    def _fire_callbacks(self, callbacks):
        def _call(event):
            for callback in callbacks:
                callback(self)
        self._thread_run(_call, self._event)

    def _wait_until_finished(self, timeout=0.0001):
        # FIXME: dummy for now
        pass

    def _load_structure_data(self, data: str, format: str = 'pdb', preset="default"):
        self._remote_call("loadStructureFromData",
                          target="Widget",
                          args=[data, format, preset])

    def _molstar_handle_message(self, widget, msg, buffers):
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

    def render_image(self):
        image = widgets.Image()
        self._js(f"this.exportImage('{image.model_id}')")
        # image.value will be updated in _molview_handle_message
        return image

    def _js(self, code, **kwargs):
        # nglview code
        self._remote_call('executeCode',
                          target='Widget',
                          args=[code],
                          **kwargs)

    def _remote_call(self,
                     method_name,
                     target='Widget',
                     args=None,
                     kwargs=None,
                     **other_kwargs):

        # adapted from nglview
        msg = self._get_remote_call_msg(method_name,
                                        target=target,
                                        args=args,
                                        kwargs=kwargs,
                                        **other_kwargs)
        def callback(widget, msg=msg):
            widget.send(msg)

        callback._method_name = method_name
        callback._msg = msg

        if self.loaded:
            self._remote_call_thread.q.append(callback)
        else:
            # send later
            # all callbacks will be called right after widget is loaded
            self._callbacks_before_loaded.append(callback)

    def _get_remote_call_msg(self,
                             method_name,
                             target='Widget',
                             args=None,
                             kwargs=None,
                             **other_kwargs):
        # adapted from nglview
        msg = {}
        msg['target'] = target
        msg['type'] = 'call_method'
        msg['methodName'] = method_name
        msg['args'] = args
        msg['kwargs'] = kwargs
        return msg

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
