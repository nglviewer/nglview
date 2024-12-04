import threading
import base64
import time
import ipywidgets as widgets
from traitlets import Bool, Dict, Integer, Unicode, observe
from .remote_thread import RemoteCallThread

class WidgetBase(widgets.DOMWidget):
    frame = Integer().tag(sync=True)
    loaded = Bool(False).tag(sync=False)
    _component_ids = []
    _trajlist = []
    _callbacks_before_loaded = []
    _event = threading.Event()

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._initialize_threads()

    def _initialize_threads(self):
        self._remote_call_thread = RemoteCallThread(self, registered_funcs=[])
        self._remote_call_thread.daemon = True
        self._remote_call_thread.start()
        self._handle_msg_thread = threading.Thread(target=self.on_msg, args=(self._handle_nglview_custom_message,))
        self._handle_msg_thread.daemon = True
        self._handle_msg_thread.start()

    def render_image(self):
        image = widgets.Image()
        self._js(f"this.exportImage('{image.model_id}')")
        return image

    def handle_resize(self):
        self._js("this.plugin.handleResize()")

    @observe('loaded')
    def on_loaded(self, change):
        if change['new']:
            self._fire_callbacks(self._callbacks_before_loaded)

    def _thread_run(self, func, *args):
        thread = threading.Thread(target=func, args=args)
        thread.daemon = True
        thread.start()
        return thread

    def _fire_callbacks(self, callbacks):
        def _call(event):
            for callback in callbacks:
                callback(self)
        self._thread_run(_call, self._event)

    def _wait_until_finished(self, timeout=0.0001):
        self._event.clear()
        while True:
            # idle to make room for waiting for
            # "finished" event sent from JS
            time.sleep(timeout)
            if self._event.is_set():
                # if event is set from another thread
                # break while True
                break

    def _js(self, code, **kwargs):
        self._remote_call('executeCode', target='Widget', args=[code], **kwargs)

    def _remote_call(self, method_name, target='Widget', args=None, kwargs=None, **other_kwargs):
        msg = self._get_remote_call_msg(method_name, target=target, args=args, kwargs=kwargs, **other_kwargs)
        def callback(widget, msg=msg):
            widget.send(msg)
        callback._method_name = method_name
        callback._msg = msg
        if self.loaded:
            self._remote_call_thread.q.append(callback)
        else:
            self._callbacks_before_loaded.append(callback)

    def _get_remote_call_msg(self, method_name, target='Widget', args=None, kwargs=None, **other_kwargs):
        msg = {'target': target, 'type': 'call_method', 'methodName': method_name, 'args': args, 'kwargs': kwargs}
        msg.update(other_kwargs)
        return msg

    def _update_max_frame(self):
        self.max_frame = max(int(traj.n_frames) for traj in self._trajlist if hasattr(traj, 'n_frames')) - 1

    def _set_coordinates(self, index):
        if self._trajlist:
            coordinates_dict = {}
            for trajectory in self._trajlist:
                traj_index = self._component_ids.index(trajectory.id)
                try:
                    coordinates_dict[traj_index] = trajectory.get_coordinates(index)
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
        msg = {'type': 'binary_single', 'data': coords_indices}
        self.send(msg, buffers=buffers)

    @observe('frame')
    def _on_frame_changed(self, change):
        self._set_coordinates(self.frame)