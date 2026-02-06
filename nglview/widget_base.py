import threading
import time
import ipywidgets as widgets
from traitlets import Bool, Integer, observe, Unicode
from .remote_thread import RemoteCallThread
from IPython.display import display
import numpy as np

from ._frontend import __frontend_version__

class WidgetBase(widgets.DOMWidget):
    _view_module = Unicode('nglview-js-widgets').tag(sync=True)
    _model_module = Unicode('nglview-js-widgets').tag(sync=True)
    _view_module_version = Unicode(__frontend_version__).tag(sync=True)
    _model_module_version = Unicode(__frontend_version__).tag(sync=True)

    frame = Integer().tag(sync=True)
    loaded = Bool(False).tag(sync=False)

    def __init__(self, **kwargs):
        # Extract recognized arguments
        recognized_kwargs = {k: v for k, v in kwargs.items() if k in self.trait_names()}
        super().__init__(**recognized_kwargs)
        self._view_component_ids = []
        self._callbacks_before_loaded = []
        self._event = threading.Event()
        self._trajlist = []
        self._initialize_threads()

    def _initialize_threads(self):
        self._remote_call_thread = RemoteCallThread(self, registered_funcs=[])
        self._remote_call_thread.daemon = True
        self._remote_call_thread.start()
        self._handle_msg_thread = threading.Thread(target=self.on_msg, args=(self._handle_nglview_custom_message,))
        self._handle_msg_thread.daemon = True
        self._handle_msg_thread.start()

    def _handle_nglview_custom_message(self, widget, msg, buffers):
        raise NotImplementedError()

    def render_image(self):
        image = widgets.Image()
        self._js(f"this.exportImage('{image.model_id}')")
        return image

    def handle_resize(self):
        self._js("this.plugin.handleResize()")

    @observe('loaded')
    def on_loaded(self, change):
        # trick for firefox on Linux
        time.sleep(0.1)
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

    def _update_max_frame(self):
        self.max_frame = max(
            int(traj.n_frames) for traj in self._trajlist
            if hasattr(traj, 'n_frames')) - 1 # index starts from 0

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

    def _set_coordinates(self, index, movie_making=False, render_params=None):
        '''update coordinates for all trajectories at index-th frame'''
        render_params = render_params or {}
        if self._trajlist:
            coordinates_dict = {}
            for trajectory in self._trajlist:
                traj_index = self._view_component_ids.index(trajectory.id)

                try:
                    if trajectory.shown:
                        coordinates_dict[traj_index] = trajectory.get_coordinates(index)
                    else:
                        coordinates_dict[traj_index] = np.empty((0), dtype='f4')
                except (IndexError, ValueError):
                    coordinates_dict[traj_index] = np.empty((0), dtype='f4')

            self.set_coordinates(coordinates_dict,
                                render_params=render_params,
                                movie_making=movie_making)
        else:
            print("no trajectory available")

    def set_coordinates(self, arr_dict, movie_making=False, render_params=None):
        """Used for update coordinates of a given trajectory
        >>> # arr: numpy array, ndim=2
        >>> # update coordinates of 1st trajectory
        >>> view.set_coordinates({0: arr})# doctest: +SKIP
        """
        render_params = render_params or {}
        self._coordinates_dict = arr_dict

        buffers = []
        coordinates_meta = dict()
        for index, arr in self._coordinates_dict.items():
            buffers.append(arr.astype('f4').tobytes())
            coordinates_meta[index] = index
        msg = {
            'type': 'binary_single',
            'data': coordinates_meta,
        }
        if movie_making:
            msg['movie_making'] = movie_making
            msg['render_params'] = render_params

        self.send(msg, buffers=buffers)

    @observe('frame')
    def _on_frame_changed(self, change):
        self._set_coordinates(change['new'])

    def _ipython_display_(self, **kwargs):
        try:
            # ipywidgets < 8
            super()._ipython_display_(**kwargs)
        except AttributeError:
            display(super()._repr_mimebundle_(), raw=True)