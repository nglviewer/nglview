from __future__ import print_function, absolute_import
import threading
import time
import base64
import uuid
import json
import numpy as np
from IPython.display import display
from ipywidgets import Box, DOMWidget, widget_image
from traitlets import (Unicode, Bool, Dict, List, Int, Integer, observe,
                       CaselessStrEnum)

from .utils import py_utils, js_utils, widget_utils
from .utils.py_utils import (seq_to_string, _camelize_dict, FileManager,
                             get_repr_names_from_dict, encode_base64,
                             _update_url)
from .player import TrajectoryPlayer
from . import interpolate
from .stage import Stage
from .shape import Shape
from .viewer_control import ViewerControl
from .representation import RepresentationControl

from .adaptor import Structure, Trajectory
from .config import BACKENDS
from .parameters import REPRESENTATION_NAME_PAIRS
from .remote_thread import RemoteCallThread

__all__ = ['NGLWidget', 'ComponentViewer']
__frontend_version__ = '0.5.4-dev.8' # must match to js/package.json and js/src/widget_ngl.js


def _add_repr_method_shortcut(self, other):
    from types import MethodType

    def make_func_add(rep):
        """return a new function object
        """

        def func(this, selection='all', **kwargs):
            """
            """
            self.add_representation(
                repr_type=rep[1], selection=selection, **kwargs)

        func.__doc__ = """Shortcut for `add_representation` method

        Examples
        --------
        >>> view.add_{name}()
        >>> # is equal to
        >>> view.add_representation('{name}')
        """.format(name=rep[0])
        return func

    def make_func_remove(rep):
        """return a new function object
        """

        def func(this, **kwargs):
            """
            """
            self._remove_representations_by_name(repr_name=rep[1], **kwargs)

        return func

    def make_func_update(rep):
        """return a new function object
        """

        def func(this, **kwargs):
            """
            """
            self._update_representations_by_name(repr_name=rep[1], **kwargs)

        return func

    for rep in REPRESENTATION_NAME_PAIRS:
        for make_func, root_fn in [(make_func_add, 'add'), (make_func_update,
                                                            'update'),
                                   (make_func_remove, 'remove')]:
            func = make_func(rep)
            fn = '_'.join((root_fn, rep[0]))
            setattr(self, fn, MethodType(func, other))


class NGLWidget(DOMWidget):
    _view_name = Unicode("NGLView").tag(sync=True)
    _view_module = Unicode("nglview-js-widgets").tag(sync=True)
    _view_module_version = Unicode(__frontend_version__).tag(sync=True)
    _ngl_version = Unicode().tag(sync=True)
    # _model_name = Unicode("NGLView").tag(sync=True)
    # _model_module = Unicode("nglview-js-widgets").tag(sync=True)
    _image_data = Unicode().tag(sync=True)
    # use Integer here, because mdtraj uses a long datatype here on Python-2.7
    frame = Integer().tag(sync=True)
    count = Integer(1).tag(sync=True)
    background = Unicode('white').tag(sync=True)
    loaded = Bool(False).tag(sync=False)
    picked = Dict().tag(sync=True)
    n_components = Int(0).tag(sync=True)
    orientation = List().tag(sync=True)
    _scene_position = Dict().tag(sync=True)
    _scene_rotation = Dict().tag(sync=True)
    _first_time_loaded = Bool(True).tag(sync=False)
    # hack to always display movie
    _n_dragged_files = Int().tag(sync=True)
    # TODO: remove _parameters?
    _parameters = Dict().tag(sync=False)
    _full_stage_parameters = Dict().tag(sync=True)
    _original_stage_parameters = Dict().tag(sync=True)
    _coordinates_dict = Dict().tag(sync=False)
    _camera_str = CaselessStrEnum(
        ['perspective', 'orthographic'], default_value='orthographic').tag(
            sync=True)
    _repr_dict = Dict().tag(sync=False)
    _ngl_component_ids = List().tag(sync=False)
    _ngl_component_names = List().tag(sync=False)
    _already_constructed = Bool(False).tag(sync=False)
    _ngl_msg = None
    _send_binary = Bool(True).tag(sync=False)
    _init_gui = Bool(False).tag(sync=False)
    _hold_image = Bool(False).tag(sync=False)

    def __init__(self,
                 structure=None,
                 representations=None,
                 parameters=None,
                 **kwargs):
        super(NGLWidget, self).__init__(**kwargs)

        self._gui = None
        self._init_gui = kwargs.pop('gui', False)
        self._theme = kwargs.pop('theme', 'default')
        self._widget_image = widget_image.Image()
        self._widget_image.width = 900.
        self._image_array = []
        # do not use _displayed_callbacks since there is another Widget._display_callbacks
        self._event = threading.Event()
        self._ngl_displayed_callbacks_before_loaded = []
        self._ngl_displayed_callbacks_after_loaded = []
        _add_repr_method_shortcut(self, self)
        self.shape = Shape(view=self)
        self.stage = Stage(view=self)
        self.control = ViewerControl(view=self)
        self._handle_msg_thread = threading.Thread(
            target=self.on_msg, args=(self._ngl_handle_msg, ))
        # # register to get data from JS side
        self._handle_msg_thread.daemon = True
        self._handle_msg_thread.start()
        self._remote_call_thread = RemoteCallThread(
            self, registered_funcs=['loadFile', 'replaceStructure'])
        self._remote_call_thread.start()
        self._trajlist = []
        self._ngl_component_ids = []
        if parameters:
            self.parameters = parameters
        if isinstance(structure, Trajectory):
            name = py_utils.get_name(structure, kwargs)
            self.add_trajectory(structure, name=name)
        elif isinstance(structure, (list, tuple)):
            trajectories = structure
            for trajectory in trajectories:
                name = py_utils.get_name(trajectory, kwargs)
                self.add_trajectory(trajectory, name=name)
        else:
            if structure is not None:
                self.add_structure(structure, **kwargs)

        if representations:
            self._init_representations = representations
        else:
            self._init_representations = [{
                "type": "cartoon",
                "params": {
                    "sele": "polymer"
                }
            }, {
                "type": "ball+stick",
                "params": {
                    "sele": "hetero OR mol"
                }
            }, {
                "type": "ball+stick",
                "params": {
                    "sele": "not protein and not nucleic"
                }
            }]

        # keep track but making copy
        if structure is not None:
            self._representations = self._init_representations[:]

        self._set_unsync_camera()
        selector = 'nglviewHolder' + str(id(self))
        self._remote_call(
            'setSelector', target='Widget', args=[
                selector,
            ])
        self.player = TrajectoryPlayer(self)
        self._already_constructed = True

    @property
    def parameters(self):
        return self._parameters

    @parameters.setter
    def parameters(self, params):
        params = _camelize_dict(params)
        self._parameters = params
        self._remote_call(
            'setParameters', target='Widget', args=[
                params,
            ])

    @property
    def camera(self):
        return self._camera_str

    @camera.setter
    def camera(self, value):
        """
        
        Parameters
        ----------
        value : str, {'perspective', 'orthographic'}
        """
        self._camera_str = value
        # use _remote_call so this function can be called right after
        # self is displayed
        self._remote_call(
            "setParameters",
            target='Stage',
            kwargs=dict(cameraType=self._camera_str))

    def _request_stage_parameters(self):
        self._remote_call('requestUpdateStageParameters', target='Widget')

    @observe('picked')
    def _on_picked(self, change):
        picked = change['new']
        if self.player.widget_picked is not None:
            self.player.widget_picked.value = json.dumps(picked)

    @observe('background')
    def _update_background_color(self, change):
        color = change['new']
        self.parameters = dict(background_color=color)

    @observe('_n_dragged_files')
    def on_update_dragged_file(self, change):
        if change['new'] - change['old'] == 1:
            self._ngl_component_ids.append(uuid.uuid4())

    @observe('n_components')
    def _handle_n_components_changed(self, change):
        if self.player.widget_repr is not None:
            component_slider = widget_utils.get_widget_by_name(
                self.player.widget_repr, 'component_slider')

            if change['new'] - 1 >= component_slider.min:
                component_slider.max = change['new'] - 1

            component_dropdown = widget_utils.get_widget_by_name(
                self.player.widget_repr, 'component_dropdown')
            component_dropdown.options = tuple(self._ngl_component_names)

            if change['new'] == 0:
                component_dropdown.options = tuple([' '])
                component_dropdown.value = ' '

                component_slider.max = 0

                reprlist_choices = widget_utils.get_widget_by_name(
                    self.player.widget_repr, 'reprlist_choices')
                reprlist_choices.options = tuple([' '])

                repr_slider = widget_utils.get_widget_by_name(
                    self.player.widget_repr, 'repr_slider')
                repr_slider.max = 0

                repr_name_text = widget_utils.get_widget_by_name(
                    self.player.widget_repr, 'repr_name_text')
                repr_selection = widget_utils.get_widget_by_name(
                    self.player.widget_repr, 'repr_selection')
                repr_name_text.value = ' '
                repr_selection.value = ' '

    @observe('_repr_dict')
    def _handle_repr_dict_changed(self, change):
        if self.player.widget_repr is not None:
            repr_slider = widget_utils.get_widget_by_name(
                self.player.widget_repr, 'repr_slider')
            component_slider = widget_utils.get_widget_by_name(
                self.player.widget_repr, 'component_slider')
            repr_name_text = widget_utils.get_widget_by_name(
                self.player.widget_repr, 'repr_name_text')
            repr_selection = widget_utils.get_widget_by_name(
                self.player.widget_repr, 'repr_selection')
            reprlist_choices = widget_utils.get_widget_by_name(
                self.player.widget_repr, 'reprlist_choices')
            repr_names = get_repr_names_from_dict(self._repr_dict,
                                                  component_slider.value)

            if change['new'] == {'c0': {}}:
                repr_selection.value = ''
            else:
                options = tuple(
                    str(i) + '-' + name for (i, name) in enumerate(repr_names))
                reprlist_choices.options = options

                try:
                    value = reprlist_choices.options[repr_slider.value]
                    if isinstance(value, tuple):
                        # https://github.com/jupyter-widgets/ipywidgets/issues/1512
                        value = value[0]
                    reprlist_choices.value = value
                except IndexError:
                    if repr_slider.value == 0:
                        # works fine with ipywidgets 5.2.2
                        reprlist_choices.options = tuple([
                            ' ',
                        ])
                        reprlist_choices.value = ' '
                    else:
                        reprlist_choices.value = reprlist_choices.options[
                            repr_slider.value - 1]

                # e.g: 0-cartoon
                repr_name_text.value = reprlist_choices.value.split('-')[-1]

                repr_slider.max = len(repr_names) - 1 if len(
                    repr_names) >= 1 else len(repr_names)

    def _update_count(self):
        self.count = max(
            traj.n_frames for traj in self._trajlist
            if hasattr(traj, 'n_frames'))

    def _wait_until_finished(self, timeout=0.0001):
        # NGL need to send 'finished' signal to
        # backend
        self._event.clear()
        while True:
            # idle to make room for waiting for
            # "finished" event sent from JS
            time.sleep(timeout)
            if self._event.is_set():
                # if event is set from another thread
                # break while True
                break

    def _run_on_another_thread(self, func, *args):
        # use `event` to singal
        # func(*args)
        thread = threading.Thread(
            target=func,
            args=args, )
        thread.daemon = True
        thread.start()
        return thread

    @observe('loaded')
    def on_loaded(self, change):
        # trick for firefox on Linux
        time.sleep(0.1)

        if change['new']:
            self._fire_callbacks(self._ngl_displayed_callbacks_before_loaded)

    def _fire_callbacks(self, callbacks):
        def _call(event):
            for callback in callbacks:
                callback(self)
                if callback._method_name == 'loadFile':
                    self._wait_until_finished()
        self._run_on_another_thread(_call, self._event)

    def _refresh_render(self):
        """useful when you update coordinates for a single structure.

        Notes
        -----
        If you are visualizing a trajectory with more than 1 frame, you can use the
        player slider to trigger the refreshing.
        """
        current_frame = self.frame
        self.frame = int(1E6)
        self.frame = current_frame

    def sync_view(self):
        """call this if you want to sync multiple views of a single viewer

        Note: unstable feature
        """
        self._fire_callbacks(self._ngl_displayed_callbacks_after_loaded)

    def _ipython_display_(self, **kwargs):
        super(NGLWidget, self)._ipython_display_(**kwargs)
        if self._first_time_loaded:
            self._first_time_loaded = False
        else:
            time.sleep(0.1)
            self.sync_view()
        if self._init_gui:
            if self._gui is None:
                self._gui = self.player._display()
            display(self._gui)

        if self._theme in ['dark', 'oceans16']:
            from nglview import theme
            display(theme.oceans16())
            self._remote_call('cleanOutput', target='Widget')

    def display(self, gui=False, use_box=False):
        if gui:
            if use_box:
                box = Box([self, self.player._display()])
                box._gui_style = 'row'
                return box
            else:
                display(self)
                display(self.player._display())
                return None
        else:
            return self

    def _set_size(self, w, h):
        '''

        Parameters
        ----------
        w, h : float or str

        Examples
        --------
        >>> import nglview; view = nglview.demo()
        >>> view._set_size(100, 100)
        >>> view._set_size('100px', '100px')
        >>> view._set_size('50%', '50%')
        '''
        self._remote_call('setSize', target='Widget', args=[w, h])

    def _set_draggable(self, yes=True):
        if yes:
            self._remote_call(
                'setDraggable', target='Widget', args=[
                    '',
                ])

        else:
            self._remote_call(
                'setDraggable', target='Widget', args=[
                    'destroy',
                ])

    def _set_sync_frame(self):
        self._remote_call("setSyncFrame", target="Widget")

    def _set_unsync_frame(self):
        self._remote_call("setUnSyncFrame", target="Widget")

    def _set_sync_camera(self):
        self._remote_call("setSyncCamera", target="Widget")

    def _set_unsync_camera(self):
        self._remote_call("setUnSyncCamera", target="Widget")

    def _set_delay(self, delay):
        """unit of millisecond
        """
        self._remote_call(
            "setDelay", target="Widget", args=[
                delay,
            ])

    def _set_spin(self, axis, angle):
        self._remote_call('setSpin', target='Stage', args=[axis, angle])

    def _set_selection(self, selection, component=0, repr_index=0):
        self._remote_call(
            "setSelection",
            target='Representation',
            args=[selection],
            kwargs=dict(component_index=component, repr_index=repr_index))

    def _set_color_by_residue(self, colors, component_index=0, repr_index=0):
        self._remote_call(
            'setColorByResidue',
            target='Widget',
            args=[colors, component_index, repr_index])

    def color_by(self, color_scheme, component=0):
        '''update color for all representations of given component

        Notes
        -----
        Unstable feature

        Parameters
        ----------
        color_scheme : str
        component : int, default 0
            component index

        Examples
        --------
        >>> import nglview
        >>> view = nglview.demo()
        >>> # component 0
        >>> view.color_by('atomindex')

        >>> # component 1
        >>> view.color_by('atomindex', component=1)
        '''
        repr_names = get_repr_names_from_dict(self._repr_dict, component)

        for index, _ in enumerate(repr_names):
            self.update_representation(
                component=component,
                repr_index=index,
                color_scheme=color_scheme)

    @property
    def representations(self):
        return self._representations

    @representations.setter
    def representations(self, reps):
        self._representations = reps[:]
        for index in range(len(self._ngl_component_ids)):
            self.set_representations(reps)

    def update_representation(self, component=0, repr_index=0, **parameters):
        """

        Parameters
        ----------
        component : int, default 0
            component index
        repr_index : int, default 0
            representation index for given component
        parameters : dict
        """
        parameters = _camelize_dict(parameters)
        kwargs = dict(component_index=component, repr_index=repr_index)
        kwargs.update(parameters)

        self._remote_call(
            'setParameters', target='Representation', kwargs=kwargs)
        self._update_repr_dict()

    def _update_repr_dict(self):
        """ Send a request to fronend to send representation parameters
        back.

        # TODO: sync or async
        """
        self._remote_call('request_repr_dict', target='Widget')

    def set_representations(self, representations, component=0):
        """
        
        Parameters
        ----------
        representations : list of dict
        """
        self.clear_representations(component=component)

        for params in representations:
            assert isinstance(params, dict), 'params must be a dict'
            kwargs = params['params']
            kwargs.update({'component_index': component})
            self._remote_call(
                'addRepresentation',
                target='compList',
                args=[
                    params['type'],
                ],
                kwargs=kwargs)

    def _remove_representation(self, component=0, repr_index=0):
        self._remote_call(
            'removeRepresentation',
            target='Widget',
            args=[component, repr_index])

    def _remove_representations_by_name(self, repr_name, component=0):
        self._remote_call(
            'removeRepresentationsByName',
            target='Widget',
            args=[repr_name, component])

    def _update_representations_by_name(self, repr_name, component=0,
                                        **kwargs):
        kwargs = _camelize_dict(kwargs)
        self._remote_call(
            'updateRepresentationsByName',
            target='Widget',
            args=[repr_name, component],
            kwargs=kwargs)

    def _display_repr(self, component=0, repr_index=0, name=None):
        c = 'c' + str(component)
        r = str(repr_index)

        try:
            name = self._repr_dict[c][r]['type']
        except KeyError:
            name = ''

        return RepresentationControl(self, component, repr_index, name=name)

    def _set_coordinates(self, index):
        '''update coordinates for all trajectories at index-th frame
        '''
        if self._trajlist:
            coordinates_dict = {}
            for trajectory in self._trajlist:
                traj_index = self._ngl_component_ids.index(trajectory.id)

                try:
                    if trajectory.shown:
                        if self.player.interpolate:
                            t = self.player.iparams.get('t', 0.5)
                            step = self.player.iparams.get('step', 1)
                            coordinates_dict[traj_index] = interpolate.linear(
                                index, t=t, traj=trajectory, step=step)
                        else:
                            coordinates_dict[
                                traj_index] = trajectory.get_coordinates(index)
                    else:
                        coordinates_dict[traj_index] = np.empty(
                            (0), dtype='f4')
                except (IndexError, ValueError):
                    coordinates_dict[traj_index] = np.empty((0), dtype='f4')

            self.coordinates_dict = coordinates_dict
        else:
            print("no trajectory available")

    @property
    def coordinates_dict(self):
        """

        Returns
        -------
        out : dict of numpy 3D-array, dtype='f4'
            coordinates of trajectories at current frame
        """
        return self._coordinates_dict

    @coordinates_dict.setter
    def coordinates_dict(self, arr_dict):
        self._coordinates_dict = arr_dict

        if not self._send_binary:
            # send base64
            encoded_coordinates_dict = dict(
                (k, encode_base64(v))
                for (k, v) in self._coordinates_dict.items())
            mytime = time.time() * 1000
            self.send({
                'type': 'base64_single',
                'data': encoded_coordinates_dict,
                'mytime': mytime
            })
        else:
            # send binary
            buffers = []
            coordinates_meta = dict()
            for index, arr in self._coordinates_dict.items():
                buffers.append(arr.astype('f4').tobytes())
                coordinates_meta[index] = index
            mytime = time.time() * 1000
            self.send(
                {
                    'type': 'binary_single',
                    'data': coordinates_meta,
                    'mytime': mytime
                },
                buffers=buffers)

    @observe('frame')
    def on_frame_changed(self, change):
        """set and send coordinates at current frame
        """
        self._set_coordinates(self.frame)

    def clear(self, *args, **kwargs):
        '''shortcut of `clear_representations`
        '''

        self.clear_representations(*args, **kwargs)

    def clear_representations(self, component=0):
        '''clear all representations for given component

        Parameters
        ----------
        component : int, default 0 (first model)
            You need to keep track how many components you added.
        '''
        self._remote_call(
            "removeAllRepresentations",
            target='compList',
            kwargs={'component_index': component})

    @_update_url
    def _add_shape(self, shapes, name='shape'):
        """add shape objects

        TODO: update doc, caseless shape keyword

        Parameters
        ----------
        shapes : list of tuple
        name : str, default 'shape'
            name of given shape

        Notes
        -----
        Supported shape: 'mesh', 'sphere', 'ellipsoid', 'cylinder', 'cone', 'arrow'.
        
        See also
        --------
        {ngl_url}

        Examples
        --------
        >>> import nglview
        >>> view = nglview.demo()
        >>> sphere = ('sphere', [0, 0, 9], [1, 0, 0], 1.5)
        >>> arrow = ('arrow', [1, 2, 7 ], [30, 3, 3], [1, 0, 1], 1.0)
        >>> view._add_shape([sphere, arrow], name='my_shape')
        """

        self._remote_call('addShape', target='Widget', args=[name, shapes])

    @_update_url
    def add_representation(self, repr_type, selection='all', **kwargs):
        '''Add structure representation (cartoon, licorice, ...) for given atom selection.

        Parameters
        ----------
        repr_type : str
            type of representation. Please see {ngl_url} for further info.
        selection : str or 1D array (atom indices) or any iterator that returns integer, default 'all'
            atom selection
        **kwargs: additional arguments for representation

        Example
        -------
        >>> import nglview as nv
        >>> import pytraj 
        >>> t = pytraj.datafiles.load_tz2()
        >>> w = nv.show_pytraj(t)
        >>> w.add_representation('cartoon', selection='protein', color='blue')
        >>> w.add_representation('licorice', selection=[3, 8, 9, 11], color='red')
        >>> w # doctest: +SKIP

        Notes
        -----
        User can also use shortcut

        >>> selection = 'protein'
        >>> w.add_cartoon(selection) # w.add_representation('cartoon', selection)
        '''
        if repr_type == 'surface':
            if 'useWorker' not in kwargs:
                kwargs['useWorker'] = False

        # avoid space sensitivity
        repr_type = repr_type.strip()
        # overwrite selection
        selection = seq_to_string(selection).strip()

        # make copy
        kwargs2 = _camelize_dict(kwargs)

        if 'component' in kwargs2:
            component = kwargs2.pop('component')
        else:
            component = 0

        for k, v in kwargs2.items():
            try:
                kwargs2[k] = v.strip()
            except AttributeError:
                # e.g.: opacity=0.4
                kwargs2[k] = v

        d = {'params': {'sele': selection}}
        d['type'] = repr_type
        d['params'].update(kwargs2)

        params = d['params']
        params.update({'component_index': component})
        self._remote_call(
            'addRepresentation',
            target='compList',
            args=[
                d['type'],
            ],
            kwargs=params)

    def center(self, *args, **kwargs):
        """alias of `center_view`
        """
        self.center_view(*args, **kwargs)

    def center_view(self, selection='*', duration=0, component=0):
        """center view for given atom selection

        Examples
        --------
        view.center_view(selection='1-4')
        """
        self._remote_call(
            'autoView',
            target='compList',
            args=[selection, duration],
            kwargs={'component_index': component})

    @observe('_image_data')
    def _on_render_image(self, change):
        '''update image data to widget_image

        Notes
        -----
        method name might be changed
        '''
        self._widget_image._b64value = change['new']
        if self._hold_image:
            self._image_array.append(change['new'])

    def render_image(self,
                     frame=None,
                     factor=4,
                     antialias=True,
                     trim=False,
                     transparent=False):
        """render and get image as ipywidgets.widget_image.Image

        Parameters
        ----------
        frame : int or None, default None
            if None, use current frame
            if specified, use this number.
        factor : int, default 4
            quality of the image, higher is better
        antialias : bool, default True
        trim : bool, default False
        transparent : bool, default False

        Examples
        --------
            # tell NGL to render send image data to notebook.
            view.render_image()
            
            # make sure to call `get_image` method
            view.get_image()

        Notes
        -----
        You need to call `render_image` and `get_image` in different notebook's Cells
        """
        if frame is not None:
            self.frame = frame
        params = dict(
            factor=factor,
            antialias=antialias,
            trim=trim,
            transparent=transparent)
        self._remote_call('_exportImage', target='Widget', kwargs=params)

    def download_image(self,
                       filename='screenshot.png',
                       factor=4,
                       antialias=True,
                       trim=False,
                       transparent=False):
        """render and download scence at current frame

        Parameters
        ----------
        filename : str, default 'screenshot.png'
        factor : int, default 4
            quality of the image, higher is better
        antialias : bool, default True
        trim : bool, default False
        transparent : bool, default False
        """
        params = dict(
            factor=factor,
            antialias=antialias,
            trim=trim,
            transparent=transparent)
        self._remote_call(
            '_downloadImage',
            target='Widget',
            args=[
                filename,
            ],
            kwargs=params)

    def _get_movie_maker(self, in_memory=True, **kwargs):
        ''' create MovieMaker object

        Examples
        --------
        >>> movie = view._get_movie_maker(output='my.gif') # doctest: +SKIP
        ... movie.make()


        Notes
        -----
        We only test with imageio 1.6 and moviepy 0.2.2.11
        Good luck.
        '''
        from nglview.contrib.movie import MovieMaker

        if 'in_memory' not in kwargs:
            kwargs['in_memory'] = in_memory

        movie_maker = MovieMaker(self, **kwargs)
        return movie_maker

    def _ngl_handle_msg(self, widget, msg, buffers):
        """store message sent from Javascript.

        How? use view.on_msg(get_msg)
        """
        self._ngl_msg = msg

        msg_type = self._ngl_msg.get('type')
        if msg_type == 'request_frame':
            self.frame += self.player.step
            if self.frame >= self.count:
                self.frame = 0
            elif self.frame < 0:
                self.frame = self.count - 1
        elif msg_type == 'repr_parameters':
            data_dict = self._ngl_msg.get('data')
            name = data_dict.pop('name') + '\n'
            selection = data_dict.get('sele', '') + '\n'
            # json change True to true
            data_dict_json = json.dumps(data_dict).replace(
                'true', 'True').replace('false', 'False')
            data_dict_json = data_dict_json.replace('null', '"null"')

            if self.player.widget_repr is not None:
                # TODO: refactor
                repr_name_text = widget_utils.get_widget_by_name(
                    self.player.widget_repr, 'repr_name_text')
                repr_selection = widget_utils.get_widget_by_name(
                    self.player.widget_repr, 'repr_selection')
                repr_name_text.value = name
                repr_selection.value = selection

        elif msg_type == 'request_loaded':
            if not self.loaded:
                # trick to trigger observe loaded
                # so two viewers can have the same representations
                self.loaded = False
            self.loaded = msg.get('data')
        elif msg_type == 'request_repr_dict':
            self._repr_dict = self._ngl_msg.get('data')
        elif msg_type == 'stage_parameters':
            self._full_stage_parameters = msg.get('data')
        elif msg_type == 'async_message':
            if msg.get('data') == 'ok':
                self._event.set()

    def _request_repr_parameters(self, component=0, repr_index=0):
        self._remote_call(
            'requestReprParameters',
            target='Widget',
            args=[component, repr_index])

    def add_structure(self, structure, **kwargs):
        '''add structure to view

        Parameters
        ----------
        structure : nglview.Structure object

        Examples
        --------
        >>> view.add_trajectory(traj0) # doctest: +SKIP
        ... view.add_trajectory(traj1)
        ... # then add Structure
        ... view.add_structure(s)

        See Also
        --------
        nglview.NGLWidget.add_component
        '''
        if not isinstance(structure, Structure):
            raise ValueError(
                '{} is not an instance of Structure'.format(structure))
        self._load_data(structure, **kwargs)
        self._ngl_component_ids.append(structure.id)
        if self.n_components > 1:
            self.center_view(component=len(self._ngl_component_ids) - 1)
        self._update_component_auto_completion()

    def add_trajectory(self, trajectory, **kwargs):
        '''add new trajectory to `view`

        Parameters
        ----------
        trajectory: nglview.Trajectory or its derived class or 
            pytraj.Trajectory-like, mdtraj.Trajectory or MDAnalysis objects

        See Also
        --------
        nglview.NGLWidget.add_component

        Examples
        --------
        >>> import nglview as nv, pytraj as pt
        >>> traj = pt.load(nv.datafiles.TRR, nv.datafiles.PDB)
        >>> view = nv.show_pytraj(traj)
        >>> # show view first
        >>> view # doctest: +SKIP
        >>> # add new Trajectory
        >>> traj2 = pt.datafiles.load_tz2()
        >>> view.add_trajectory(traj2)
        '''
        backends = BACKENDS

        package_name = trajectory.__module__.split('.')[0]

        if package_name in backends:
            trajectory = backends[package_name](trajectory)
        else:
            trajectory = trajectory

        self._load_data(trajectory, **kwargs)
        setattr(trajectory, 'shown', True)
        self._trajlist.append(trajectory)
        self._update_count()
        self._ngl_component_ids.append(trajectory.id)
        self._update_component_auto_completion()

    def add_pdbid(self, pdbid):
        '''add new Structure view by fetching pdb id from rcsb

        Examples
        --------
        >>> import nglview
        >>> view = nglview.NGLWidget()
        >>> view.add_pdbid('1tsu')
        >>> # which is equal to 
        >>> # view.add_component('rcsb://1tsu.pdb')
        '''
        self.add_component('rcsb://{}.pdb'.format(pdbid))

    def add_component(self, filename, **kwargs):
        '''add component from file/trajectory/struture

        Parameters
        ----------
        filename : str or Trajectory or Structure or their derived class or url
        **kwargs : additional arguments, optional

        Examples
        --------
        >>> import nglview
        >>> view = nglview.NGLWidget()
        >>> view # doctest: +SKIP
        ... filename = 'somefile.ccp4'
        ... view.add_component(filename)

        Notes
        -----
        If you want to load binary file such as density data, mmtf format, it is
        faster to load file from current or subfolder.
        '''
        self._load_data(filename, **kwargs)
        # assign an ID
        self._ngl_component_ids.append(str(uuid.uuid4()))
        self._update_component_auto_completion()

    def _load_data(self, obj, **kwargs):
        '''

        Parameters
        ----------
        obj : nglview.Structure or any object having 'get_structure_string' method or
              string buffer (open(fn).read())
        '''
        kwargs2 = _camelize_dict(kwargs)

        try:
            is_url = FileManager(obj).is_url
        except NameError:
            is_url = False

        if 'defaultRepresentation' not in kwargs2:
            kwargs2['defaultRepresentation'] = True

        if not is_url:
            if hasattr(obj, 'get_structure_string'):
                blob = obj.get_structure_string()
                kwargs2['ext'] = obj.ext
                passing_buffer = True
                binary = False
            else:
                fh = FileManager(
                    obj,
                    ext=kwargs.get('ext'),
                    compressed=kwargs.get('compressed'))
                # assume passing string
                blob = fh.read()
                passing_buffer = not fh.use_filename

                if fh.ext is None and passing_buffer:
                    raise ValueError('must provide extension')

                kwargs2['ext'] = fh.ext
                binary = fh.is_binary
                use_filename = fh.use_filename

            if binary and not use_filename:
                # send base64
                blob = base64.b64encode(blob).decode('utf8')
            blob_type = 'blob' if passing_buffer else 'path'
            args = [{'type': blob_type, 'data': blob, 'binary': binary}]
        else:
            # is_url
            blob_type = 'url'
            url = obj
            args = [{'type': blob_type, 'data': url, 'binary': False}]

        name = py_utils.get_name(obj, kwargs2)
        self._ngl_component_names.append(name)
        self._remote_call(
            "loadFile", target='Stage', args=args, kwargs=kwargs2)

    def remove_component(self, component_id):
        """remove component by its uuid

        Examples
        --------
        >>> view.add_trajectory(traj0) # doctest: +SKIP
        ... view.add_trajectory(traj1)
        ... view.add_struture(structure)
        ... # remove last component
        ... view.remove_component(view._ngl_component_ids[-1])
        """
        self._clear_component_auto_completion()
        if self._trajlist:
            for traj in self._trajlist:
                if traj.id == component_id:
                    self._trajlist.remove(traj)
        component_index = self._ngl_component_ids.index(component_id)
        self._ngl_component_ids.remove(component_id)
        self._ngl_component_names.pop(component_index)

        self._remote_call(
            'removeComponent', target='Stage', args=[
                component_index,
            ])

        self._update_component_auto_completion()

    def _remote_call(self,
                     method_name,
                     target='Widget',
                     args=None,
                     kwargs=None):
        """call NGL's methods from Python.
        
        Parameters
        ----------
        method_name : str
        target : str, {'Stage', 'Viewer', 'compList', 'StructureComponent'}
        args : list
        kwargs : dict
            if target is 'compList', "component_index" could be passed
            to specify which component will call the method.

        Examples
        --------
        view._remote_call('loadFile', args=['1L2Y.pdb'],
                          target='Stage', kwargs={'defaultRepresentation': True})

        # perform autoView for 1st component
        # JS code
        # component = Stage.compList[1];
        # component.autoView('*', 200)

        # python
        view._remote_call('autoView',
                          target='component',
                          args=['*', 200],
                          kwargs={'component_index': 1})
        """
        args = [] if args is None else args
        kwargs = {} if kwargs is None else kwargs

        msg = {}

        if 'component_index' in kwargs:
            msg['component_index'] = kwargs.pop('component_index')
        if 'repr_index' in kwargs:
            msg['repr_index'] = kwargs.pop('repr_index')

        msg['target'] = target
        msg['type'] = 'call_method'
        msg['methodName'] = method_name
        msg['args'] = args
        msg['kwargs'] = kwargs

        def callback(widget, msg=msg):
            widget.send(msg)

        callback._method_name = method_name

        if self.loaded:
            self._remote_call_thread.q.append(callback)
        else:
            # send later
            # all callbacks will be called right after widget is loaded
            self._ngl_displayed_callbacks_before_loaded.append(callback)

        self._ngl_displayed_callbacks_after_loaded.append(callback)

    def _get_traj_by_id(self, itsid):
        """return nglview.Trajectory or its derived class object
        """
        for traj in self._trajlist:
            if traj.id == itsid:
                return traj
        return None

    def hide(self, indices):
        """set invisibility for given component/struture/trajectory (by their indices)
        """
        traj_ids = set(traj.id for traj in self._trajlist)

        for index in indices:
            comp_id = self._ngl_component_ids[index]
            if comp_id in traj_ids:
                traj = self._get_traj_by_id(comp_id)
                traj.shown = False
            self._remote_call(
                "setVisibility",
                target='compList',
                args=[
                    False,
                ],
                kwargs={'component_index': index})

    def show(self, **kwargs):
        """shortcut of `show_only`
        """
        self.show_only(**kwargs)

    def show_only(self, indices='all'):
        """set visibility for given components (by their indices)

        Parameters
        ----------
        indices : {'all', array-like}, component index, default 'all'
        """
        traj_ids = set(traj.id for traj in self._trajlist)

        if indices == 'all':
            indices_ = set(range(self.n_components))
        else:
            indices_ = set(indices)

        for index, comp_id in enumerate(self._ngl_component_ids):
            if comp_id in traj_ids:
                traj = self._get_traj_by_id(comp_id)
            else:
                traj = None
            if index in indices_:
                args = [
                    True,
                ]
                if traj is not None:
                    traj.shown = True
            else:
                args = [
                    False,
                ]
                if traj is not None:
                    traj.shown = False

            self._remote_call(
                "setVisibility",
                target='compList',
                args=args,
                kwargs={'component_index': index})

    def _js_console(self):
        self.send(dict(type='get', data='any'))

    def _get_full_params(self):
        self.send(dict(type='get', data='parameters'))

    def _display_image(self):
        '''for testing
        '''
        from IPython import display
        return display.Image(self._image_data)

    def _clear_component_auto_completion(self):
        for index, _ in enumerate(self._ngl_component_ids):
            name = 'component_' + str(index)
            delattr(self, name)

    def _update_component_auto_completion(self):
        trajids = [traj.id for traj in self._trajlist]

        for index, cid in enumerate(self._ngl_component_ids):
            comp = ComponentViewer(self, index)
            name = 'component_' + str(index)
            setattr(self, name, comp)

            if cid in trajids:
                traj_name = 'trajectory_' + str(trajids.index(cid))
                setattr(self, traj_name, comp)

    def __getitem__(self, index):
        """return ComponentViewer
        """
        postive_index = py_utils.get_positive_index(
            index, len(self._ngl_component_ids))
        return ComponentViewer(self, postive_index)

    def __iter__(self):
        """return ComponentViewer
        """
        for i, _ in enumerate(self._ngl_component_ids):
            yield self[i]

    def detach(self, split=False):
        """detach player from its original container.

        Parameters
        ----------
        split : bool, default False
            if True, resize notebook then move it to the right of its container
        """
        if not self.loaded:
            raise RuntimeError("must display view first")

        # resize notebook first
        # width of the dialog will be calculated based on notebook container offset
        if split:
            # rename
            js_utils._move_notebook_to_the_right()
        self._remote_call('setDialog', target='Widget')


class ComponentViewer(object):
    """Convenient attribute for NGLWidget. See example below.

    Examples
    --------
    >>> view = nv.NGLWidget() # doctest: +SKIP
    ... view.add_trajectory(traj) # traj is a component 0
    ... view.add_component(filename) # component 1
    ... view.component_0.clear_representations()
    ... view.component_0.add_cartoon()
    ... view.component_1.add_licorice()
    ... view.remove_component(view.comp1.id)
    """

    def __init__(self, view, index):
        self._view = view
        self._index = index
        _add_repr_method_shortcut(self, self._view)
        self._borrow_attribute(self._view, [
            'clear_representations', '_remove_representations_by_name',
            '_update_representations_by_name', 'center_view', 'center',
            'clear', 'set_representations'
        ], ['get_structure_string', 'get_coordinates', 'n_frames'])

    @property
    def id(self):
        return self._view._ngl_component_ids[self._index]

    def hide(self):
        """set invisibility for given components (by their indices)
        """
        self._view._remote_call(
            "setVisibility",
            target='compList',
            args=[
                False,
            ],
            kwargs={'component_index': self._index})
        traj = self._view._get_traj_by_id(self.id)
        if traj is not None:
            traj.shown = False

    def show(self):
        """set invisibility for given components (by their indices)
        """
        self._view._remote_call(
            "setVisibility",
            target='compList',
            args=[
                True,
            ],
            kwargs={'component_index': self._index})

        traj = self._view._get_traj_by_id(self.id)
        if traj is not None:
            traj.shown = True

    def add_representation(self, repr_type, selection='all', **kwargs):
        kwargs['component'] = self._index
        self._view.add_representation(
            repr_type=repr_type, selection=selection, **kwargs)

    def _borrow_attribute(self, view, attributes, trajectory_atts=None):
        from functools import partial
        from types import MethodType

        traj = view._get_traj_by_id(self.id)

        for attname in attributes:
            view_att = getattr(view, attname)
            setattr(self, '_' + attname, MethodType(view_att, view))
            self_att = partial(getattr(view, attname), component=self._index)
            setattr(self, attname, self_att)

        if traj is not None and trajectory_atts is not None:
            for attname in trajectory_atts:
                traj_att = getattr(traj, attname)
                setattr(self, attname, traj_att)
