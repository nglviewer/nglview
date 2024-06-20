from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional
try:
    import moviepy.editor as mpy
except ImportError:
    print("You have to install moviepy, imageio and ffmeg")
    print("pip install moviepy==0.2.2.11")
    print("pip install imageio==1.6")

import os
import threading
import time
from ipywidgets import Button, Output, IntProgress
from itertools import tee


@dataclass
class MovieMaker:
    """

    Parameters
    ----------
    view : NGLWidget
    prefix : str, default 'movie'
        prefix name of rendered image.
    output : str, default 'my_movie.gif'
        output filename of the movie.
        if output has '.gif', call `write_gif`, otherwise calling `write_videofile`
    fps : 8
        frame per second
    start, stop, step : int, default (0, -1, 1)
        how many frames you want to render.
    timeout : a number (second), default 0.1
        The waiting time between rendering two consecutive frames.
        This option should be only used with `perframe_hook` option.
    render_params : dict or None, default None
        NGL rendering params. see NGLWidget.download_image.
        If None, use default values
    moviepy_params : dict or None, default None
        moviepy params for `write_gif` method.
        if None, use default values
    in_memory : bool, default False
        if False, save rendered images to disk first
        if True, keep all image data in memory (good for small video)
    perframe_hook : callable with `view` as a single argument, default None
        if given, update the `view` by `perframe_hook`.

    Examples
    --------
    >>> import nglview as nv
    >>> import pytraj as pt
    >>> traj = pt.load(nv.datafiles.XTC, top=nv.datafiles.PDB)
    >>> view = nv.show_pytraj(traj)
    >>> from nglview.contrib.movie import MovieMaker
    >>> output = 'my.gif'
    >>> mov = MovieMaker(view, output=output)
    >>> mov.make()

    >>> # write avi format
    >>> from nglview.contrib.movie import MovieMaker

    >>> moviepy_params = {
    ...     'codec': 'mpeg4'
    ... }
    >>> movie = MovieMaker(view, output='my.avi', in_memory=True, moviepy_params=moviepy_params)
    >>> movie.make()

    Notes
    -----
    unstable API. Currently supports .gif format
    If you are using remote notebook, make sure to set in_memory=True

    Requires
    --------

    moviepy. e.g:

        pip install moviepy
        conda install freeimage

    """
    view: Any
    prefix: str = 'movie'
    output: str = 'my_movie.gif'
    fps: int = 8
    start: int = 0
    stop: int = -1
    step: int = 1
    timeout: float = 0.1
    in_memory: bool = False
    perframe_hook: Optional[Callable] = None
    render_params: Dict = field(default_factory=lambda: dict(factor=4, antialias=True, trim=False, transparent=False))
    moviepy_params: Dict = field(default_factory=dict)

    def __post_init__(self):
        if self.stop < 0:
            self.stop = self.view.max_frame + 1
        self._time_range = range(self.start, self.stop, self.step)
        self._iframe = iter(self._time_range)
        self._progress = IntProgress(max=len(self._time_range) - 1)
        self._woutput = Output()
        self._event = threading.Event()
        self._thread = None
        self._image_array = []

    def sleep(self):
        time.sleep(self.timeout)

    def _set_frame_and_hook(self, frame):
        """Set frame and hook for the view."""
        self.view.frame = frame
        time.sleep(self.timeout)
        self.perframe_hook(self.view)
        time.sleep(self.timeout)
        self.view._set_coordinates(
            frame, movie_making=True, render_params=self.render_params)

    def _handle_msg(self, widget, msg, buffers, image_array, iframe, movie, keep_data):
        """Handle messages from the view."""
        if msg['type'] == 'movie_image_data':
            image_array.append(msg.get('data'))
            try:
                frame = next(iframe)
                self.perframe_hook and self._set_frame_and_hook(frame)
                self._progress.value = frame
            except StopIteration:
                if movie:
                    self._progress.description = 'Making movie...'
                    with self._woutput:
                        # suppress moviepy's log
                        self._make_from_array(image_array)
                    if not os.path.exists(self.output):
                        self._progress.description = "ERROR: Check the maker's log"
                    else:
                        self._progress.description = 'Done'
                self._remove_on_msg()
                if keep_data:
                    self._image_array = image_array

    def make(self, movie=True, keep_data=False):
        """
        Parameters
        ----------
        keep_data: bool
            if True, save the image data in self._image_array
        movie: bool
            if True, make the movie
            else, only do the rendering (make sure keep_data=True in this case)
        """
        self._woutput.clear_output()
        image_array = []
        iframe = iter(self._time_range)
        frame = next(iframe)

        # trigger movie making communication between backend and frontend
        self.perframe_hook and self._set_frame_and_hook(frame)
        self._progress.description = 'Rendering ...'

        self._on_msg = lambda widget, msg, buffers: self._handle_msg(widget, msg, buffers, image_array, iframe, movie, keep_data)
        # FIXME: if exception happens, the on_msg callback will be never removed
        # from `self.view`
        self.view.on_msg(self._on_msg)
        return self._progress

    def _remove_on_msg(self):
        self.view.on_msg(self._on_msg, remove=True)

    def _make_from_array(self, image_array: List[str]):
        image_files = [self._base64_to_ndarray(a) for a in image_array]
        clip = mpy.ImageSequenceClip(image_files, fps=self.fps)
        with self._woutput:
            if self.output.endswith('.gif'):
                clip.write_gif(
                    self.output,
                    fps=self.fps,
                    verbose=False,
                    **self.moviepy_params)
            else:

                clip.write_videofile(
                    self.output, fps=self.fps, **self.moviepy_params)

    def interupt(self):
        """ Stop making process """
        if self._event is not None:
            self._event.set()

    @classmethod
    def _base64_to_ndarray(cls, value):
        import io
        import base64
        from PIL import Image
        import numpy as np
        im_bytes = base64.b64decode(value)
        im_bytes = io.BytesIO(im_bytes)
        # convert to numpy RGB value (for moviepy.editor.VideoClip)
        return np.array(Image.open(im_bytes))

