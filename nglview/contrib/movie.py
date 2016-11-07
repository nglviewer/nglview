import os
import time
import moviepy.editor as mpy
import threading

class MovieMaker(object):
    """

    Parameters
    ----------
    view : NGLWidget
    download_folder : str or None
        Folder that stores images. You can not arbitarily set this folder. It must be
        the download directory of the web browser you are using.
        If None, $HOME/Downloads/ will be used.
    prefix : str, default 'movie'
        prefix name of rendered image.
        Note that if calling `make(in_memory=True)`, there won't be any images on disk
    output : str, default 'my_movie.gif'
        output filename of the movie.
    fps : 8
        frame per second
    start, stop, step : int, default (0, -1, 1)
        how many frames you want to render.
    skip_render : bool, default False
        if True, do not render any frame and uses existings images in `download_folder`
        for movie making.
        if False, perform rendering first.
    render_params : dict or None, default None
        NGL rendering params. see NGLWidget.download_image.
        If None, use default values
    mpy_params : dict or None, default None
        moviepy params for `write_gif` method.
        if None, use default values

    Examples
    --------
    >>> import nglview as nv
    >>> import pytraj as pt
    >>> traj = pt.load(nv.datafiles.XTC, top=nv.datafiles.PDB)
    >>> view = nv.show_pytraj(traj)
    >>> from nglview.contrib.movie import MovieMaker
    >>> download_folder = '/Users/xxx/Downloads'
    >>> output = 'my.gif'
    >>> mov = MovieMaker(view, download_folder=download_folder, output=output)
    >>> mov.make()

    Notes
    -----
    unstable API. Currently supports .gif format
    Do not work with remote cluster since the rendered images will be downloaded to your local
    computer.

    Requires
    --------
    
    moviepy. e.g:
    
        pip install moviepy
        conda install freeimage

    """
    
    def __init__(self,
                 view,
                 download_folder=None,
                 prefix='movie',
                 output='my_movie.gif',
                 fps=8,
                 start=0,
                 stop=-1,
                 step=1,
                 skip_render=False,
                 timeout=1.,
                 render_params=None,
                 mpy_params=None):
        if download_folder is None:
            download_folder = os.getenv('HOME', '') + '/Downloads/'
        self.view = view
        self.skip_render = skip_render
        self.prefix = prefix
        self.download_folder = download_folder
        self.timeout = timeout
        self.fps = fps
        self.render_params = render_params if render_params is not None else {}
        self.mpy_params = mpy_params if mpy_params is not None else {}
        if self.render_params is not None:
            assert isinstance(self.render_params, dict)
        if self.mpy_params is not None:
            assert isinstance(self.mpy_params, dict)
        assert os.path.exists(download_folder), '{} must exists'.format(download_folder)
        self.output = output
        if stop < 0:
            stop = self.view.count
        self._time_range = range(start, stop, step)
        self._event = threading.Event()
        self._thread = None
        self._image_array = []

    def make(self, in_memory=False):
        # TODO : make base class so we can reuse this with sandbox/base.py
        self._event = threading.Event()

        def _make(event):
            if not self.skip_render:
                for i in self._time_range:
                    if not event.is_set():
                        self.view.frame = i
                        time.sleep(self.timeout)
                        if not in_memory:
                            self.view.download_image(self.prefix + '.' + str(i) + '.png',
                                    **self.render_params)
                        else:
                            self.view.render_image(**self.render_params)
                        time.sleep(self.timeout)
                        if in_memory:
                            rgb = self._base64_to_ndarray(self.view._image_data)
                            self._image_array.append(rgb)
                if not in_memory:
                    template = "{}/{}.{}.png"
                    image_files = [image_dir for image_dir in 
                                      (template.format(self.download_folder,
                                                               self.prefix,
                                                               str(i))
                                       for i in self._time_range)
                                   if os.path.exists(image_dir)]
                else:
                    image_files = self._image_array
            if not self._event.is_set():
                clip = mpy.ImageSequenceClip(image_files, fps=self.fps)
                clip.write_gif(self.output, fps=self.fps, **self.mpy_params)
        self.thread = threading.Thread(target=_make, args=(self._event,))
        self.thread.daemon = True
        self.thread.start()

    def interupt(self):
        """ Stop making process """
        if self._event is not None:
            self._event.set()

    def _base64_to_ndarray(self, value):
        import io
        import base64
        import matplotlib.image as mpimg
        im_bytes = base64.b64decode(value)
        im_bytes = io.BytesIO(im_bytes)
        # convert to numpy RGB value (for moviepy.editor.VideoClip)
        return mpimg.imread(im_bytes, format='JPG')
