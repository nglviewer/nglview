import os
import time
import moviepy.editor as mpy

class MovieMaking(object):
    """

    Parameters
    ----------
    view : NGLWidget
    folder_dir : str
        Folder that stores images. You can not arbitarily set this folder. It must be
        the download directory of the web browser you are using.
        Normally this will be $HOME/Downloads/
    prefix : str, default 'movie'
        prefix name of rendered image
    output : str, default 'my_movie.gif'
        output filename of the movie.
    fps : 8
        frame per second
    start, stop, step : int, default (0, -1, 1)
        how many frames you want to render.
    skip_render : bool, default False
        if True, do not render any frame and uses existings images in `folder_dir`
        for movie making.
        if False, perform rendering first.

    Examples
    --------
    >>> import nglview as nv
    >>> import pytraj as pt
    >>> traj = pt.load(nv.datafiles.XTC, top=nv.datafiles.PDB)
    >>> view = nv.show_pytraj(traj)
    >>> from nglview.contrib.movie import MovieMaking
    >>> folder_dir = '/Users/xxx/Downloads'
    >>> output = 'my.gif'
    >>> mov = MovieMaking(view, folder_dir=folder_dir, output=output)
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
                 folder_dir,
                 prefix='movie',
                 output='my_movie.gif',
                 fps=8,
                 start=0,
                 stop=-1,
                 step=1,
                 skip_render=False,
                 timeout=1.):
        self.view = view
        self.skip_render = skip_render
        self.prefix = prefix
        self.folder_dir = folder_dir
        self.timeout = timeout
        self.fps = fps
        assert os.path.exists(folder_dir), '{} must exists'.format(folder_dir)
        self.output = output
        if stop < 0:
            stop = self.view.count
        self._time_range = range(start, stop, step)

    def _render(self):
        for i in self._time_range:
            self.view.frame = i
            time.sleep(self.timeout)
            self.view.download_image(self.prefix + '.' + str(i) + '.png')
            time.sleep(self.timeout)
    
    def make(self):
        """
    
        Parameters
        ----------
        view : nglview.NGLWidget
        """

        if not self.skip_render:
            self._render()
        template = "{}/{}.{}.png"
        image_files = [template.format(self.folder_dir,
                                           self.prefix,
                                           str(i))
                       for i in self._time_range]
        im = mpy.ImageSequenceClip(image_files, fps=self.fps)
        im.write_gif(self.output, fps=self.fps)
