
import warnings, os
import ipywidgets as widgets
from traitlets import Unicode, Bool, Dict, List

from IPython.display import display, Javascript
try:
    from notebook.nbextensions import install_nbextension
except ImportError:
    from IPython.html.nbextensions import install_nbextension

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from pkg_resources import resource_filename

from simpletraj import trajectory


TRAJ_CACHE = trajectory.TrajectoryCache()


class NGLWidget(widgets.DOMWidget):
    _view_name = Unicode("NGLView", sync=True)
    _view_module = Unicode("nbextensions/nglview/widget_ngl", sync=True)
    selection = Unicode("*", sync=True)
    structure = Dict(sync=True)
    coordinates = List(sync=True)

    def load_file( self, path, ext ):
        with open(path, "r") as f:
            self.structure = {
                "data": f.read(),
                "ext": ext
            }

    def load_pdb_file( self, path ):
        self.load_file( path, "pdb" )

    def load_gro_file( self, path ):
        self.load_file( path, "gro" )

    def set_coordinates( self, path, index ):
        traj = TRAJ_CACHE.get( os.path.abspath( path ) )
        frame = traj.get_frame( int( index ) )
        self.coordinates = frame[ "coords" ].flatten().tolist()


staticdir = resource_filename('nglview', os.path.join('html', 'static'))
install_nbextension(staticdir, destination='nglview', user=True)
