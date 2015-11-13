
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

try:
    from urllib.request import urlopen
except ImportError:
    from urllib2 import urlopen


class Structure():
    def __init__( self, ext="pdb" ):
        self.ext = ext

    def get_string( self ):
        raise NotImplementedError()


class FileStructure(Structure):
    def __init__( self, path, ext="pdb" ):
        self.path = path
        self.ext = ext

    def get_string( self ):
        with open(self.path, "r") as f:
            return f.read()


class PdbIdStructure(Structure):
    def __init__( self, pdbid ):
        self.pdbid = pdbid
        self.ext = "cif"

    def get_string( self ):
        url = "http://www.rcsb.org/pdb/files/" + self.pdbid + ".cif"
        return urlopen( url ).read()


class Trajectory():
    def __init__( self ):
        pass

    def get_coordinates( self, index ):
        # [ 1,1,1, 2,2,2 ]
        raise NotImplementedError()


class SimpletrajTrajectory(Trajectory):
    def __init__( self, path ):
        import simpletraj
        self.traj_cache = simpletraj.trajectory.TrajectoryCache()
        self.path = path

    def get_coordinates( self, index ):
        traj = self.traj_cache.get( os.path.abspath( self.path ) )
        frame = traj.get_frame( int( index ) )
        return frame[ "coords" ].flatten().tolist()


class NGLWidget(widgets.DOMWidget):
    _view_name = Unicode("NGLView", sync=True)
    _view_module = Unicode("nbextensions/nglview/widget_ngl", sync=True)
    selection = Unicode("*", sync=True)
    structure = Dict(sync=True)
    representations = List(sync=True)
    coordinates = List(sync=True)

    def __init__( self, structure=None, trajectory=None, representations=None, **kwargs ):
        super(NGLWidget, self).__init__(**kwargs)
        if structure:
            self.set_structure( structure )
        if trajectory:
            self.trajectory = trajectory
        if representations:
            self.representations = representations
        else:
            self.representations = [
                { "type": "cartoon", "params": {
                    "sele": "protein"
                } },
                { "type": "ball+stick", "params": {
                    "sele": "hetero"
                } }
            ]

    def set_representations( self, representations ):
        self.representations = representations

    def set_structure( self, structure ):
        self.structure = {
            "data": structure.get_string(),
            "ext": structure.ext
        }

    def set_frame( self, index ):
        if self.trajectory:
            coordinates = self.trajectory.get_coordinates( index )
            self.coordinates = coordinates
        else:
            print "no trajectory available"


staticdir = resource_filename('nglview', os.path.join('html', 'static'))
install_nbextension(staticdir, destination='nglview', user=True)
