
from __future__ import print_function

import os
import warnings
import tempfile
import ipywidgets as widgets
from traitlets import Unicode, Bool, Dict, List, Int

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


class Structure(object):

    def __init__(self, text, ext='pdb', params={}):
        self._buffer = text
        self.ext = ext
        self.params = params

    def get_structure_string(self):
        return self._buffer


def load_file(path):
    '''return a Structure
    '''
    with open(path, "r") as f:
        return Structure(f.read())

def fetch_pdb(pdbid):
    '''return a Structure
    '''
    url = "http://www.rcsb.org/pdb/files/" + pdbid + ".cif"
    return Structure(urlopen(url).read(), ext='cif')


class Trajectory(object):

    def __init__(self, xyz, topology):
        '''
        Parameters
        ----------
        xyz : ndarray, shape=(n_frames, n_atoms, 3)
        topology : Topology (Structure?) or a dictionary?

        Proposed idea
        ------------
        >>> import nglview as nv
        >>> traj = nv.fetch_pdb('1l2y')
        >>> traj = nv.load_file('my_loca_file')

        >>> # load from mdtraj
        >>> import mdtraj as md
        >>> mtraj = md.load('x.trr', 'x.gro')
        >>> from nglview import convert_topology
        >>> nv_top = convert_topology(mtraj.top)
        >>> traj = nm.Trajectory(xyz=mtraj.xyz*10, topology=nv_top)

        >>> # load from pytraj
        >>> import pytraj as pt
        >>> ptraj = pt.load('amber.nc', 'amber.prmtop')
        >>> nv_top = convert_topology(ptraj.top)
        >>> traj = nm.Trajectory(xyz=ptraj.xyz, topology=nv_top)

        >>> # create viewer
        >>> from nglview import TrajectoryViewer
        >>> viwer = TrajectoryViewer(traj=traj, representations=rep, **kwd)
        >>> viewer
        '''
        self.xyz = xyz
        self.topology = topology
        self.ext = "pdb"

    def get_coordinates(self, index):
        '''return coordinate for index-th frame, length=n_atoms*3
        '''
        return self.xyz[index].flatten().tolist()

    @property
    def n_frames(self):
        return self.xyz.shape[0]



class TrajectoryViewer(widgets.DOMWidget):

    # NGLWidget is a weird name (vs TrajectoryViewer) for general users.
    _view_name = Unicode("NGLView", sync=True)
    _view_module = Unicode("nbextensions/nglview/widget_ngl", sync=True)
    selection = Unicode("*", sync=True)
    structure = Dict(sync=True)
    representations = List(sync=True)
    coordinates = List(sync=True)
    picked = Dict(sync=True)
    frame = Int(sync=True)
    count = Int(sync=True)
    clip = Dict(sync=True)
    fog = Dict(sync=True)

    def __init__(self, trajectory, representations=None, **kwargs):
        super(TrajectoryViewer, self).__init__(**kwargs)
        self.set_structure(trajectory.topology)
        # should we consider 'structure' as a Trajectory?
        self.trajectory = trajectory
        self.count = self.trajectory.n_frames
        if representations:
            self.representations = representations
        else:
            self.representations = [
                {"type": "cartoon", "params": {
                    "sele": "polymer"
                }},
                {"type": "ball+stick", "params": {
                    "sele": "hetero OR mol"
                }}
            ]

    def set_representations(self, representations):
        self.representations = representations

    def set_structure(self, structure):
        self.structure = {
            "data": structure.get_structure_string(),
            "ext": 'pdb'
            "params": structure.params
        }

    def _set_coordinates(self, index):
        if self.trajectory:
            coordinates = self.trajectory.get_coordinates(index)
            self.coordinates = coordinates
        else:
            print("no trajectory available")

    def _frame_changed(self):
        self._set_coordinates(self.frame)


staticdir = resource_filename('nglview', os.path.join('html', 'static'))
install_nbextension(staticdir, destination='nglview', user=True, verbose=0)
