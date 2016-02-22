
from __future__ import print_function

import os
import os.path
import warnings
import tempfile
import ipywidgets as widgets
from traitlets import Unicode, Bool, Dict, List, Int, Float

from IPython.display import display, Javascript
from notebook.nbextensions import install_nbextension
from notebook.services.config import ConfigManager

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from pkg_resources import resource_filename

try:
    from urllib.request import urlopen
except ImportError:
    from urllib2 import urlopen

MODULE_DIR = os.path.split( os.path.abspath( __file__ ) )[0]
DATA_DIR = os.path.join( MODULE_DIR, "data" )
GRO_FILE = os.path.join( DATA_DIR, "md_1u19.gro" )
XTC_FILE = os.path.join( DATA_DIR, "md_1u19.xtc" )


##############
# Simple API


def show_pdbid(pdbid, **kwargs):
    '''Show PDB entry.

    Example
    -------
    >>> import nglview as nv
    >>> w = nv.show_pdbid("3pqr")
    >>> w
    '''
    structure = PdbIdStructure(pdbid)
    return NGLWidget(structure, **kwargs)


def show_structure_file(path, **kwargs):
    '''Show structure file.

    Example
    -------
    >>> import nglview as nv
    >>> w = nv.show_structure_file(nv.GRO_FILE)
    >>> w
    '''
    extension = os.path.splitext(path)[1][1:]
    structure = FileStructure(path, ext=extension)
    return NGLWidget(structure, **kwargs)


def show_simpletraj(structure_path, trajectory_path, **kwargs):
    '''Show simpletraj trajectory and structure file.

    Example
    -------
    >>> import nglview as nv
    >>> w = nv.show_simpletraj(nv.GRO_FILE, nv.XTC_FILE)
    >>> w
    '''
    extension = os.path.splitext(structure_path)[1][1:]
    structure = FileStructure(structure_path, ext=extension)
    trajectory = SimpletrajTrajectory(trajectory_path)
    return NGLWidget(structure, trajectory, **kwargs)


def show_mdtraj(mdtraj_trajectory, **kwargs):
    '''Show mdtraj trajectory.

    Example
    -------
    >>> import nglview as nv
    >>> import mdtraj as md
    >>> t = md.load(nv.XTC_FILE, top=nv.GRO_FILE)
    >>> w = nv.show_mdtraj(t)
    >>> w
    '''
    structure_trajectory = MDTrajTrajectory(mdtraj_trajectory)
    return NGLWidget(structure_trajectory, **kwargs)


def show_pytraj(pytraj_trajectory, **kwargs):
    '''Show pytraj trajectory.

    Example
    -------
    >>> import nglview as nv
    >>> import pytraj as pt
    >>> t = pt.load(nv.XTC_FILE, nv.GRO_FILE)
    >>> w = nv.show_pytraj(t)
    >>> w
    '''
    structure_trajectory = PyTrajTrajectory(pytraj_trajectory)
    return NGLWidget(structure_trajectory, **kwargs)


def show_mdanalysis(atomgroup, **kwargs):
    '''Show NGL widget with MDAnalysis AtomGroup.

    Can take either a Universe or AtomGroup as its data input.

    Example
    -------
    >>> import nglview as nv
    >>> import MDAnalysis as mda
    >>> u = mda.Universe(nv.GRO_FILE, nv.XTC_FILE)
    >>> prot = u.select_atoms('protein')
    >>> w = nv.show_mdanalysis(prot)
    >>> w
    '''
    structure_trajectory = MDAnalysisTrajectory(atomgroup)
    return NGLWidget(structure_trajectory, **kwargs)


###################
# Adaptor classes


class Structure(object):

    def __init__(self):
        self.ext = "pdb"
        self.params = {}

    def get_structure_string(self):
        raise NotImplementedError()


class FileStructure(Structure):

    def __init__(self, path, ext="pdb"):
        self.path = path
        self.ext = ext
        self.params = {}
        if not os.path.isfile(path):
            raise IOError("Not a file: " + path)

    def get_structure_string(self):
        with open(self.path, "r") as f:
            return f.read()


class PdbIdStructure(Structure):

    def __init__(self, pdbid):
        self.pdbid = pdbid
        self.ext = "cif"
        self.params = {}

    def get_structure_string(self):
        url = "http://www.rcsb.org/pdb/files/" + self.pdbid + ".cif"
        return urlopen(url).read()


class Trajectory(object):

    def __init__(self):
        pass

    def get_coordinates_list(self, index):
        # [ 1,1,1, 2,2,2 ]
        raise NotImplementedError()

    def get_frame_count(self):
        raise NotImplementedError()


class SimpletrajTrajectory(Trajectory):
    '''simpletraj adaptor.

    Example
    -------
    >>> import nglview as nv
    >>> t = nv.SimpletrajTrajectory(nv.XTC_FILE)
    >>> w = nv.NGLWidget(t)
    >>> w
    '''
    def __init__(self, path):
        try:
            import simpletraj
        except ImportError as e:
            raise "'SimpletrajTrajectory' requires the 'simpletraj' package"
        self.traj_cache = simpletraj.trajectory.TrajectoryCache()
        self.path = path
        try:
            self.traj_cache.get(os.path.abspath(self.path))
        except Exception as e:
            raise e

    def get_coordinates_list(self, index):
        traj = self.traj_cache.get(os.path.abspath(self.path))
        frame = traj.get_frame(int(index))
        return frame["coords"].flatten().tolist()

    def get_frame_count(self):
        traj = self.traj_cache.get(os.path.abspath(self.path))
        return traj.numframes


class MDTrajTrajectory(Trajectory, Structure):
    '''mdtraj adaptor.

    Example
    -------
    >>> import nglview as nv
    >>> import mdtraj as md
    >>> t = md.load(nv.XTC_FILE, nv.GRO_FILE)
    >>> w = nv.NGLWidget(t)
    >>> w
    '''
    def __init__(self, trajectory):
        self.trajectory = trajectory
        self.ext = "pdb"
        self.params = {}

    def get_coordinates_list(self, index):
        frame = self.trajectory[index].xyz * 10  # convert from nm to A
        return frame.flatten().tolist()

    def get_frame_count(self):
        return len(self.trajectory.xyz)

    def get_structure_string(self):
        fd, fname = tempfile.mkstemp()
        self.trajectory[0].save_pdb(fname)
        pdb_string = os.fdopen(fd).read()
        # os.close( fd )
        return pdb_string


class PyTrajTrajectory(Trajectory, Structure):

    def __init__(self, trajectory):
        self.trajectory = trajectory
        self.ext = "pdb"
        self.params = {}

    def get_coordinates_list(self, index):
        # use trajectory[index] to use both in-memory
        #   (via pytraj.load method)
        # and out-of-core trajectory
        #   (via pytraj.iterload method)
        frame = self.trajectory[index].xyz
        return frame.flatten().tolist()

    def get_frame_count(self):
        return self.trajectory.n_frames

    def get_structure_string(self):
        fd, fname = tempfile.mkstemp(suffix=".pdb")
        self.trajectory[:1].save(fname, format="pdb", overwrite=True, options='conect')
        # self.trajectory[:1].save(fname, format="pdb", overwrite=True)
        pdb_string = os.fdopen(fd).read()
        # os.close( fd )
        return pdb_string


class MDAnalysisTrajectory(Trajectory, Structure):
    '''MDAnalysis adaptor.

    Can take either a Universe or AtomGroup.

    Example
    -------
    >>> import nglview as nv
    >>> import MDAnalysis as mda
    >>> u = mda.Universe(nv.GRO_FILE, nv.XTC_FILE)
    >>> prot = u.select_atoms('protein')
    >>> t = nv.MDAnalysisTrajectory(prot)
    >>> w = nv.NGLWidget(t)
    >>> w
    '''
    def __init__(self, atomgroup):
        self.atomgroup = atomgroup
        self.ext = "pdb"
        self.params = {}

    def get_coordinates_list(self, index):
        self.atomgroup.universe.trajectory[index]
        frame = self.atomgroup.atoms.positions
        return frame.flatten().tolist()

    def get_frame_count(self):
        return self.atomgroup.universe.trajectory.n_frames

    def get_structure_string(self):
        import MDAnalysis as mda
        import cStringIO
        u = self.atomgroup.universe
        u.trajectory[0]
        f = mda.lib.util.NamedStream(cStringIO.StringIO(), 'tmp.pdb')
        # add PDB output to the named stream
        with mda.Writer(f, self.atomgroup.atoms.n_atoms, multiframe=False) as W:
            W.write(self.atomgroup.atoms)
        # extract from the stream
        pdb_string = f.read()
        return pdb_string


###########################
# Jupyter notebook widget


class NGLWidget(widgets.DOMWidget):
    _view_name = Unicode("NGLView", sync=True)
    _view_module = Unicode("nbextensions/nglview/widget_ngl", sync=True)
    selection = Unicode("*", sync=True)
    structure = Dict(sync=True)
    representations = List(sync=True)
    coordinates = List(sync=True)
    picked = Dict(sync=True)
    frame = Int(sync=True)
    count = Int(sync=True)
    parameters = Dict(sync=True)

    def __init__(self, structure, trajectory=None,
                 representations=None, parameters=None, **kwargs):
        super(NGLWidget, self).__init__(**kwargs)
        if parameters:
            self.parameters = parameters
        self.set_structure(structure)
        if trajectory:
            self.trajectory = trajectory
        elif hasattr(structure, "get_coordinates_list"):
            self.trajectory = structure
        if hasattr(self, "trajectory") and \
                hasattr(self.trajectory, "get_frame_count"):
            self.count = self.trajectory.get_frame_count()
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
            "ext": structure.ext,
            "params": structure.params
        }

    def _set_coordinates(self, index):
        if self.trajectory:
            coordinates = self.trajectory.get_coordinates_list(index)
            self.coordinates = coordinates
        else:
            print("no trajectory available")

    def _frame_changed(self):
        self._set_coordinates(self.frame)

    def add_representation(self, repr_type, selection='all', **kwd):
        '''Add representation.

        Parameters
        ----------
        repr_type : str
            type of representation. Please see:
            http://arose.github.io/ngl/doc/#User_manual/Usage/Molecular_representations
        selection : str, default 'all'
            atom selection
        **kwd: additional arguments for representation

        Example
        -------
        >>> import nglview as nv
        >>> import pytraj as pt
        >>> t = (pt.datafiles.load_dpdp()[:].superpose('@CA'))
        >>> w = nv.show_pytraj(t)
        >>> w.add_representation('cartoon', selection='protein', color='blue')
        >>> w
        '''
        # avoid space sensitivity
        repr_type = repr_type.strip()
        selection = selection.strip()

        for k, v in kwd.items():
            kwd[k] = v.strip()

        rep = self.representations[:]
        d = {'params': {'sele': selection}}
        d['type'] = repr_type
        d['params'].update(kwd)
        rep.append(d)
        # reassign representation to trigger change
        self.representations = rep


def install(user=True, symlink=False):
    """Install the widget nbextension.

    Parameters
    ----------
    user: bool
        Install for current user instead of system-wide.
    symlink: bool
        Symlink instead of copy (for development).
    """
    staticdir = resource_filename('nglview', os.path.join('html', 'static'))
    install_nbextension(staticdir, destination='nglview',
                        user=user, symlink=symlink)

    cm = ConfigManager()
    cm.update('notebook', {
        "load_extensions": {
            "widgets/notebook/js/extension": True,
        }
    })

install(symlink=True)
