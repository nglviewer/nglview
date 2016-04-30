
from __future__ import print_function, absolute_import

from . import datafiles
from .utils import seq_to_string

import os
import os.path
import warnings
import tempfile
import ipywidgets as widgets
from traitlets import Unicode, Bool, Dict, List, Int, Float, Any, Bytes

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


import base64

def encode_numpy(arr, dtype='f4'):
    arr = arr.astype(dtype)
    return base64.b64encode(arr.data).decode('utf8')

def decode_base64(data, shape, dtype='f4'):
    import numpy as np
    decoded_str = base64.b64decode(data)
    return np.frombuffer(decoded_str, dtype=dtype).reshape(shape)


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
    '''Show structure file. Allowed are text-based structure
    file formats that are by supported by NGL, including pdb,
    gro, mol2, sdf.

    Example
    -------
    >>> import nglview as nv
    >>> w = nv.show_structure_file(nv.datafiles.GRO)
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
    >>> w = nv.show_simpletraj(nv.datafiles.GRO, nv.datafiles.XTC)
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
    >>> t = md.load(nv.datafiles.XTC, top=nv.datafiles.GRO)
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
    >>> t = pt.load(nv.datafiles.TRR, nv.datafiles.PDB)
    >>> w = nv.show_pytraj(t)
    >>> w
    '''
    structure_trajectory = PyTrajTrajectory(pytraj_trajectory)
    return NGLWidget(structure_trajectory, **kwargs)


def show_parmed(parmed_structure, **kwargs):
    '''Show pytraj trajectory.

    Example
    -------
    >>> import nglview as nv
    >>> import parmed as pmd
    >>> t = pt.load_file(nv.datafiles.PDB)
    >>> w = nv.show_parmed(t)
    >>> w
    '''
    structure_trajectory = ParmEdTrajectory(parmed_structure)
    return NGLWidget(structure_trajectory, **kwargs)


def show_mdanalysis(atomgroup, **kwargs):
    '''Show NGL widget with MDAnalysis AtomGroup.

    Can take either a Universe or AtomGroup as its data input.

    Example
    -------
    >>> import nglview as nv
    >>> import MDAnalysis as mda
    >>> u = mda.Universe(nv.datafiles.GRO, nv.datafiles.XTC)
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

    def get_coordinates_dict(self):
        raise NotImplementedError()

    def get_coordinates(self, index):
        raise NotImplementedError()

    @property
    def n_frames(self):
        raise NotImplementedError()


class SimpletrajTrajectory(Trajectory):
    '''simpletraj adaptor.

    Example
    -------
    >>> import nglview as nv
    >>> t = nv.SimpletrajTrajectory(nv.datafiles.XTC)
    >>> w = nv.NGLWidget(t)
    >>> w
    '''
    def __init__(self, path):
        try:
            import simpletraj
        except ImportError as e:
            raise ImportError(
                "'SimpletrajTrajectory' requires the 'simpletraj' package"
            )
        self.traj_cache = simpletraj.trajectory.TrajectoryCache()
        self.path = path
        try:
            self.traj_cache.get(os.path.abspath(self.path))
        except Exception as e:
            raise e

    def get_coordinates(self, index):
        traj = self.traj_cache.get(os.path.abspath(self.path))
        frame = traj.get_frame(int(index))
        return frame["coords"]

    @property
    def n_frames(self):
        traj = self.traj_cache.get(os.path.abspath(self.path))
        return traj.numframes


class MDTrajTrajectory(Trajectory, Structure):
    '''mdtraj adaptor.

    Example
    -------
    >>> import nglview as nv
    >>> import mdtraj as md
    >>> traj = md.load(nv.datafiles.XTC, nv.datafiles.GRO)
    >>> t = MDTrajTrajectory(traj)
    >>> w = nv.NGLWidget(t)
    >>> w
    '''
    def __init__(self, trajectory):
        self.trajectory = trajectory
        self.ext = "pdb"
        self.params = {}

    def get_coordinates_dict(self):
        return dict((index, encode_numpy(xyz))
                    for index, xyz in enumerate(self.trajectory.xyz))

    def get_coordinates(self, index):
        return self.trajectory.xyz[index]

    @property
    def n_frames(self):
        return self.trajectory.n_frames

    def get_structure_string(self):
        fd, fname = tempfile.mkstemp()
        self.trajectory[0].save_pdb(fname)
        pdb_string = os.fdopen(fd).read()
        # os.close( fd )
        return pdb_string


class PyTrajTrajectory(Trajectory, Structure):
    '''PyTraj adaptor.

    Example
    -------
    >>> import nglview as nv
    >>> import pytraj as pt
    >>> traj = pt.load(nv.datafiles.TRR, nv.datafiles.PDB)
    >>> t = nv.PyTrajTrajectory(traj)
    >>> w = nv.NGLWidget(t)
    >>> w
    '''
    def __init__(self, trajectory):
        self.trajectory = trajectory
        self.ext = "pdb"
        self.params = {}

    def get_coordinates_dict(self):
        return dict((index, encode_numpy(xyz))
                    for index, xyz in enumerate(self.trajectory.xyz))

    def get_coordinates(self, index):
        return self.trajectory[index].xyz

    @property
    def n_frames(self):
        return self.trajectory.n_frames

    def get_structure_string(self):
        fd, fname = tempfile.mkstemp(suffix=".pdb")
        self.trajectory[:1].save(
            fname, format="pdb", overwrite=True, options='conect'
        )
        pdb_string = os.fdopen(fd).read()
        # os.close( fd )
        return pdb_string


class ParmEdTrajectory(Trajectory, Structure):
    '''ParmEd adaptor.
    '''
    def __init__(self, trajectory):
        self.trajectory = trajectory
        self.ext = "pdb"
        self.params = {}
        # only call get_coordinates once
        self._xyz = trajectory.get_coordinates()

    def get_coordinates_dict(self):
        return dict((index, encode_numpy(xyz))
                    for index, xyz in enumerate(self._xyz))

    def get_coordinates(self, index):
        return self._xyz[index]

    @property
    def n_frames(self):
        return len(self._xyz)

    def get_structure_string(self):
        fd, fname = tempfile.mkstemp(suffix=".pdb")
        self.trajectory.save(
            fname, overwrite=True)
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
    >>> u = mda.Universe(nv.datafiles.GRO, nv.datafiles.XTC)
    >>> prot = u.select_atoms('protein')
    >>> t = nv.MDAnalysisTrajectory(prot)
    >>> w = nv.NGLWidget(t)
    >>> w
    '''
    def __init__(self, atomgroup):
        self.atomgroup = atomgroup
        self.ext = "pdb"
        self.params = {}

    def get_coordinates_dict(self):

        return dict((index, encode_numpy(self.atomgroup.atoms.positions))
                    for index, _ in enumerate(self.atomgroup.universe.trajectory))

    def get_coordinates(self, index):
        self.atomgroup.universe.trajectory[index]
        xyz = self.atomgroup.atoms.positions
        return xyz

    @property
    def n_frames(self):
        return self.atomgroup.universe.trajectory.n_frames

    def get_structure_string(self):
        try:
            import MDAnalysis as mda
        except ImportError as e:
            raise ImportError(
                "'MDAnalysisTrajectory' requires the 'MDAnalysis' package"
            )
        import cStringIO
        u = self.atomgroup.universe
        u.trajectory[0]
        f = mda.lib.util.NamedStream(cStringIO.StringIO(), 'tmp.pdb')
        atoms = self.atomgroup.atoms
        # add PDB output to the named stream
        with mda.Writer(f, atoms.n_atoms, multiframe=False) as W:
            W.write(atoms)
        # extract from the stream
        pdb_string = f.read()
        return pdb_string


###########################
# Jupyter notebook widget


class NGLWidget(widgets.DOMWidget):
    _view_name = Unicode("NGLView").tag(sync=True)
    _view_module = Unicode("nbextensions/nglview/widget_ngl").tag(sync=True)
    selection = Unicode("*").tag(sync=True)
    cache = Bool().tag(sync=True)
    frame = Int().tag(sync=True)
    count = Int().tag(sync=True)
    representations = List().tag(sync=True)
    structure = Dict().tag(sync=True)
    parameters = Dict().tag(sync=True)
    _coordinates_meta = Dict().tag(sync=True)
    coordinates_dict = Dict().tag(sync=True)
    picked = Dict().tag(sync=True)

    _tmp_msg = None

    def __init__(self, structure, trajectory=None,
                 representations=None, parameters=None, **kwargs):
        try:
            self.cache = kwargs.pop('cache')
        except KeyError:
            self.cache = False
        super(NGLWidget, self).__init__(**kwargs)
        if parameters:
            self.parameters = parameters
        self.set_structure(structure)
        if trajectory:
            self.trajectory = trajectory
        elif hasattr(structure, "get_coordinates_dict"):
            self.trajectory = structure
        if hasattr(self, "trajectory") and \
                hasattr(self.trajectory, "n_frames"):
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

        self._add_repr_method_shortcut()
        self.on_msg(self._ngl_get_msg)

    @property
    def coordinates(self):
        if self.cache:
            return
        else:
            data = self._coordinates_meta['data']
            dtype = self._coordinates_meta['dtype']
            shape = self._coordinates_meta['shape']
            return decode_base64(data, dtype=dtype, shape=shape)

    @coordinates.setter
    def coordinates(self, arr):
        """return current coordinate

        Parameters
        ----------
        arr : 2D array, shape=(n_atoms, 3)
        """
        dtype = 'f4'
        coordinates_meta = dict(data=encode_numpy(arr, dtype=dtype),
                                dtype=dtype,
                                shape=arr.shape)
        self._coordinates_meta = coordinates_meta
                                      
    def _add_repr_method_shortcut(self):
        # dynamically add method for NGLWidget
        repr_names  = [
                ('point', 'point'),
                ('line', 'line'),
                ('rope', 'rope'),
                ('tube', 'tube'),
                ('trace', 'trace'),
                ('label', 'label'),
                ('cartoon', 'cartoon'),
                ('licorice', 'licorice'),
                ('ribbon', 'ribbon'),
                ('surface', 'surface'),
                ('backbone', 'backbone'),
                ('contact', 'contact'),
                ('crossing', 'crossing'),
                ('hyperball', 'hyperball'),
                ('rocket', 'rocket'),
                ('helixorient', 'helixorient'),
                ('simplified_base', 'base'),
                ('ball_and_stick', 'ball+stick'),
                ]

        def make_func(rep):
            """return a new function object
            """
            def func(this, selection='all', **kwd):
                """
                """
                self.add_representation(repr_type=rep[1], selection=selection, **kwd)
            return func

        for rep in repr_names:
            func = make_func(rep)
            fn = 'add_' + rep[0]
            from types import MethodType
            setattr(self, fn, MethodType(func, self))

    def caching(self):
        if hasattr(self.trajectory, "get_coordinates_dict"):
            # should use use coordinates_dict to sync?
            # my molecule disappear.
            # self.coordinates_dict = self.trajectory.get_coordinates_dict()

            self.cache = True
            msg = dict(type='base64',
                       cache=self.cache,
                       data=self.trajectory.get_coordinates_dict())
            self.send(msg)
        else:
            print('warning: does not have get_coordinates_dict method, turn off cache') 
            self.cache = False

    def uncaching(self):
        self.cache = False

    def set_representations(self, representations):
        self.representations = representations

    def set_structure(self, structure):
        self.structure = {
            "data": structure.get_structure_string(),
            "ext": structure.ext,
            "params": structure.params
        }

    def _set_coordinates(self, index):
        if self.trajectory and not self.cache:
            coordinates = self.trajectory.get_coordinates(index)
            self.coordinates = coordinates
        else:
            print("no trajectory available")

    def _frame_changed(self):
        if not self.cache:
            self._set_coordinates(self.frame)


    def add_representation(self, repr_type, selection='all', **kwd):
        '''Add representation.

        Parameters
        ----------
        repr_type : str
            type of representation. Please see:
            http://arose.github.io/ngl/doc/#User_manual/Usage/Molecular_representations
        selection : str or 1D array (atom indices) or any iterator that returns integer, default 'all'
            atom selection
        **kwd: additional arguments for representation

        Example
        -------
        >>> import nglview as nv
        >>> import pytraj as pt
        >>> t = (pt.datafiles.load_dpdp()[:].superpose('@CA'))
        >>> w = nv.show_pytraj(t)
        >>> w.add_representation('cartoon', selection='protein', color='blue')
        >>> w.add_representation('licorice', selection=[3, 8, 9, 11], color='red')
        >>> w
        '''
        # avoid space sensitivity
        repr_type = repr_type.strip()
        # overwrite selection
        selection = seq_to_string(selection).strip()

        for k, v in kwd.items():
            try:
                kwd[k] = v.strip()
            except AttributeError:
                # e.g.: opacity=0.4
                kwd[k] = v

        rep = self.representations[:]
        d = {'params': {'sele': selection}}
        d['type'] = repr_type
        d['params'].update(kwd)
        rep.append(d)
        # reassign representation to trigger change
        self.representations = rep

    def _ngl_get_msg(self, widget, msg, buffers):
        """store message sent from Javascript.

        How? use view.on_msg(get_msg)
        """
        import json
        self._tmp_msg = json.loads(msg)
        
    def _remote_call(self, method_name, target='stage', args=None, kwargs=None):
        """call NGL's methods from Python.
        
        Parameters
        ----------
        method_name : str
        target : str, {'stage', 'viewer', 'component'}
        args : list
        kwargs : dict
            if target is 'component', "component_index" could be passed
            to specify which component will call the method.

        Examples
        --------
        view._remote_call('loadFile', args=['1L2Y.pdb'],
                          target='stage', kwargs={'defaultRepresentation': True})

        # perform centerView for 1-th component
        # component = stage.compList[1];
        # component.centerView(true, "1-12");
        view._remote_call('centerView',
                          target='component',
                          args=[True, "1-12"],
                          kwargs={'component_index': 1})
        """
        args = [] if args is None else args
        kwargs = {} if kwargs is None else kwargs

        msg = {}

        if 'component_index' in kwargs:
            msg['component_index'] = kwargs.pop('component_index')

        msg['target'] = target
        msg['type'] = 'call_method'
        msg['methodName'] = method_name
        msg['args'] = args
        msg['kwargs'] = kwargs

        self.send(msg)

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
                        user=user, symlink=symlink,
                        verbose=0)

    cm = ConfigManager()
    cm.update('notebook', {
        "load_extensions": {
            "widgets/notebook/js/extension": True,
        }
    })

install(symlink=False)

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
