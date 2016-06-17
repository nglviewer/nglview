
from __future__ import print_function, absolute_import

from . import datafiles
from .utils import seq_to_string, string_types, _camelize, _camelize_dict
from .utils import FileManager
from .player import TrajectoryPlayer
import time

import os
import os.path
import uuid
import warnings
import tempfile
import ipywidgets as widgets
from traitlets import (Unicode, Bool, Dict, List, Int, Float, Any, Bytes, observe,
                       CaselessStrEnum)
from ipywidgets import widget_image

try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO

from IPython.display import display, Javascript
from notebook.nbextensions import install_nbextension
from notebook.services.config import ConfigManager

import numpy as np

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from pkg_resources import resource_filename

try:
    from urllib.request import urlopen
except ImportError:
    from urllib2 import urlopen


import base64

def encode_base64(arr, dtype='f4'):
    arr = arr.astype(dtype)
    return base64.b64encode(arr.data).decode('utf8')

def decode_base64(data, shape, dtype='f4'):
    import numpy as np
    decoded_str = base64.b64decode(data)
    return np.frombuffer(decoded_str, dtype=dtype).reshape(shape)

def _add_repr_method_shortcut(self, other):
    from types import MethodType

    repr_names  = [
            ('point', 'point'),
            ('line', 'line'),
            ('rope', 'rope'),
            ('tube', 'tube'),
            ('trace', 'trace'),
            ('label', 'label'),
            ('unitcell', 'unitcell'),
            ('cartoon', 'cartoon'),
            ('licorice', 'licorice'),
            ('ribbon', 'ribbon'),
            ('surface', 'surface'),
            ('backbone', 'backbone'),
            ('contact', 'contact'),
            ('hyperball', 'hyperball'),
            ('rocket', 'rocket'),
            ('helixorient', 'helixorient'),
            ('simplified_base', 'base'),
            ('ball_and_stick', 'ball+stick'),
            ]

    def make_func_add(rep):
        """return a new function object
        """
        def func(this, selection='all', **kwargs):
            """
            """
            self.add_representation(repr_type=rep[1], selection=selection, **kwargs)
        return func

    def make_func_remove(rep):
        """return a new function object
        """
        def func(this, **kwargs):
            """
            """
            self._remove_representations_by_name(repr_name=rep[1], **kwargs)
        return func

    for rep in repr_names:
        func_add = make_func_add(rep)
        fn_add = 'add_' + rep[0]

        func_remove = make_func_remove(rep)
        fn_remove = '_remove_' + rep[0]

        setattr(self, fn_add, MethodType(func_add, other))
        setattr(self, fn_remove, MethodType(func_remove, other))


##############
# Simple API

def show_pdbid(pdbid, **kwargs):
    '''Show PDB entry.

    Examples
    --------
    >>> import nglview as nv
    >>> w = nv.show_pdbid("3pqr")
    >>> w
    '''
    structure = PdbIdStructure(pdbid)
    return NGLWidget(structure, **kwargs)

def show_url(url, **kwargs):
    kwargs2 = dict((k, v) for k, v in kwargs.items())
    kwargs2['url'] = True
    view = NGLWidget()
    view.add_component(url, **kwargs2)
    return view

def show_text(text, **kwargs):
    """for development
    """
    structure = TextStructure(text)
    return NGLWidget(structure, **kwargs)

def show_structure_file(path, **kwargs):
    '''Show structure file. Allowed are text-based structure
    file formats that are by supported by NGL, including pdb,
    gro, mol2, sdf.

    Examples
    --------
    >>> import nglview as nv
    >>> w = nv.show_structure_file(nv.datafiles.GRO)
    >>> w
    '''
    structure = FileStructure(path)
    return NGLWidget(structure, **kwargs)


def show_simpletraj(traj, **kwargs):
    '''Show simpletraj trajectory and structure file.

    Examples
    --------
    >>> import nglview as nv
    >>> w = nv.show_simpletraj(nv.datafiles.GRO, nv.datafiles.XTC)
    >>> w
    '''
    return NGLWidget(traj, **kwargs)


def show_mdtraj(mdtraj_trajectory, **kwargs):
    '''Show mdtraj trajectory.

    Examples
    --------
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

    Examples
    --------
    >>> import nglview as nv
    >>> import pytraj as pt
    >>> t = pt.load(nv.datafiles.TRR, nv.datafiles.PDB)
    >>> w = nv.show_pytraj(t)
    >>> w
    '''
    trajlist = pytraj_trajectory if isinstance(pytraj_trajectory, (list, tuple)) else [pytraj_trajectory,]

    trajlist = [PyTrajTrajectory(traj) for traj in trajlist]
    return NGLWidget(trajlist, **kwargs)


def show_parmed(parmed_structure, **kwargs):
    '''Show pytraj trajectory.

    Examples
    --------
    >>> import nglview as nv
    >>> import parmed as pmd
    >>> t = pt.load_file(nv.datafiles.PDB)
    >>> w = nv.show_parmed(t)
    >>> w
    '''
    structure_trajectory = ParmEdTrajectory(parmed_structure)
    return NGLWidget(structure_trajectory, **kwargs)

def show_rdkit(rdkit_mol, **kwargs):
    '''Show rdkit's Mol.

    Parameters
    ----------
    rdkit_mol : rdkit.Chem.rdchem.Mol
    kwargs : additional keyword argument

    Examples
    --------
    >>> import nglview as nv
    >>> from rdkit import Chem
    >>> from rdkit.Chem import AllChem
    >>> m = Chem.AddHs(Chem.MolFromSmiles('COc1ccc2[C@H](O)[C@@H](COc2c1)N3CCC(O)(CC3)c4ccc(F)cc4'))
    >>> AllChem.EmbedMultipleConfs(m, useExpTorsionAnglePrefs=True, useBasicKnowledge=True)
    >>> view = nv.show_rdkit(m)
    >>> view

    >>> # add component m2
    >>> # create file-like object
    >>> fh = StringIO(Chem.MolToPDBBlock(m2))
    >>> view.add_component(fh, ext='pdb')

    >>> # load as trajectory, need to have ParmEd
    >>> view = nv.show_rdkit(m, parmed=True)
    '''
    from rdkit import Chem
    fh = StringIO(Chem.MolToPDBBlock(rdkit_mol))

    try:
        use_parmed = kwargs.pop("parmed")
    except KeyError:
        use_parmed = False

    if not use_parmed:
        view = NGLWidget()
        view.add_component(fh, ext='pdb', **kwargs)
        return view
    else:
        import parmed as pmd
        parm = pmd.load_rdkit(rdkit_mol)
        parm_nv = ParmEdTrajectory(parm)

        # set option for ParmEd
        parm_nv.only_save_1st_model = False

        # set option for NGL
        # wait for: https://github.com/arose/ngl/issues/126
        # to be fixed in NGLView
        # parm_nv.params = dict(firstModelOnly=True)
        return NGLWidget(parm_nv, **kwargs)

def show_mdanalysis(atomgroup, **kwargs):
    '''Show NGL widget with MDAnalysis AtomGroup.

    Can take either a Universe or AtomGroup as its data input.

    Examples
    --------
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
        self.id = str(uuid.uuid4())

    def get_structure_string(self):
        raise NotImplementedError()


class FileStructure(Structure):

    def __init__(self, path):
        super(FileStructure, self).__init__()
        self.fm = FileManager(path)
        self.ext = self.fm.ext
        self.params = {}
        if not self.fm.is_filename:
            raise IOError("Not a file: " + path)

    def get_structure_string(self):
        return self.fm.read(force_buffer=True)

class TextStructure(Structure):

    def __init__(self, text, ext='pdb', params={}):
        super(TextStructure, self).__init__()
        self.path = ''
        self.ext = ext
        self.params = params
        self._text = text

    def get_structure_string(self):
        return self._text

class RdkitStructure(Structure):

    def __init__(self, rdkit_mol, ext="pdb"):
        super(RdkitStructure, self).__init__()
        self.path = ''
        self.ext = ext
        self.params = {}
        self._rdkit_mol = rdkit_mol

    def get_structure_string(self):
        from rdkit import Chem
        fh = StringIO(Chem.MolToPDBBlock(self._rdkit_mol))
        return fh.read()

class PdbIdStructure(Structure):

    def __init__(self, pdbid):
        super(PdbIdStructure, self).__init__()
        self.pdbid = pdbid
        self.ext = "cif"
        self.params = {}

    def get_structure_string(self):
        url = "http://www.rcsb.org/pdb/files/" + self.pdbid + ".cif"
        return urlopen(url).read()


class Trajectory(object):

    def __init__(self):
        self.id = str(uuid.uuid4())
        self.shown = True

    def get_coordinates(self, index):
        raise NotImplementedError()

    @property
    def n_frames(self):
        raise NotImplementedError()


class SimpletrajTrajectory(Trajectory, Structure):
    '''simpletraj adaptor.

    Examples
    --------
    >>> import nglview as nv
    >>> t = nv.SimpletrajTrajectory(nv.datafiles.XTC, nv.datafiles.GRO)
    >>> w = nv.NGLWidget(t)
    >>> w
    '''
    def __init__(self, path, structure_path):
        try:
            import simpletraj
        except ImportError as e:
            raise ImportError(
                "'SimpletrajTrajectory' requires the 'simpletraj' package"
            )
        self.traj_cache = simpletraj.trajectory.TrajectoryCache()
        self.path = path
        self._structure_path = structure_path
        self.ext = os.path.splitext(structure_path)[1][1:]
        self.params = {}
        self.trajectory = None
        self.id = str(uuid.uuid4())
        try:
            self.traj_cache.get(os.path.abspath(self.path))
        except Exception as e:
            raise e

    def get_coordinates(self, index):
        traj = self.traj_cache.get(os.path.abspath(self.path))
        frame = traj.get_frame(index)
        return frame["coords"]

    def get_structure_string(self):
        return open(self._structure_path).read()

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
        self.id = str(uuid.uuid4())

    def get_coordinates(self, index):
        return 10*self.trajectory.xyz[index]

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
        self.id = str(uuid.uuid4())

    def get_coordinates(self, index):
        return self.trajectory[index].xyz

    @property
    def n_frames(self):
        return self.trajectory.n_frames

    def get_structure_string(self):
        fd, fname = tempfile.mkstemp(suffix=".pdb")
        self.trajectory[:1].save(fname, format="pdb", overwrite=True)
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
        self.id = str(uuid.uuid4())
        self.only_save_1st_model = True

    def get_coordinates(self, index):
        return self._xyz[index]

    @property
    def n_frames(self):
        return len(self._xyz)

    def get_structure_string(self):
        fd, fname = tempfile.mkstemp(suffix=".pdb")
        # only write 1st model
        if self.only_save_1st_model:
            self.trajectory.save(
                fname, overwrite=True,
                coordinates=self.trajectory.coordinates)
        else:
            self.trajectory.save(fname, overwrite=True)
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
        self.id = str(uuid.uuid4())

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
        u = self.atomgroup.universe
        u.trajectory[0]
        f = mda.lib.util.NamedStream(StringIO(), 'tmp.pdb')
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
    _image_data = Unicode().tag(sync=True)
    background = Unicode().tag(sync=True)
    loaded = Bool(False).tag(sync=True)
    frame = Int().tag(sync=True)
    # hack to always display movie
    count = Int(1).tag(sync=True)
    _n_dragged_files = Int().tag(sync=True)
    _init_representations = List().tag(sync=True)
    _init_structure_list = List().tag(sync=True)
    _parameters = Dict().tag(sync=True)
    picked = Dict().tag(sync=True)
    _coordinates_dict = Dict().tag(sync=False)
    _camera_str = CaselessStrEnum(['perspective', 'orthographic'],
        default_value='orthographic').tag(sync=True)
    orientation = List().tag(sync=True)

    displayed = False
    _ngl_msg = None
    _send_binary = Bool(True).tag(sync=False)

    def __init__(self, structure=None, representations=None, parameters=None, **kwargs):
        super(NGLWidget, self).__init__(**kwargs)


        # do not use _displayed_callbacks since there is another Widget._display_callbacks
        self._ngl_displayed_callbacks = []
        _add_repr_method_shortcut(self, self)

        # create after initilizing _ngl_displayed_callbacks
        self.player = TrajectoryPlayer(self)

        # register to get data from JS side
        self.on_msg(self._ngl_handle_msg)

        self._trajlist = []
        self._ngl_component_ids = []
        self._init_structures = []

        if parameters:
            self.parameters = parameters

        if isinstance(structure, Trajectory):
            self.add_trajectory(structure)
        elif isinstance(structure, (list, tuple)):
            trajectories = structure
            for trajectory in trajectories:
                self.add_trajectory(trajectory)
        else:
            if structure is not None:
                self.add_structure(structure)

        # initialize _init_structure_list
        # hack to trigger update on JS side
        if structure is not None:
            self._set_initial_structure(self._init_structures)

        if representations:
            self._ini_representations = representations
        else:
            self._init_representations = [
                {"type": "cartoon", "params": {
                    "sele": "polymer"
                }},
                {"type": "ball+stick", "params": {
                    "sele": "hetero OR mol"
                }},
                {"type": "ball+stick", "params": {
                    "sele": "not protein and not nucleic"
                }}
            ]

        # keep track but making copy
        if structure is not None:
            self._representations = self._init_representations[:]

    @property
    def parameters(self):
        return self._parameters

    @parameters.setter
    def parameters(self, params):
        params = _camelize_dict(params)
        self._parameters = params

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
        self._remote_call("setParameters",
                target='Stage',
                kwargs=dict(cameraType=self._camera_str))

    @observe('background')
    def _update_background_color(self, change):
        color = change['new']
        self.parameters = dict(background_color=color)

    @observe('_n_dragged_files')
    def on_update_dragged_file(self, change):
        if change['new'] - change['old'] == 1:
            self._ngl_component_ids.append(uuid.uuid4())

    def _update_count(self):
         self.count = max(traj.n_frames for traj in self._trajlist if hasattr(traj,
                         'n_frames'))

    @observe('loaded')
    def on_loaded(self, change):
        if change['new']:
            [callback(self) for callback in self._ngl_displayed_callbacks]

    def _ipython_display_(self, **kwargs):
        super(NGLWidget, self)._ipython_display_(**kwargs)
        self.displayed = True

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
        self._remote_call("setDelay", target="Widget", args=[delay,])

    def _set_spin(self, axis, angle):
        self._remote_call('setSpin',
                          target='Stage',
                          args=[axis, angle])
        
    @property
    def representations(self):
        return self._representations

    @representations.setter
    def representations(self, reps):
        self._representations = reps[:]
        for index in range(len(self._ngl_component_ids)):
            self.set_representations(reps)

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
            self._remote_call('addRepresentation',
                              target='compList',
                              args=[params['type'],],
                              kwargs=kwargs)

    def _remove_representations_by_name(self, repr_name, component=0):
        self._remote_call('removeRepresentationsByName',
                          target='Widget',
                          args=[repr_name, component])

    def _set_initial_structure(self, structures):
        """initialize structures for Widget

        Parameters
        ----------
        structures : list
            list of Structure or Trajectory
        """
        _init_structure_list = structures if isinstance(structures, (list, tuple)) else [structures,]
        self._init_structure_list = [{"data": _structure.get_structure_string(),
                                "ext": _structure.ext,
                                "params": _structure.params,
                                "id": _structure.id
                                } for _structure in _init_structure_list]

    def _set_coordinates(self, index):
        '''update coordinates for all trajectories at index-th frame
        '''
        if self._trajlist:
            coordinates_dict = {}
            for trajectory in self._trajlist:
                traj_index = self._ngl_component_ids.index(trajectory.id)

                try:
                    if trajectory.shown:
                        coordinates_dict[traj_index] = trajectory.get_coordinates(index)
                    else:
                        coordinates_dict[traj_index] = np.empty((0), dtype='f4')
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
            encoded_coordinates_dict = dict((k, encode_base64(v))
                                 for (k, v) in self._coordinates_dict.items())
            mytime = time.time() * 1000
            self.send({'type': 'base64_single', 'data': encoded_coordinates_dict,
                'mytime': mytime})
        else:
            # send binary
            buffers = []
            coordinates_meta = dict()
            for index, arr in self._coordinates_dict.items():
                buffers.append(arr.astype('f4').tobytes())
                coordinates_meta[index] = index
            mytime = time.time() * 1000
            self.send({'type': 'binary_single', 'data': coordinates_meta,
                'mytime': mytime}, buffers=buffers)

    @observe('frame')
    def on_frame(self, change):
        """set and send coordinates at current frame
        """
        self._set_coordinates(self.frame)

    def clear(self, *args, **kwargs):
        self.clear_representations(*args, **kwargs)

    def clear_representations(self, component=0):
        '''clear all representations for given component

        Parameters
        ----------
        component : int, default 0 (first model)
            You need to keep track how many components you added.
        '''
        self._clear_repr(component=component)

    def _clear_repr(self, component=0):
        self._remote_call("clearRepresentations",
                target='compList',
                kwargs={'component_index': component})

    def add_representation(self, repr_type, selection='all', **kwargs):
        '''Add representation.

        Parameters
        ----------
        repr_type : str
            type of representation. Please see:
            http://arose.github.io/ngl/doc/#User_manual/Usage/Molecular_representations
        selection : str or 1D array (atom indices) or any iterator that returns integer, default 'all'
            atom selection
        **kwargs: additional arguments for representation

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
        self._remote_call('addRepresentation',
                          target='compList',
                          args=[d['type'],],
                          kwargs=params)


    def center(self, *args, **kwargs):
        """alias of `center_view`
        """
        self.center_view(*args, **kwargs)

    def center_view(self, zoom=True, selection='*', component=0):
        """center view

        Examples
        --------
        view.center_view(selection='1-4')
        """
        self._remote_call('centerView', target='compList',
                          args=[zoom, selection],
                          kwargs={'component_index': component})

    @observe('_image_data')
    def get_image(self, change=""):
        '''get rendered image. Make sure to call `render_image` first

        Notes
        -----
        method name might be changed
        '''
        image = widget_image.Image()
        image._b64value = self._image_data
        return image

    def render_image(self, frame=None,
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
        params = dict(factor=factor,
                      antialias=antialias,
                      trim=trim,
                      transparent=transparent)
        self._remote_call('_exportImage',
                          target='Widget',
                          kwargs=params)

    def download_image(self, filename='screenshot.png',
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
        params = dict(factor=factor,
                      antialias=antialias,
                      trim=trim,
                      transparent=transparent)
        self._remote_call('_downloadImage',
                          target='Widget',
                          args=[filename,],
                          kwargs=params)

    def _ngl_handle_msg(self, widget, msg, buffers):
        """store message sent from Javascript.

        How? use view.on_msg(get_msg)
        """
        import json
        if isinstance(msg, string_types):
            self._ngl_msg = json.loads(msg)
        else:
            self._ngl_msg = msg

            msg_type = self._ngl_msg.get('type')
            if msg_type == 'request_frame':
                self.frame += self.player.step
                if self.frame >= self.count:
                    self.frame = 0
                elif self.frame < 0:
                    self.frame = self.count - 1

    def add_structure(self, structure, **kwargs):
        '''

        Parameters
        ----------
        structure : nglview.Structure object

        Notes
        -----
        If you combine both Structure and Trajectory, make sure
        to load all trajectories first.

        >>> view.add_trajectory(traj0)
        >>> view.add_trajectory(traj1)
        >>> # then add Structure
        >>> view.add_structure(...)
        '''
        if self.loaded:
            self._load_data(structure, **kwargs)
        else:
            # update via structure_list
            self._init_structures.append(structure)
        self._ngl_component_ids.append(structure.id)
        self.center_view(component=len(self._ngl_component_ids)-1)
        self._update_component_auto_completion()

    def add_trajectory(self, trajectory, **kwargs):
        '''

        Parameters
        ----------
        trajectory: nglview.Trajectory or its derived class or 
            pytraj.Trajectory-like, mdtraj.Trajectory or MDAnalysis objects

        Notes
        -----
        `add_trajectory` is just a special case of `add_component`
        '''
        backends = dict(pytraj=PyTrajTrajectory,
                       mdtraj=MDTrajTrajectory,
                       MDAnalysis=MDAnalysisTrajectory,
                       parmed=ParmEdTrajectory)

        package_name = trajectory.__module__.split('.')[0]

        if package_name in backends:
            trajectory = backends[package_name](trajectory)
        else:
            trajectory = trajectory

        if self.loaded:
            self._load_data(trajectory, **kwargs)
        else:
            # update via structure_list
            self._init_structures.append(trajectory)
        setattr(trajectory, 'shown', True)
        self._trajlist.append(trajectory)
        self._update_count()
        self._ngl_component_ids.append(trajectory.id)
        self._update_component_auto_completion()

    def add_component(self, filename, **kwargs):
        '''add component from file/trajectory/struture

        Parameters
        ----------
        filename : str or Trajectory or Structure or their derived class or url
            if you specify url, you must specify `url=True` in kwargs
        **kwargs : additional arguments, optional

        Notes
        -----
        `add_component` should be always called after Widget is loaded

        Examples
        --------
        >>> view = nglview.Widget()
        >>> view
        >>> view.add_component(filename)
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
            is_url = kwargs2.pop('url')
        except KeyError:
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
                fh = FileManager(obj,
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
            args=[{'type': blob_type, 'data': blob, 'binary': binary}]
        else:
            # is_url
            blob_type = 'url'
            url = obj
            args=[{'type': blob_type, 'data': url, 'binary': False}]

        self._remote_call("loadFile",
                target='Stage',
                args=args,
                kwargs=kwargs2)

    def remove_component(self, component_id):
        """remove component by its uuid

        Examples
        --------
        >>> view.add_trajectory(traj0)
        >>> view.add_trajectory(traj1)
        >>> view.add_struture(structure)
        >>> # remove last component
        >>> view.remove_component(view._ngl_component_ids[-1])
        """
        self._clear_component_auto_completion()
        if self._trajlist:
            for traj in self._trajlist:
                if traj.id == component_id:
                    self._trajlist.remove(traj)
        component_index = self._ngl_component_ids.index(component_id)
        self._ngl_component_ids.remove(component_id)

        self._remove_component(component=component_index)
        self._update_component_auto_completion()

    def _remove_component(self, component):
        """tell NGL.Stage to remove component from Stage.compList
        """
        self._remote_call('removeComponent',
                target='Stage',
                args=[component,])
        
    def _remote_call(self, method_name, target='Stage', args=None, kwargs=None):
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

        # perform centerView for 1-th component
        # component = Stage.compList[1];
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

        if self.displayed is True:
            self.send(msg)
        else:
            # send later
            def callback(widget, msg=msg):
                widget.send(msg)

            callback._method_name = method_name

            # all callbacks will be called right after widget is loaded
            self._ngl_displayed_callbacks.append(callback)

    @property
    def n_components(self):
        return len(self._ngl_component_ids)

    def _get_traj_by_id(self, itsid):
        """return nglview.Trajectory or its derived class object
        """
        for traj in self._trajlist:
            if traj.id == itsid:
                return traj
        return None

    def hide(self, indices):
        """set invisibility for given components (by their indices)
        """
        traj_ids = set(traj.id for traj in self._trajlist)

        for index in indices:
            assert index < self.n_components
            comp_id = self._ngl_component_ids[index]
            if comp_id in traj_ids:
                traj = self._get_traj_by_id(comp_id)
                traj.shown = False
            self._remote_call("setVisibility",
                    target='compList',
                    args=[False,],
                    kwargs={'component_index': index})

    def show(self, *args):
        self.show_only(*args)

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
                args = [True,]
                if traj is not None:
                    traj.shown = True
            else:
                args = [False,]
                if traj is not None:
                    traj.shown = False

            self._remote_call("setVisibility",
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
        assert index < len(self._ngl_component_ids)
        return ComponentViewer(self, index) 

    def __iter__(self):
        for i, _ in enumerate(self._ngl_component_ids):
            yield self[i]

    def _play(self, start=0, stop=-1, step=1, delay=0.08, n_times=1):
        '''for testing. Might be removed in the future

        Notes
        -----
        To stop, you need to choose 'Kernel' --> 'Interupt' in your notebook tab (top)
        '''
        from itertools import repeat
        from time import sleep

        if stop == -1:
            stop = self.count

        for indices in repeat(range(start, stop, step), n_times):
            for frame in indices:
                self.frame = frame
                sleep(delay)


class ComponentViewer(object):
    """Convenient attribute for NGLWidget. See example below.

    Examples
    --------
    >>> view = nv.NGLWidget()
    >>> view.add_trajectory(traj) # traj is a component 0
    >>> view.add_component(filename) # component 1
    >>> view.component_0.clear_representations()
    >>> view.component_0.add_cartoon()
    >>> view.component_1.add_licorice()
    >>> view.remove_component(view.comp1.id)
    """

    def __init__(self, view, index):
        self._view = view
        self._index = index
        _add_repr_method_shortcut(self, self._view)
        self._borrow_attribute(self._view, ['clear_representations',
                                            '_remove_representations_by_name',
                                            'center_view',
                                            'center',
                                            'clear',
                                            'set_representations'],

                                            ['get_structure_string',
                                             'get_coordinates',
                                             'n_frames'])

    @property
    def id(self):
        return self._view._ngl_component_ids[self._index]

    def hide(self):
        """set invisibility for given components (by their indices)
        """
        self._view._remote_call("setVisibility",
                target='compList',
                args=[False,],
                kwargs={'component_index': self._index})
        traj = self._view._get_traj_by_id(self.id)
        if traj is not None:
            traj.shown = False

    def show(self):
        """set invisibility for given components (by their indices)
        """
        self._view._remote_call("setVisibility",
                target='compList',
                args=[True,],
                kwargs={'component_index': self._index})

        traj = self._view._get_traj_by_id(self.id)
        if traj is not None:
            traj.shown = True

    def add_representation(self, repr_type, selection='all', **kwargs):
        kwargs['component'] = self._index
        self._view.add_representation(repr_type=repr_type, selection=selection, **kwargs)

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
        
def install(user=True, symlink=False):
    """Install the widget nbextension.

    Parameters
    ----------
    user: bool
        Install for current user instead of system-wide.
    symlink: bool
        Symlink instead of copy (for development).
    """
    staticdir = resource_filename('nglview', 'js')
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

def _get_notebook_info():
    import notebook, ipywidgets, traitlets

    print('notebook', notebook.__version__)
    print('ipywidgets', ipywidgets.__version__)
    print('traitlets', traitlets.__version__)
    print('nglview', __version__)
