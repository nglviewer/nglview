
from __future__ import print_function, absolute_import

from . import datafiles
from .utils import seq_to_string, string_types, _camelize, _camelize_dict
from .utils import FileManager

import os
import os.path
import uuid
import warnings
import tempfile
import ipywidgets as widgets
from traitlets import Unicode, Bool, Dict, List, Int, Float, Any, Bytes, observe
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

    Examples
    --------
    >>> import nglview as nv
    >>> w = nv.show_pdbid("3pqr")
    >>> w
    '''
    structure = PdbIdStructure(pdbid)
    return NGLWidget(structure, **kwargs)

def show_rdkit(rdkit_mol, **kwargs):
    '''Show rdkit

    Examples
    --------
    >>> import nglview as nv
    >>> w = nv.show_rdkit(rdkit_mol)
    >>> w
    '''
    structure = RdkitStructure(rdkit_mol)
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

class RdkitStructure(Structure):

    def __init__(self, rdkit_mol2, ext="mol2"):
        super(RdkitStructure, self).__init__()
        self.path = ''
        self.ext = ext
        self.params = {}
        self._data = rdkit_mol2

    def get_structure_string(self):
        return StringIO(self._data).read()

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
        pass

    def get_coordinates_dict(self):
        raise NotImplementedError()

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

    def get_coordinates_dict(self):
        traj = self.traj_cache.get(os.path.abspath(self.path))

        coordinates_dict = {}
        for i in range(self.n_frames):
            frame = traj.get_frame(i)
            coordinates_dict[i] = encode_numpy(frame['coords'])
        return coordinates_dict

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

    def get_coordinates_dict(self):
        return dict((index, encode_numpy(xyz*10))
                    for index, xyz in enumerate(self.trajectory.xyz))

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
        self.id = str(uuid.uuid4())

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
        # only write 1st model
        self.trajectory.save(
            fname, overwrite=True,
            coordinates=self.trajectory.coordinates)
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
    cache = Bool().tag(sync=True)
    loaded = Bool(False).tag(sync=True)
    _finish_caching = Bool(False).tag(sync=True)
    frame = Int().tag(sync=True)
    # hack to always display movie
    count = Int(1).tag(sync=True)
    _init_representations = List().tag(sync=True)
    _init_structure_list = List().tag(sync=True)
    parameters = Dict().tag(sync=True)
    coordinates_dict = Dict().tag(sync=True)
    picked = Dict().tag(sync=True)
    _coordinate_dict2 = Dict().tag(sync=False)
    camera_str = Unicode().tag(sync=True)
    orientation = List().tag(sync=True)

    displayed = False
    _ngl_msg = None

    def __init__(self, structure=None, representations=None, parameters=None, **kwargs):
        try:
            self.cache = kwargs.pop('cache')
        except KeyError:
            self.cache = False
        super(NGLWidget, self).__init__(**kwargs)

        # do not use _displayed_callbacks since there is another Widget._display_callbacks
        self._ngl_displayed_callbacks = []
        self._add_repr_method_shortcut()

        # register to get data from JS side
        self.on_msg(self._ngl_handle_msg)

        self._trajlist = []
        self._ngl_model_ids = []
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

    def _update_count(self):
         self.count = max(traj.n_frames for traj in self._trajlist if hasattr(traj,
                         'n_frames'))
    @observe('_finish_caching')
    def on_finish_caching(self, change):
        if self._finish_caching:
            print('finish caching. Enjoy')

    @observe('loaded')
    def on_loaded(self, change):
        if change['new']:
            [callback(self) for callback in self._ngl_displayed_callbacks]

    def _ipython_display_(self, **kwargs):
        super(NGLWidget, self)._ipython_display_(**kwargs)
        self.displayed = True

    @property
    def representations(self):
        '''return list of dict
        '''
        return self._representations

    @representations.setter
    def representations(self, params_list):
        '''

        Parameters
        ----------
        params_list : list of dict
        '''

        if params_list is not self.representations:
            assert isinstance(params_list, list), 'must provide list of dict'

            if not params_list:
                for index in range(10):
                    self._clear_repr(model=index)
            else:
                for index, params in enumerate(params_list):
                    assert isinstance(params, dict), 'params must be a dict'
                    kwargs = params['params']
                    kwargs.update({'component_index': index})
                    self._representations.append(params)
                    self._remote_call('addRepresentation',
                                      target='compList',
                                      args=[params['type'],],
                                      kwargs=kwargs)


    def _add_repr_method_shortcut(self):
        # dynamically add method for NGLWidget
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

        def make_func(rep):
            """return a new function object
            """
            def func(this, selection='all', **kwargs):
                """
                """
                self.add_representation(repr_type=rep[1], selection=selection, **kwargs)
            return func

        for rep in repr_names:
            func = make_func(rep)
            fn = 'add_' + rep[0]
            from types import MethodType
            setattr(self, fn, MethodType(func, self))

    def caching(self):
        """sending all coordinates to Javascript's side. Caching makes trajectory play smoother
        but doubling your memory. If you using cache, you can not update coordinates.
        Use `view.uncaching()` then update your coordinates, then `view.caching()` again.

        This method is experimental and its name can be changed.
        """
        if self._trajlist:
            # do not use traitlets to sync. slow.
            self.cache = True
            import json
            data = json.dumps([trajectory.get_coordinates_dict() for trajectory in
                self._trajlist])
            msg = dict(type='base64',
                       cache=self.cache,
                       data=data)
            self.send(msg)
        else:
            print('does not have trajlist. skip caching') 
            self.cache = False
            self._finish_caching = False

    def uncaching(self):
        self.cache = False
        self._finish_caching = False

    def set_representations(self, representations):
        self.representations = representations

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
            if not self.cache or (self.cache and not self._finish_caching):
                coordinate_dict = {}
                for trajectory in self._trajlist:
                    traj_index = self._ngl_model_ids.index(trajectory.id)

                    try:
                        coordinate_dict[traj_index] = trajectory.get_coordinates(index)
                    except (IndexError, ValueError):
                        coordinate_dict[traj_index] = np.empty((0), dtype='f4')

                dtype = 'f4'

                # reset
                self._coordinate_dict2 = dict()
                for index, arr in coordinate_dict.items():
                    coordinates_meta = dict(data=encode_numpy(arr, dtype=dtype),
                                            dtype=dtype,
                                            shape=arr.shape)
                    self._coordinate_dict2[index] = coordinates_meta
                self.send({'type': 'base64_single', 'data': self._coordinate_dict2})
        else:
            print("no trajectory available")

    @observe('frame')
    def on_frame(self, change):
        if not self.cache or (self.cache and not self._finish_caching):
            self._set_coordinates(self.frame)

    def clear_representations(self, model=0):
        '''clear all representations for given model

        Parameters
        ----------
        model : int, default 0 (first model)
            You need to keep track how many models you added.
        '''
        self._clear_repr(model=model)

    def _clear_repr(self, model=0):
        self._remote_call("clearRepresentations",
                target='compList',
                kwargs={'component_index': model})

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

        if 'model' in kwargs2:
            model = kwargs2.pop('model')
        else:
            model = 0

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
        params.update({'component_index': model})
        self._remote_call('addRepresentation',
                          target='compList',
                          args=[d['type'],],
                          kwargs=params)


    def center_view(self, zoom=True, selection='*', model=0):
        """center view

        Examples
        --------
        view.center_view(selection='1-4')
        """
        self._remote_call('centerView', target='compList',
                          args=[zoom, selection],
                          kwargs={'component_index': model})

    @observe('_image_data')
    def get_image(self, change=""):
        '''get rendered image. Make sure to call `export_image` first

        Notes
        -----
        method name might be changed
        '''
        image = widget_image.Image()
        image._b64value = self._image_data
        return image

    def export_image(self, frame=None,
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
            view.export_image()
            
            # make sure to call `get_image` method
            view.get_image()

        Notes
        -----
        You need to call `export_image` and `get_image` in different notebook's Cells
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
        self._ngl_model_ids.append(structure.id)
        self.center_view(model=len(self._ngl_model_ids)-1)

    def add_trajectory(self, trajectory, **kwargs):
        '''

        Parameters
        ----------
        trajectory: nglview.Trajectory or derived class

        Notes
        -----
        If you combine both Structure and Trajectory, make sure
        to load all trajectories first.
        '''
        if self.loaded:
            self._load_data(trajectory, **kwargs)
        else:
            # update via structure_list
            self._init_structures.append(trajectory)
        self._trajlist.append(trajectory)
        self._update_count()
        self._ngl_model_ids.append(trajectory.id)

    def add_model(self, filename, **kwargs):
        '''add model from file/trajectory/struture

        Parameters
        ----------
        filename : str or Trajectory or Structure or their derived class
        **kwargs : additional arguments, optional

        Notes
        -----
        `add_model` should be always called after Widget is loaded

        Examples
        --------
        >>> view = nglview.Widget()
        >>> view
        >>> view.add_model(filename)
        '''
        self._load_data(filename, **kwargs)
        # assign an ID
        self._ngl_model_ids.append(str(uuid.uuid4()))

    def _load_data(self, obj, **kwargs):
        '''

        Parameters
        ----------
        obj : nglview.Structure or any object having 'get_structure_string' method or
              string buffer (open(fn).read())
        '''
        kwargs2 = _camelize_dict(kwargs)

        if 'defaultRepresentation' not in kwargs2:
            kwargs2['defaultRepresentation'] = True

        if hasattr(obj, 'get_structure_string'):
            blob = obj.get_structure_string()
            kwargs2['ext'] = obj.ext
            passing_buffer = True
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

        blob_type = 'blob' if passing_buffer else 'path'
        args=[{'type': blob_type, 'data': blob}]

        self._remote_call("loadFile",
                target='Stage',
                args=args,
                kwargs=kwargs2)

    def remove_model(self, model_id):
        """remove model by its uuid

        Examples
        --------
        >>> view.add_trajectory(traj0)
        >>> view.add_trajectory(traj1)
        >>> view.add_struture(structure)
        >>> # remove last component
        >>> view.remove_model(view._ngl_model_ids[-1])
        """
        if self._trajlist:
            for traj in self._trajlist:
                if traj.id == model_id:
                    self._trajlist.remove(traj)
        model_index = self._ngl_model_ids.index(model_id)
        self._ngl_model_ids.remove(model_id)

        self._remove_component(model=model_index)

    def _remove_component(self, model):
        """tell NGL.Stage to remove component from Stage.compList
        """
        self._remote_call('removeComponent',
                target='Stage',
                args=[model,])
        
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

    def _js_console(self):
        self.send(dict(type='get', data='any'))
        

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
