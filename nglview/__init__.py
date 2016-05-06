
from __future__ import print_function, absolute_import

from . import datafiles
from .utils import seq_to_string, string_types

import os
import os.path
import warnings
import tempfile
import ipywidgets as widgets
from traitlets import Unicode, Bool, Dict, List, Int, Float, Any, Bytes, observe

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
    if 'ext' in kwargs:
        extension = kwargs.pop('ext')
    else:
        extension = os.path.splitext(path)[1][1:]
    structure = FileStructure(path, ext=extension)
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
        try:
            from cStringIO import StringIO
        except ImportError:
            from io import StringIO
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
    cache = Bool().tag(sync=True)
    loaded = Bool(False).tag(sync=True)
    frame = Int().tag(sync=True)
    count = Int().tag(sync=True)
    _init_representations = List().tag(sync=True)
    structure_list = List().tag(sync=True)
    parameters = Dict().tag(sync=True)
    coordinates_dict = Dict().tag(sync=True)
    picked = Dict().tag(sync=True)
    _coordinate_dict2 = Dict().tag(sync=False)
    camera_str = Unicode().tag(sync=True)
    orientation = List().tag(sync=True)

    displayed = False
    _ngl_msg = None

    def __init__(self, structure, trajectory=None,
                 representations=None, parameters=None, **kwargs):
        try:
            self.cache = kwargs.pop('cache')
        except KeyError:
            self.cache = False
        super(NGLWidget, self).__init__(**kwargs)
        self.trajlist = []

        if parameters:
            self.parameters = parameters
        self.set_structure(structure)
        if trajectory:
            self.trajectory = trajectory
        elif hasattr(structure, "get_coordinates_dict"):
            self.trajlist = [structure,]
        elif isinstance(structure, (list, tuple)):
            self.trajlist = structure
        if self.trajlist:
            self.count = max(traj.n_frames for traj in self.trajlist if hasattr(traj,
            'n_frames'))

        # use _init_representations so we can view representations right after view is made.
        # self.representations is only have effect if we already call `view`

        if representations:
            self._ini_representations = representations
        else:
            self._init_representations = [
                {"type": "cartoon", "params": {
                    "sele": "polymer"
                }},
                {"type": "ball+stick", "params": {
                    "sele": "hetero OR mol"
                }}
            ]

        # keep track but making copy
        self._representations = self._init_representations[:]
        self._add_repr_method_shortcut()

        # do not use _displayed_callbacks since there is another Widget._display_callbacks
        self._ngl_displayed_callbacks = []

        # register to get data from JS side
        self.on_msg(self._ngl_handle_msg)

    @observe('loaded')
    def on_loaded(self, change):
        [callback(self) for callback in self._ngl_displayed_callbacks]

        if self.trajlist:
            self._set_coordinates(0)

    def _ipython_display_(self, **kwargs):
        super(NGLWidget, self)._ipython_display_(**kwargs)
        self.displayed = True

    @property
    def coordinates(self):
        if self.cache:
            return
        else:
            clist = []
            for index, traj in enumerate(self.trajlist):
                data = self._coordinate_dict2[index]['data']
                dtype = self._coordinate_dict2[index]['dtype']
                shape = self._coordinate_dict2[index]['shape']
                clist.append(decode_base64(data, dtype=dtype, shape=shape))
            return clist

    @coordinates.setter
    def coordinates(self, arrlist):
        """return current coordinate

        Parameters
        ----------
        arr : list of 2D array (shape=(n_atoms, 3))
        """
        dtype = 'f4'

        for index, arr in enumerate(arrlist): 
            coordinates_meta = dict(data=encode_numpy(arr, dtype=dtype),
                                    dtype=dtype,
                                    shape=arr.shape)
            self._coordinate_dict2[index] = coordinates_meta
        self.send({'type': 'base64_single', 'data': self._coordinate_dict2})

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
        if self.trajlist:
            # do not use traitlets to sync. slow.
            self.cache = True
            import json
            data = json.dumps([trajectory.get_coordinates_dict() for trajectory in
                self.trajlist])
            msg = dict(type='base64',
                       cache=self.cache,
                       data=data)
            self.send(msg)
        else:
            print('does not have trajlist. skip caching') 
            self.cache = False

    def uncaching(self):
        self.cache = False

    def set_representations(self, representations):
        self.representations = representations

    def set_structure(self, structures):
        structure_list = structures if isinstance(structures, (list, tuple)) else [structures,]
        self.structure_list = [{"data": _structure.get_structure_string(),
                           "ext": _structure.ext,
                           "params": _structure.params
                           } for _structure in structure_list]

    def _set_coordinates(self, index):
        if self.trajlist and not self.cache:
            coordinate_list = []
            for trajectory in self.trajlist:
                try:
                    coordinate_list.append(trajectory.get_coordinates(index))
                except IndexError:
                    coordinate_list.append(np.empty((0), dtype='f4'))
            self.coordinates = coordinate_list
        else:
            print("no trajectory available")

    @observe('frame')
    def on_frame(self, change):
        if not self.cache:
            self._set_coordinates(self.frame)

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
        # avoid space sensitivity
        repr_type = repr_type.strip()
        # overwrite selection
        selection = seq_to_string(selection).strip()

        if 'model' in kwargs:
            model = kwargs.pop('model')
        else:
            model = 0

        for k, v in kwargs.items():
            try:
                kwargs[k] = v.strip()
            except AttributeError:
                # e.g.: opacity=0.4
                kwargs[k] = v

        d = {'params': {'sele': selection}}
        d['type'] = repr_type
        d['params'].update(kwargs)

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

    def export_image(self, factor=2,
                     antialias=True,
                     trim=False,
                     transparent=False):
        """render and download scence at current frame
        """
        onProgress = False
        self._remote_call('exportImage',
                          target='Stage',
                          args=[factor, antialias, trim, transparent, onProgress])

    def _ngl_handle_msg(self, widget, msg, buffers):
        """store message sent from Javascript.

        How? use view.on_msg(get_msg)
        """
        import json
        if isinstance(msg, string_types):
            self._ngl_msg = json.loads(msg)
        else:
            self._ngl_msg = msg

    def _load_data(self, obj, **kwargs):
        '''

        Parameters
        ----------
        obj : nglview.Structure or any object having 'get_structure_string' method or
              string buffer (open(fn).read())
        '''
        if hasattr(obj, 'get_structure_string'):
            blob = obj.get_structure_string()
            kwargs['ext'] = obj.ext
            obj_is_file = False
        else:
            obj_is_file = os.path.isfile(obj)
            # assume passing string
            blob = obj
            if 'ext' not in kwargs:
                assert obj_is_file, 'must be a filename if ext is None'

        blob_type = 'path' if obj_is_file else 'blob'
        args=[{'type': blob_type, 'data': blob}]

        self._remote_call("loadFile",
                target='Stage',
                args=args,
                kwargs=kwargs)

    def _remove_component(self, model):
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

            # all callbacks will be called right after widget is loaded
            self._ngl_displayed_callbacks.append(callback)

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
