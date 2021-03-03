import os
from contextlib import contextmanager
import warnings
import os.path
import uuid
from functools import partial
from io import StringIO
from tempfile import mkstemp
from urllib.request import urlopen

import numpy as np

from . import config
from .base_adaptor import Structure, Trajectory
from .utils.py_utils import FileManager, tempfolder

__all__ = [
    'FileStructure',
    'TextStructure',
    'RdkitStructure',
    'PdbIdStructure',
    'ASEStructure',
    'BiopythonStructure',
    'IOTBXStructure',
    'IODataStructure',
    'QCElementalStructure',
    'Psi4Structure',
    'OpenbabelStructure',
    'RosettaStructure',
    'ProdyStructure',
    'SimpletrajTrajectory',
    'ProdyTrajectory',
    'MDTrajTrajectory',
    'PyTrajTrajectory',
    'ParmEdTrajectory',
    'MDAnalysisTrajectory',
    'HTMDTrajectory',
    'ASETrajectory',
    'SchrodingerStructure',
    'SchrodingerTrajectory',
    'register_backend',
]


@contextmanager
def mkstemp_wrapper(*args, **kwargs):
    # NamedTemporaryFile cannot be used here because it makes it impossible
    # on Windows to the file for writing. mkstemp is a bit less restrictive
    # in this regard.
    fd, fname = mkstemp(*args, **kwargs)
    yield fname
    # On windows, the file must be closed before it can be removed.
    os.close(fd)
    os.remove(fname)


def _get_structure_string(write_method, suffix='.pdb'):
    with mkstemp_wrapper(suffix=suffix) as fname:
        write_method(fname)
        with open(fname) as fh:
            return fh.read()


class register_backend:
    def __init__(self, package_name):
        # package_name must match exactly to your Python package
        self.package_name = package_name

    def __call__(self, cls):
        config.BACKENDS[self.package_name] = cls
        return cls


class FileStructure(Structure):
    def __init__(self, path):
        super().__init__()
        self.fm = FileManager(path)
        self.ext = self.fm.ext
        self.params = {}
        self.path = path

    def get_structure_string(self):
        return self.fm.read(force_buffer=True)


class TextStructure(Structure):
    def __init__(self, text, ext='pdb', params={}):
        super().__init__()
        self.path = ''
        self.ext = ext
        self.params = params
        self._text = text

    def get_structure_string(self):
        return self._text


@register_backend('rdkit')
class RdkitStructure(Structure):
    def __init__(self, rdkit_mol, ext="pdb", conf_id=-1):
        super().__init__()
        self.path = ''
        self.ext = ext
        self._conf_id = conf_id
        self.params = {}
        self._rdkit_mol = rdkit_mol

    def get_structure_string(self):
        from rdkit import Chem

        reader = (self.ext == "pdb") and \
                 Chem.MolToPDBBlock or Chem.MolToMolBlock
        return reader(self._rdkit_mol, confId=self._conf_id)


class PdbIdStructure(Structure):
    def __init__(self, pdbid):
        super().__init__()
        self.pdbid = pdbid
        self.ext = "cif"
        self.params = {}

    def get_structure_string(self):
        url = "http://www.rcsb.org/pdb/files/" + self.pdbid + ".cif"
        return urlopen(url).read().decode('utf-8')


class ASEStructure(Structure):
    def __init__(self, ase_atoms, ext='pdb', params={}):
        super().__init__()
        self.path = ''
        self.ext = ext
        self.params = params
        self._ase_atoms = ase_atoms

    def get_structure_string(self):
        return _get_structure_string(self._ase_atoms.write)


class IODataStructure(Structure):
    def __init__(self, obj):
        super().__init__()
        self._obj = obj

    def get_structure_string(self):
        """Require `ase` package
        """
        import ase.io
        with mkstemp_wrapper(suffix='.xyz') as fname:
            self._obj.to_file(fname)
            return ASEStructure(ase.io.read(fname)).get_structure_string()


class QCElementalStructure(Structure):
    def __init__(self, obj):
        super().__init__()
        self._obj = obj
        self.ext = 'sdf'

    def get_structure_string(self):
        return self._obj.to_string('nglview-sdf')


class Psi4Structure(QCElementalStructure):
    pass


class OpenbabelStructure(Structure):
    def __init__(self, obj):
        super().__init__()
        self._obj = obj

    def get_structure_string(self):
        try:
            # Open Babel >= 3.0.0
            from openbabel import openbabel
        except ImportError:
            import openbabel
        oc = openbabel.OBConversion()
        oc.SetOutFormat('pdb')
        write = partial(oc.WriteFile, self._obj)
        return _get_structure_string(write)


class BiopythonStructure(Structure):
    def __init__(self, entity, ext='pdb', params={}):
        super().__init__()
        self.path = ''
        self.ext = ext
        self.params = params
        self._entity = entity

    def get_structure_string(self):
        from Bio.PDB import PDBIO
        from io import StringIO
        io_pdb = PDBIO()
        io_pdb.set_structure(self._entity)
        io_str = StringIO()
        io_pdb.save(io_str)
        return io_str.getvalue()


class IOTBXStructure(Structure):
    def __init__(self, obj, ext='pdb', params={}):
        """
        obj must have as_pdb_string method
        """
        super().__init__()
        self.path = ''
        self.ext = ext
        self.params = params
        self._mol = obj

    def get_structure_string(self):
        return self._mol.as_pdb_string()


class RosettaStructure(Structure):
    def __init__(self, pose, ext='pdb', params={}):
        # type: (pyrosetta.rosetta.core.pose.Pose, str, Dict) -> None
        super().__init__()
        self.path = ''
        self.ext = ext
        self.params = params
        self._mol = pose

    def get_structure_string(self):
        return _get_structure_string(self._mol.dump_pdb)


@register_backend('prody')
class ProdyStructure(Structure):
    def __init__(self, obj):
        super().__init__()
        self._obj = obj
        self.ext = 'pdb'
        self.params = {}

    def get_structure_string(self):
        import prody

        def write(fname):
            if isinstance(self._obj, prody.Ensemble):
                st = self._obj[0]
            else:
                st = self._obj
            prody.writePDB(fname, st)

        return _get_structure_string(write)


@register_backend('prody')
class ProdyTrajectory(Trajectory, ProdyStructure):
    def __init__(self, obj):
        ProdyStructure.__init__(self, obj)

    @property
    def n_frames(self):
        return self._obj.numConfs()

    def get_coordinates(self, index):
        return self._obj.getConformation(index).getCoords()


@register_backend('simpletraj')
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
                "'SimpletrajTrajectory' requires the 'simpletraj' package")
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


@register_backend('mdtraj')
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
        return 10 * self.trajectory.xyz[index]

    @property
    def n_frames(self):
        return self.trajectory.n_frames

    def get_structure_string(self):
        return _get_structure_string(self.trajectory[0].save_pdb)


@register_backend('pytraj')
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

    def get_structure_string(self, index=0):
        return _get_structure_string(self.trajectory[index:index + 1].save)


@register_backend('parmed')
class ParmEdStructure(Structure):
    def __init__(self, structure):
        self._structure = structure
        self.only_save_1st_model = True

    def get_structure_string(self):
        # only write 1st model
        with mkstemp_wrapper(suffix='.pdb') as fname:
            if self.only_save_1st_model:
                self._structure.write_pdb(
                    fname, coordinates=self._structure.coordinates)
            else:
                self._structure.write_pdb(fname)
            with open(fname) as fh:
                return fh.read()


@register_backend('parmed')
class ParmEdTrajectory(Trajectory, ParmEdStructure):
    '''ParmEd adaptor.
    '''

    def __init__(self, trajectory):
        ParmEdStructure.__init__(self, trajectory)
        self.trajectory = self._structure = trajectory
        self.ext = "pdb"
        self.params = {}
        # only call get_coordinates once
        self._xyz = trajectory.get_coordinates()
        self.id = str(uuid.uuid4())

    def get_coordinates(self, index):
        return self._xyz[index]

    @property
    def n_frames(self):
        return len(self._xyz)


@register_backend('MDAnalysis')
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
        except ImportError:
            raise ImportError(
                "'MDAnalysisTrajectory' requires the 'MDAnalysis' package")
        u = self.atomgroup.universe
        u.trajectory[0]
        f = mda.lib.util.NamedStream(StringIO(), 'tmp.pdb')
        atoms = self.atomgroup.atoms
        # add PDB output to the named stream
        with mda.Writer(f, atoms.n_atoms, multiframe=False) as w:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                w.write(atoms)
        # extract from the stream
        pdb_string = f.read()
        return pdb_string


@register_backend('htmd')
class HTMDTrajectory(Trajectory):
    '''HTMD adaptor.

    Takes a Molecule object.

    Example
    -------
    >>> import nglview as nv
    >>> from htmd import Molecule
    >>> mol = Molecule(nv.datafiles.PDB)
    >>> t = nv.HTMDTrajectory(mol)
    >>> w = nv.NGLWidget(t)
    >>> w
    '''

    def __init__(self, mol):
        self.mol = mol
        self.ext = "pdb"
        self.params = {}
        self.id = str(uuid.uuid4())

    def get_coordinates(self, index):
        return np.squeeze(self.mol.coords[:, :, index])

    @property
    def n_frames(self):
        return self.mol.numFrames

    def get_structure_string(self):
        return _get_structure_string(self.mol.write)


@register_backend('ase')
class ASETrajectory(Trajectory, Structure):
    '''asetraj adaptor.

    Examples
    --------
    >>> import nglview as nv
    >>> from ase.io.trajectory import Trajectory
    >>> traj = Trajectory(nv.datafiles.ASE_Traj)
    >>> t = nv.ASETrajectory(traj)
    >>> w = nv.NGLWidget(t)
    >>> w.add_spacefill()
    >>> w
    '''

    def __init__(self, ase_traj):
        self.ext = 'pdb'
        self.params = {}
        self.trajectory = ase_traj
        self.id = str(uuid.uuid4())

    def get_coordinates(self, index):
        return self.trajectory[index].positions

    def get_structure_string(self):
        return _get_structure_string(self.trajectory[0].write)

    @property
    def n_frames(self):
        return len(self.trajectory)


@register_backend('schrodinger')
class SchrodingerStructure(Structure):
    '''Only read first structure
    '''

    def __init__(self, structure, ext="pdb"):
        super().__init__()
        self.path = ''
        self.ext = ext
        self.params = {}
        self._schrodinger_structure = structure

    def get_structure_string(self):
        # NOTE: For some reasons, _get_structure_string always return
        # empty string in this case.
        with tempfolder():
            pdb_fn = 'tmp.pdb'
            self._schrodinger_structure.write(pdb_fn)
            with open(pdb_fn) as fh:
                content = fh.read()
        return content

        return _get_structure_string(self._schrodinger_structure.write)


@register_backend('schrodinger')
class SchrodingerTrajectory(SchrodingerStructure, Trajectory):
    """Require `parmed` package.
    """

    def __init__(self, structure, traj):
        super().__init__(structure)
        self._traj = traj

    @property
    def n_frames(self):
        return len(self._traj) if self._traj else 1

    def get_coordinates(self, index):
        return self._traj[index].pos()

    def get_structure_string(self):
        """Require `parmed` package.
        """
        import parmed as pmd
        c = self._schrodinger_structure
        s = pmd.Structure()
        fsys = c.fsys_ct if hasattr(c, 'fsys_ct') else c
        for atom in fsys.atom:
            parm_atom = pmd.Atom(name=atom.pdbname.strip(),
                                 atomic_number=atom.atomic_number)
            s.add_atom(parm_atom,
                       atom.pdbres.strip(),
                       atom.resnum,
                       chain=atom.chain)
        s.coordinates = fsys.getXYZ()
        return ParmEdStructure(s).get_structure_string()

    @classmethod
    def from_files(cls, cms_fname, traj_fname):
        from schrodinger.application.desmond.packages import topo, traj
        _, model = topo.read_cms(cms_fname)
        traj = traj.read_traj(traj_fname)
        return cls(model, traj)
