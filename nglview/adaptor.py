from __future__ import print_function, absolute_import

import os
import os.path
import uuid
import numpy as np

try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO

try:
    from urllib.request import urlopen
except ImportError:
    from urllib2 import urlopen

from .base_adaptor import Structure, Trajectory
from .utils.py_utils import FileManager, tempfolder
from . import config

__all__ = [
    'FileStructure', 'TextStructure', 'RdkitStructure', 'PdbIdStructure',
    'ASEStructure', 'BiopythonStructure', 'IOTBXStructure',
    'RosettaStructure',
    'SimpletrajTrajectory',
    'MDTrajTrajectory', 'PyTrajTrajectory', 'ParmEdTrajectory',
    'MDAnalysisTrajectory', 'HTMDTrajectory', 'ASETrajectory',
    'register_backend',
]


class register_backend(object):
    def __init__(self, package_name):
        # package_name must match exactly to your Python package
        self.package_name = package_name

    def __call__(self, cls):
        config.BACKENDS[self.package_name] = cls
        return cls


class FileStructure(Structure):
    def __init__(self, path):
        super(FileStructure, self).__init__()
        self.fm = FileManager(path)
        self.ext = self.fm.ext
        self.params = {}
        self.path = path

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
        return urlopen(url).read().decode('utf-8')


class ASEStructure(Structure):
    def __init__(self, ase_atoms, ext='pdb', params={}):
        super(ASEStructure, self).__init__()
        self.path = ''
        self.ext = ext
        self.params = params
        self._ase_atoms = ase_atoms

    def get_structure_string(self):
        with tempfolder():
            self._ase_atoms.write('tmp.pdb')
            return open('tmp.pdb').read()


class BiopythonStructure(Structure):
    def __init__(self, entity, ext='pdb', params={}):
        super(BiopythonStructure, self).__init__()
        self.path = ''
        self.ext = ext
        self.params = params
        self._entity = entity

    def get_structure_string(self):
        from Bio.PDB import PDBIO
        try:
            from StringIO import StringIO
        except ImportError:
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
        super(IOTBXStructure, self).__init__()
        self.path = ''
        self.ext = ext
        self.params = params
        self._mol = obj

    def get_structure_string(self):
        return self._mol.as_pdb_string()


class RosettaStructure(Structure):
    def __init__(self, pose, ext='pdb', params={}):
        # type: (pyrosetta.rosetta.core.pose.Pose, str, Dict) -> None
        super(RosettaStructure, self).__init__()
        self.path = ''
        self.ext = ext
        self.params = params
        self._mol = pose

    def get_structure_string(self):
        with tempfolder():
            self._mol.dump_pdb('tmp.pdb')
            return open('tmp.pdb').read()


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
        with tempfolder():
            fname = 'tmp.pdb'
            self.trajectory[0].save_pdb(fname)
            with open(fname) as fh:
                pdb_string = fh.read()
        return pdb_string


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

    def get_structure_string(self):
        fname = 'tmp.pdb'
        with tempfolder():
            self.trajectory[:1].save(fname, format="pdb", overwrite=True)
            with open(fname) as fh:
                pdb_string = fh.read()
        return pdb_string


@register_backend('parmed')
class ParmEdStructure(Structure):
    def __init__(self, structure):
        self._structure = structure
        self.only_save_1st_model = True

    def get_structure_string(self):
        # only write 1st model
        fname = 'tmp.pdb'
        with tempfolder():
            if self.only_save_1st_model:
                self._structure.save(
                    fname,
                    overwrite=True,
                    coordinates=self._structure.coordinates)
            else:
                self._structure.save(fname, overwrite=True)
            with open(fname) as fh:
                pdb_string = fh.read()
        return pdb_string


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
        with mda.Writer(f, atoms.n_atoms, multiframe=False) as W:
            W.write(atoms)
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
        fname = 'tmp.pdb'
        with tempfolder():
            self.mol.write(fname)
            with open(fname) as fh:
                pdb_string = fh.read()
        return pdb_string


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
        with tempfolder():
            self.trajectory[0].write('tmp.pdb')
            return open('tmp.pdb').read()

    @property
    def n_frames(self):
        return len(self.trajectory)


@register_backend('schrodinger')
class SchrodingerStructure(Structure):
    '''Only read first structure
    '''
    def __init__(self, structure, ext="pdb"):
        super(SchrodingerStructure, self).__init__()
        self.path = ''
        self.ext = ext
        self.params = {}
        self._schrodinger_structure = structure

    def get_structure_string(self):
        with tempfolder():
            pdb_fn = 'tmp.pdb'
            self._schrodinger_structure.write(pdb_fn)
            with open(pdb_fn) as fh:
                content = fh.read()
        return content


@register_backend('schrodinger')
class SchrodingerTrajectory(SchrodingerStructure, Trajectory):
    """Require `parmed` package.
    """
    def __init__(self, structure, traj):
        super(SchrodingerTrajectory, self).__init__(structure)
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
            parm_atom = pmd.Atom(name=atom.pdbname.strip(), atomic_number=atom.atomic_number)
            s.add_atom(parm_atom, atom.pdbres.strip(), atom.resnum, chain=atom.chain)
        s.coordinates = fsys.getXYZ()
        return ParmEdStructure(s).get_structure_string()
