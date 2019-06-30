import subprocess
import sys

import numpy as np
from mock import MagicMock

import nglview
from make_dummy_comm import *

# TODO : add more show_xxx
# (check test_widgets.py)


def _write(*args, **kargs):
    # fake write method
    subprocess.check_call(['cp', nglview.datafiles.PDB, 'tmp.pdb'])


class MockStructure:
    def as_pdb_string(self):
        with open(nglview.datafiles.PDB) as fh:
            return fh.read()


class MockRosettaPose:
    def dump_pdb(self, _):
        _write()


def test_show_schrodinger():
    # Show a schrodinger.structure.Structure
    s = MagicMock()
    s.write = _write
    v0 = nglview.show_schrodinger(s)

    # Show a trajectory with a Structure as topology
    class MockAtom:
        pdbname = 'C'
        atomic_number = 6
        pdbres = 'GLU'
        resnum = 1
        chain = 'A'

    s = MagicMock()
    s.fsys_ct.getXYZ.return_value = np.zeros((3, 3))
    s.fsys_ct.atom = [MockAtom()] * 3
    traj = MagicMock()
    traj.__len__.return_value = 3
    v1 = nglview.show_schrodinger(s, traj)
    assert v1.max_frame == 2
    v1.frame = 2  # trigger `get_coordinates`


def test_show_htmd():
    mol = MagicMock()
    mol.write = _write
    n_frames = 10
    mol.coordinates = np.zeros((n_frames, 1000, 3))
    mol.numFrames = n_frames
    view = nglview.show_htmd(mol)
    view


def test_show_rosetta():
    pose = MockRosettaPose()
    view = nglview.show_rosetta(pose)


def test_show_iotbx():
    mol = MockStructure()
    nglview.show_iotbx(mol)


def test_show_iodata():
    class MockIO:
        def to_file(self, fname):
            with open(fname, 'w') as fh:
                fh.write("""3
water
H    0.7838370000   -0.4922360000   -0.0000000000
O   -0.0000000000    0.0620200000   -0.0000000000
H   -0.7838370000   -0.4922360000   -0.0000000000""")

    v = nglview.show_iodata(MockIO())
    v


def test_show_qcelemental_show_psi4():
    class MockMol:
        def to_string(self, format):
            return '1\nHe\nHe                    0.000000000000     0.000000000000     0.000000000000\n'

    v = nglview.show_qcelemental(MockMol())
    v
    nglview.show_psi4(MockMol())


def test_show_openbabel():
    import types
    openbabel = types.ModuleType('openbabel')
    sys.modules['openbabel'] = openbabel
    b = MagicMock()
    mol = MagicMock()
    openbabel.OBConversion = MagicMock(return_value=b)
    nglview.show_openbabel(mol)
    b.SetOutFormat.assert_called_with('pdb')
    assert b.WriteFile.called
    assert len(b.WriteFile.call_args_list[0]) == 2


def test_show_prody():
    import types
    prody = types.ModuleType('prody')
    sys.modules['prody'] = prody

    class MockEnsemble:
        def __getitem__(self, index):
            return 0

        def numConfs(self):
            return 1

        def getConformation(self, index):
            class Struct:
                def getCoords(self):
                    return

    prody.Ensemble = MockEnsemble
    prody.writePDB = MagicMock()
    ens = MockEnsemble()
    nglview.show_prody(ens)
    assert prody.writePDB.called
    prody.writePDB.reset_mock()
    assert not prody.writePDB.called

    st = MagicMock()
    nglview.show_prody(st)
    assert prody.writePDB.called
