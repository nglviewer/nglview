import subprocess
from mock import MagicMock
import nglview
import numpy as np

# local
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
    s = MagicMock()
    s.write = _write
    nglview.show_schrodinger(s)


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
