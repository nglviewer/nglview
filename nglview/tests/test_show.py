import sys
import subprocess

import unittest
import numpy as np
from mock import MagicMock, patch

import nglview


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

        def orient_molecule(self):
            return self

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



try:
    import rdkit
    from rdkit import Chem
    from rdkit.Chem import AllChem
    has_rdkit = True
except ImportError:
    rdkit = AllChem = None
    has_rdkit = False


@unittest.skipUnless(has_rdkit, 'must have rdkit')
def test_show_rdkit():
    import rdkit.Chem as Chem
    rdkit_mol = Chem.AddHs(
        Chem.MolFromSmiles(
            'COc1ccc2[C@H](O)[C@@H](COc2c1)N3CCC(O)(CC3)c4ccc(F)cc4'))
    AllChem.EmbedMultipleConfs(rdkit_mol,
                               useExpTorsionAnglePrefs=True,
                               useBasicKnowledge=True)
    # FIXME: create test_adaptor.py?
    # Test RdkitStructure
    structure = nglview.RdkitStructure(rdkit_mol)
    assert "HETATM" in structure.get_structure_string()
    structure = nglview.RdkitStructure(rdkit_mol, ext="sdf")
    assert "RDKit          3D" in structure.get_structure_string()
    structure2 = nglview.RdkitStructure(rdkit_mol, ext="sdf", conf_id=0)
    assert "RDKit          3D" in structure2.get_structure_string()
    assert structure.get_structure_string() == structure2.get_structure_string()
    structure3 = nglview.RdkitStructure(rdkit_mol, ext="sdf", conf_id=1)
    assert "RDKit          3D" in structure3.get_structure_string()
    assert structure.get_structure_string() != structure3.get_structure_string()

    # Test show_rdkit
    with patch.object(Chem, 'MolToPDBBlock') as mock_MolToPDBBlock,\
            patch.object(Chem, 'MolToMolBlock') as mock_MolToMolBlock:
        nglview.show_rdkit(rdkit_mol, fmt='sdf', conf_id=1)
        assert mock_MolToMolBlock.called
        assert not mock_MolToPDBBlock.called
        assert mock_MolToMolBlock.call_args_list[-1][-1] == {'confId': 1}
    with patch.object(Chem, 'MolToPDBBlock') as mock_MolToPDBBlock,\
            patch.object(Chem, 'MolToMolBlock') as mock_MolToMolBlock:
        nglview.show_rdkit(rdkit_mol)
        assert not mock_MolToMolBlock.called
        assert mock_MolToPDBBlock.called
        assert mock_MolToPDBBlock.call_args_list[-1][-1] == {'confId': -1}
