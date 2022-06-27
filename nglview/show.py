from io import StringIO

from . import datafiles
from .adaptor import (ASEStructure, ASETrajectory, BiopythonStructure,
                      FileStructure, HTMDTrajectory, IODataStructure,
                      IOTBXStructure, MDAnalysisTrajectory, MDTrajTrajectory,
                      OpenbabelStructure, ParmEdTrajectory, PdbIdStructure,
                      ProdyStructure, ProdyTrajectory, PyTrajTrajectory,
                      QCElementalStructure, RosettaStructure,
                      SchrodingerStructure, SchrodingerTrajectory,
                      RdkitStructure,
                      TextStructure)
from .widget import NGLWidget

__all__ = [
    'demo',
    'show_pdbid',
    'show_url',
    'show_text',
    'show_ase',
    'show_pymatgen',
    'show_iotbx',
    'show_iodata',
    'show_qcelemental',
    'show_openbabel',
    'show_psi4',
    'show_rosetta',
    'show_asetraj',
    'show_simpletraj',
    'show_prody',
    'show_mdtraj',
    'show_pytraj',
    'show_mdanalysis',
    'show_parmed',
    'show_rdkit',
    'show_structure_file',
    'show_file',
    'show_htmd',
    'show_schrodinger',
    'show_biopython',
]


def show_pdbid(pdbid, **kwargs):
    '''Show PDB entry.

    Examples
    --------
    >>> import nglview as nv
    >>> w = nv.show_pdbid("3pqr")
    >>> w # doctest: +SKIP
    '''
    structure = PdbIdStructure(pdbid)
    return NGLWidget(structure, **kwargs)


def show_url(url, **kwargs):
    kwargs2 = {k: v for k, v in kwargs.items()}
    view = NGLWidget()
    view.add_component(url, **kwargs2)
    return view


def show_text(text, **kwargs):
    """for development
    """
    structure = TextStructure(text)
    return NGLWidget(structure, **kwargs)


def show_ase(ase_atoms, **kwargs):
    """

    Examples
    --------
    >>> import nglview as nv
    >>> from ase import Atom, Atoms
    >>> dimer = Atoms([Atom('X', (0, 0, 0)),
    ...                Atom('X', (0, 0, 1))])
    >>> dimer.set_positions([(1, 2, 3), (4, 5, 6.2)])
    >>> w = nv.show_ase(dimer)
    >>> w # doctest: +SKIP
    """
    structure = ASEStructure(ase_atoms)
    return NGLWidget(structure, **kwargs)


def show_iodata(obj, **kwargs):
    """Show iodata's IOData (require `ase` package).

    Examples
    --------
    >>> import nglview as nv
    >>> from iodata import IOData # doctest: +SKIP
    ... obj = IOData.from_file('what.xyz')
    ... view = nv.show_iodata(obj)
    ... view
    """
    return NGLWidget(IODataStructure(obj), **kwargs)


def show_qcelemental(obj, **kwargs):
    """Show QCElemental's Molecule.

    Examples
    --------
    >>> import nglview as nv
    >>> import qcelemental as qcel # doctest: +SKIP
    ... mol = qcel.models.Molecule.from_data("He 0 0 0")
    ... view = nv.show_qcelemental(mol)
    ... view
    """
    return NGLWidget(QCElementalStructure(obj), **kwargs)


def show_psi4(obj, **kwargs):
    """Show Psi4's Molecule.

    Examples
    --------
    >>> import nglview as nv
    >>> import psi4 # doctest: +SKIP
    ... mol = psi4.geometry('xyz content here')
    ... view = nv.show_psi4(mol)
    ... view
    """
    return NGLWidget(QCElementalStructure(obj), **kwargs)


def show_openbabel(obj, **kwargs):
    """Show openbabel's Omol

    >>> import nglview
    >>> mol = openbabel.OBMol() # doctest: +SKIP
    ... obConversion = openbabel.OBConversion() 
    ... obConversion.SetInFormat('xyz')
    ... obConversion.ReadFile(mol, 'what.xyz')
    ... nglview.show_openbabel(mol)
    """
    return NGLWidget(OpenbabelStructure(obj), **kwargs)


def show_pymatgen(struct, **kwargs):
    """Require `ase` package.

    Examples
    --------
    >>> import nglview as nv # doctest: +SKIP
    ... import pymatgen as mg
    ... lattice = mg.Lattice.cubic(4.2)
    ... structure = mg.Structure(lattice, ["Cs", "Cl"],
                          [[0, 0, 0], [0.5, 0.5, 0.5]])
    ... view = nv.show_pymatgen(structure)
    ... view
    """
    from pymatgen.io.ase import AseAtomsAdaptor
    return show_ase(AseAtomsAdaptor().get_atoms(struct))


def show_iotbx(mol, **kwargs):
    """

    Examples
    --------
    >>> import iotbx.pdb
    ... x = iotbx.pdb.hierarchy.input('file.pdb')
    ... mol = x.construct_hierarchy()
    ... view = nglview.show_iotbx(hi)
    ... view # doctest: +SKIP
    """
    structure = IOTBXStructure(mol)
    return NGLWidget(structure, **kwargs)


def show_rosetta(pose, **kwargs):
    """

    Examples
    --------
    >>> from pyrosetta import pose_from_sequence, init
    ... init()
    ... pose = pose_from_sequence('AAAAAA')
    ... view = nglview.show_rosetta(pose)
    ... view # doctest: +SKIP
    """
    structure = RosettaStructure(pose)
    return NGLWidget(structure, **kwargs)


def show_asetraj(ase_traj, **kwargs):
    '''Show ase trajectory and structure file.

    Examples
    --------
    >>> import nglview as nv
    >>> from ase.io.trajectory import Trajectory
    >>> traj = Trajectory(nv.datafiles.ASE_Traj)
    >>> view = nv.show_asetraj(traj)
    >>> view.add_spacefill()
    >>> view # doctest: +SKIP
    '''
    trajectory = ASETrajectory(ase_traj)
    return NGLWidget(trajectory, **kwargs)


def show_structure_file(path, **kwargs):
    '''Show structure file. Allowed are text-based structure
    file formats that are by supported by NGL, including pdb,
    gro, mol2, sdf.

    Examples
    --------
    >>> import nglview as nv
    >>> w = nv.show_structure_file(nv.datafiles.GRO)
    >>> w # doctest: +SKIP
    '''
    structure = FileStructure(path)
    return NGLWidget(structure, **kwargs)


def _show_schrodinger_file(path, **kwargs):
    from schrodinger.structure import StructureReader

    view = NGLWidget()
    cts = list(StructureReader(path))
    for ct in cts:
        view.add_component(SchrodingerStructure(ct), **kwargs)
    return view


def show_file(path, **kwargs):
    '''Show any supported file format (e.g: .gro, .dx, ...)

    Examples
    --------
    >>> import nglview as nv
    ... w = nv.show_file('my.dx')
    ... w # doctest: +SKIP
    '''
    if isinstance(path, str) and (path.endswith(
        ('.mae', '.mae.gz', '.meagz', '.cms', '.cms.gz', '.cms.gz')) or
                                  (kwargs.get('format', None) == 'maestro')):
        view = _show_schrodinger_file(path, **kwargs)
    else:
        view = NGLWidget()
        view.add_component(path, **kwargs)
    return view


def show_simpletraj(traj, **kwargs):
    '''Show simpletraj trajectory and structure file.

    Examples
    --------
    >>> import nglview as nv
    >>> traj = nv.SimpletrajTrajectory(nv.datafiles.XTC, nv.datafiles.GRO)
    >>> view = nv.show_simpletraj(traj)
    >>> view # doctest: +SKIP
    '''
    return NGLWidget(traj, **kwargs)


def show_prody(obj, **kwargs):
    """

    Examples
    --------
    >>> import nglview as nv # doctest: +SKIP
    ... import prody
    ... structure = prody.parsePDB('what.pdb')
    ... ensemble = prody.parseDCD('what.dcd')
    ... ensemble.setAtoms(structure)
    ... nv.show_prody(ensemble)
    ... # nv.show_prody(structure)
    """
    import prody
    if isinstance(obj, prody.Ensemble):
        view_obj = ProdyTrajectory(obj)
    else:
        view_obj = ProdyStructure(obj)
    return NGLWidget(view_obj, **kwargs)


def show_mdtraj(mdtraj_trajectory, **kwargs):
    '''Show mdtraj trajectory.

    Examples
    --------
    >>> import nglview as nv # doctest: +SKIP
    ... import mdtraj as md
    ... t = md.load(nv.datafiles.XTC, top=nv.datafiles.GRO)
    ... w = nv.show_mdtraj(t)
    ... w
    '''
    structure_trajectory = MDTrajTrajectory(mdtraj_trajectory)
    return NGLWidget(structure_trajectory, **kwargs)


def show_pytraj(pytraj_trajectory, **kwargs):
    '''Show pytraj trajectory.

    Examples
    --------
    >>> import nglview as nv # doctest: +SKIP
    ... import pytraj as pt
    ... t = pt.load(nv.datafiles.TRR, nv.datafiles.PDB)
    ... w = nv.show_pytraj(t)
    ... w
    '''
    trajlist = pytraj_trajectory if isinstance(pytraj_trajectory,
                                               (list, tuple)) else [
                                                   pytraj_trajectory,
                                               ]

    trajlist = [PyTrajTrajectory(traj) for traj in trajlist]
    return NGLWidget(trajlist, **kwargs)


def show_parmed(parmed_structure, **kwargs):
    '''Show pytraj trajectory.

    Examples
    --------
    >>> import nglview as nv
    >>> import parmed as pmd
    >>> t = pmd.load_file(nv.datafiles.PDB)
    >>> w = nv.show_parmed(t)
    >>> w # doctest: +SKIP
    '''
    structure_trajectory = ParmEdTrajectory(parmed_structure)
    return NGLWidget(structure_trajectory, **kwargs)


def show_rdkit(rdkit_mol, **kwargs):
    '''Show rdkit's Mol.

    Parameters
    ----------
    rdkit_mol : rdkit.Chem.rdchem.Mol
    kwargs : additional keyword argument
    If kwargs contains the "conf_id" key, this will be passed to the
    RDKit Chem.MolToXXXBlock function as the confId parameter.
    If the "conf_id" key is not provided, -1 will be used as confId.
    If kwargs contains the "fmt" key, this will be used to decide
    whether rdkit_mol should be visualized as a PDB block (fmt == "pdb")
    or as a SDF block (fmt == "sdf").
    If the "fmt" key is not provided, a simple heuristic is used:
    if the first atom contains PDB residue information, rdkit_mol
    is visualized as a PDB block, otherwise as a SDF block.

    Examples
    --------
    >>> import nglview as nv
    >>> from rdkit import Chem # doctest: +SKIP
    ... from rdkit.Chem import AllChem
    ... m = Chem.AddHs(Chem.MolFromSmiles('COc1ccc2[C@H](O)[C@@H](COc2c1)N3CCC(O)(CC3)c4ccc(F)cc4'))
    ... _ = AllChem.EmbedMultipleConfs(m, useExpTorsionAnglePrefs=True, useBasicKnowledge=True)
    ... view = nv.show_rdkit(m)
    ... view # doctest: +SKIP

    >>> # add component m2
    >>> # create file-like object
    >>> from nglview.show import StringIO
    >>> m2 = Chem.AddHs(Chem.MolFromSmiles('N[C@H](C)C(=O)O')) # doctest: +SKIP
    ... fh = StringIO(Chem.MolToPDBBlock(m2))
    ... view.add_component(fh, ext='pdb')

    >>> # load as trajectory, need to have ParmEd
    >>> view = nv.show_rdkit(m, parmed=True) # doctest: +SKIP
    '''
    ext = kwargs.pop("fmt", "pdb")
    conf_id = kwargs.pop("conf_id", -1)
    struc = RdkitStructure(rdkit_mol, ext=ext, conf_id=conf_id)
    return  NGLWidget(struc, **kwargs)


def show_mdanalysis(atomgroup, **kwargs):
    '''Show NGL widget with MDAnalysis AtomGroup.

    Can take either a Universe or AtomGroup as its data input.

    Examples
    --------
    >>> import nglview as nv # doctest: +SKIP
    ... import MDAnalysis as mda
    ... u = mda.Universe(nv.datafiles.GRO, nv.datafiles.XTC)
    ... prot = u.select_atoms('protein')
    ... w = nv.show_mdanalysis(prot)
    ... w
    '''
    structure_trajectory = MDAnalysisTrajectory(atomgroup)
    return NGLWidget(structure_trajectory, **kwargs)


def show_htmd(mol, **kwargs):
    '''Show NGL widget with HTMD Molecule.

    Takes a Molecule object as its data input.

    Examples
    --------
    >>> import nglview as nv # doctest: +SKIP
    ... from htmd import Molecule
    ... mol = Molecule(nv.datafiles.PDB)
    ... mol.filter('protein')
    ... w = nv.show_htmd(mol)
    ... w
    '''
    structure_trajectory = HTMDTrajectory(mol)
    return NGLWidget(structure_trajectory, **kwargs)


def show_schrodinger(mol, traj=None, **kwargs):
    '''Show NGL widget with Schrodinger's Structure

    Notes
    -----
    EXPERIMENTAL

    Examples
    --------
    >>> import nglview as nv # doctest: +SKIP
    ... from schrodinger.structure import StructureReader
    ... for s in StructureReader(fn):
    ...    break
    ... w = nv.show_schrodinger(s)
    ... w
    '''
    if traj is None:
        structure_trajectory = SchrodingerStructure(mol)
    else:
        structure_trajectory = SchrodingerTrajectory(mol, traj)
    return NGLWidget(structure_trajectory, **kwargs)


def show_biopython(entity, **kwargs):
    '''Show NGL widget with Biopython structural entity.

    Takes a Structure, Model, Chain, Residue or Atom
    from Bio.PDB as its data input.

    Examples
    --------
    >>> import nglview as nv # doctest: +SKIP
    ... from Bio.PDB import PDBParser
    ... parser = PDBParser()
    ... structure = parser.get_structure("protein", "protein.pdb")
    ... w = nv.show_biopython(structure[0]["A"])
    ... w
    '''
    entity = BiopythonStructure(entity)
    return NGLWidget(entity, **kwargs)


def demo(*args, **kwargs):
    '''

    Examples
    --------
    >>> import nglview as nv
    >>> view = nv.demo()
    >>> view # doctest: +SKIP
    '''
    from nglview import show_structure_file
    return show_structure_file(datafiles.PDB, *args, **kwargs)
