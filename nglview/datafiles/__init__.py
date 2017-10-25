import os
import os.path

MODULE_DIR = os.path.split(os.path.abspath(__file__))[0]
PDB = os.path.join(MODULE_DIR, "md_1u19.pdb")
GRO = os.path.join(MODULE_DIR, "md_1u19.gro")
XTC = os.path.join(MODULE_DIR, "md_1u19.xtc")
TRR = os.path.join(MODULE_DIR, "md_1u19.trr")
ASE_Traj = os.path.join(MODULE_DIR, "md_1u19.traj")
ALA3 = os.path.join(MODULE_DIR, "ala3.pdb")
