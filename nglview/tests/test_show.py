import subprocess
from mock import MagicMock, patch
import nglview

# local
from make_dummy_comm import *
# TODO : add more show_xxx

def _write(*args, **kargs):
    # fake write method
    subprocess.check_call([
        'cp', nglview.datafiles.PDB,
        'tmp.pdb'
    ])

def test_show_schrodinger_structure():
    s = MagicMock()
    s.write = _write
    nglview.show_schrodinger_structure(s)
