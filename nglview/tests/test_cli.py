import os
import subprocess
from nglview import datafiles
import pytest

this_path = os.path.dirname(os.path.abspath(__file__))

def test_cli():
    # demo
    command = 'nglview demo --test'
    subprocess.check_call(command.split())

    # raises
    command = 'nglview demo --test --dummy-args'
    def func():
        print('raise error for "{}"'.format(command))
        subprocess.check_call(command.split())
    pytest.raises(subprocess.CalledProcessError, func)

    # single pdb
    command = 'nglview {} --test'.format(datafiles.PDB)
    subprocess.check_call(command.split())

    # single pdb, disable-autorun
    command = 'nglview {} --disable-autorun --test'.format(datafiles.PDB)
    subprocess.check_call(command.split())

    # single pdb, specify browser
    command = 'nglview {} --browser=google-chrome --test'.format(datafiles.PDB)
    subprocess.check_call(command.split())

    # single pdb, does not exists
    command = 'nglview test_hellooooo.pdb --test'
    pytest.raises(subprocess.CalledProcessError, func)

    # density data
    densityfile = os.path.join(this_path, 'data/volmap.dx')
    command = 'nglview {} --test'.format(densityfile)
    subprocess.check_call(command.split())

    # pytraj
    command = 'nglview {} -c {} --test'.format(datafiles.PDB, datafiles.XTC)
    subprocess.check_call(command.split())

    # remote
    command = 'nglview {} -c {} --remote --test'.format(datafiles.PDB, datafiles.XTC)
    subprocess.check_call(command.split())

    # python script
    pyfile = os.path.join(this_path, 'test_widget.py')
    command = 'nglview {pyfile} --test'.format(pyfile=pyfile)
    subprocess.check_call(command.split())

    # notebook (.ipynb)
    nbfile = os.path.join(this_path, 'notebooks/api/test_detach.ipynb')
    command = 'nglview {nbfile} --test'.format(nbfile=nbfile)
    subprocess.check_call(command.split())


if __name__ == '__main__':
    test_cli()
