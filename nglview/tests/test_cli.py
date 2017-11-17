import sys
import os
import unittest
import subprocess
from nglview import datafiles
import pytest
from nglview.scripts.nglview import main, get_remote_port
from mock import patch, mock_open

this_path = os.path.dirname(os.path.abspath(__file__))
PY2 = sys.version_info[0] == 2

@patch('subprocess.check_call')
def test_cli(mock_call):
    # no argument 
    command = []
    # subprocess.check_call(command.split())
    main(cmd=command)

    # only nglview and --clean
    command = ['--clean']
    main(cmd=command)

    # demo
    command = ['demo']
    main(cmd=command)

    # single pdb
    command = [datafiles.PDB]
    main(cmd=command)

    # single pdb, disable-autorun
    command = [datafiles.PDB, '--auto']
    main(cmd=command)

    # single pdb, specify browser
    command = 'nglview {} --browser=google-chrome'.format(datafiles.PDB)
    subprocess.check_call(command.split())

    # single pdb, does not exists
    command = ['test_hellooooo.pdb']
    with pytest.raises(AssertionError):
        main(cmd=command)

    # density data
    densityfile = os.path.join(this_path, 'data/volmap.dx')
    command = [densityfile]
    main(cmd=command)

    # python script
    command = ['my.py']
    if PY2:
        with patch('__builtin__.open'), patch('json.dumps'):
            main(cmd=command)
    else:
        with patch('builtins.open'), patch('json.dumps'):
            main(cmd=command)

    # pytraj
    command = 'nglview {} -c {}'.format(datafiles.PDB, datafiles.XTC)
    subprocess.check_call(command.split())

    # remote
    command = 'nglview {} -c {} --remote'.format(
        datafiles.PDB, datafiles.XTC)
    subprocess.check_call(command.split())

    # python script
    pyfile = os.path.join(this_path, 'test_widget.py')
    command = 'nglview {pyfile}'.format(pyfile=pyfile)
    subprocess.check_call(command.split())

    # notebook (.ipynb)
    nbfile = os.path.join(this_path, 'notebooks/api/test_detach.ipynb')
    command = [nbfile]
    main(cmd=command)

    # install
    command = ['install', '--symlink']
    with pytest.raises(SystemExit):
        main(cmd=command)
    mock_call.assert_called_with(
            ['jupyter', 'nbextension', 'install', '--py', '--sys-prefix', 'nglview',
             '--overwrite', '--symlink'])

    # enable
    command = ['enable']
    with pytest.raises(SystemExit):
        main(cmd=command)
    mock_call.assert_called_with(
            ['jupyter', 'nbextension', 'enable', '--py', '--sys-prefix', 'nglview'])

    # remote
    with patch('nglview.scripts.nglview.get_remote_port') as mock_remote:
        command = ['--remote']
        main(cmd=command)
        mock_remote.assert_called_with(None, 'tmpnb_ngl.ipynb')



@patch('nglview.scripts.app.NGLViewApp.get_port')
@patch('os.getlogin')
@patch('socket.gethostname')
def test_get_port(mock_gethostname, mock_getlogin, mock_get_port):
    mock_gethostname.return_value = 'myhost'
    mock_getlogin.return_value = 'mylogin'
    mock_get_port.return_value = 9999
    get_remote_port(None, 'my.ipynb')
