#!/usr/bin/env python

from __future__ import absolute_import
import os, sys, argparse, json
import subprocess
from subprocess import CalledProcessError
from .cmd_example import CMD_EXAMPLE 

bin_path = os.path.join(sys.prefix, 'bin')

simple_source = """
import nglview as nv
""".strip()

demo_source = """
import nglview as nv

view = nv.demo(gui=True)
view
""".strip()

density_source = """
import nglview as nv

view = nv.NGLWidget(gui=True)
view.add_component('filename')
view
""".strip()

def _is_density_data(filename):
    from nglview.utils import FileManager

    try:
        fm = FileManager(filename)

        return fm.ext.lower() in ['dx', 'ccp4', 'mrc', 'map', 'dxbin', 'cube']
    except ValueError:
        return False

notebook_dict = {
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 'null',
   "metadata": {
    "collapsed": True
   },
   "outputs": [],
   "source": [
    "import nglview as nv\n",
    "import pytraj as pt\n",
    "\n",
    "traj = pt.iterload('test.nc', top='prmtop')\n",
    "view = nv.show_pytraj(traj)\n",
    "view"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython",
   "version": "3.5.1"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 0
}

def install_nbextension(jupyter, user=True):
    path = os.path.dirname(__file__)
    nglview_main = os.path.join(path, 'nglview_main.js')

    local = '--user' if user else ''
    cm_install = '{jupyter} nbextension install {nglview_main} {local}'.format(jupyter=jupyter,
            nglview_main=nglview_main,
            local=local)
    cm_activate = '{jupyter} nbextension enable nglview_main'.format(jupyter=jupyter) 

    with open(os.devnull, 'wb') as devnull:
        subprocess.check_call(cm_install.split(),
            stdout=devnull,
            stderr=subprocess.STDOUT)
        subprocess.check_call(cm_activate.split(),
            stdout=devnull,
            stderr=subprocess.STDOUT)

def disable_extension(jupyter):
    print("disable nglview_main extension")
    cm = '{jupyter} nbextension disable nglview_main'.format(jupyter=jupyter)
    subprocess.check_call(cm.split())


remote_msg = """
Try to use port = {port}

\033[32m In your local machine, run: \033[0m

    {client_cm}

\033[32m NOTE: you might want to replace {hostname} by full hostname with domain name \033[0m

\033[32m Then open your web browser, copy and paste: \033[0m
    http://localhost:{port}/notebooks/{notebook_name}
"""

def get_remote_port(port=None, notebook_path=''):
    import os, socket
    from nglview.scripts.app import NGLViewApp
    port = NGLViewApp().get_port(port=port)

    username = os.getlogin()
    hostname = socket.gethostname()
    client_cm = "ssh -NL localhost:{port}:localhost:{port} {username}@{hostname}".format(username=username,
            hostname=hostname,
            port=port)
    base_notebook_name = os.path.basename(notebook_path)
    print(remote_msg.format(client_cm=client_cm, port=port,
                            hostname=hostname,
                            notebook_name=base_notebook_name))
    return port

def main(notebook_dict=notebook_dict, cmd_arg=sys.argv[1:]):
    PY3 = sys.version_info[0] == 3
    pyv_full_string = ','.join(str(i) for i in sys.version_info)
    pyv_short_string = str(sys.version_info[0])
    default_jexe = ' '.join((sys.executable, '-m jupyter_core'))

    parser = argparse.ArgumentParser(description='NGLView: An IPython/Jupyter widget to '
                                     'interactively view molecular structures and trajectories.',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=CMD_EXAMPLE)
    parser.add_argument('command', nargs='?',
            help='command could be a topology filename (.pdb, .mol2, .parm7, ...) or \n'
                          'could be a python script (.py), a notebook (.ipynb). '
                          'If not given, a notebook will be created with only nglview imported') 
    parser.add_argument('traj', nargs='?', help='coordinate filename, optional')
    parser.add_argument('-c', '--crd', help='coordinate filename')
    parser.add_argument('--browser', help='web browser')
    parser.add_argument('-j', '--jexe', default=default_jexe, help='jupyter path')
    parser.add_argument('--notebook-name', default='tmpnb_ngl.ipynb', help='notebook name')
    parser.add_argument('--port', type=int, help='port number')
    parser.add_argument('--remote', action='store_true', help='create remote notebook')
    parser.add_argument('--clean', action='store_true', help='delete temp file after closing notebook')
    parser.add_argument('--disable-autorun', action='store_true', help='do not run 1st cell right after openning notebook')
    parser.add_argument('--test', action='store_true', help='test')
    args = parser.parse_args(cmd_arg)

    command = args.command

    crd = args.traj if args.traj is not None else args.crd
    if crd is None:
        crd = command

    browser = '--browser ' + args.browser if args.browser else ''

    create_new_nb = False

    if command is not None and command.endswith('.ipynb'):
        notebook_name = command
    else:
        notebook_name = args.notebook_name
        if not PY3:
            kernelspec = notebook_dict['metadata']['kernelspec']
            kernelspec['display_name'] = 'Python 2'
            kernelspec['name'] = 'python2'
            notebook_dict['metadata']['kernelspec'] = kernelspec

            notebook_dict['metadata']['language_info']['codemirror_mode']['version'] = pyv_short_string
            notebook_dict['metadata']['version'] = pyv_full_string
        if command is None:
            # create a notebook and import nglview
            notebook_dict['cells'][0]['source'] = simple_source
            nb_json = json.dumps(notebook_dict)
            nb_json = nb_json.replace('"null"', 'null')
        elif command.endswith('.py'):
            # a Python script
            pycontent = open(command).read().strip()
            notebook_dict['cells'][0]['source'] = pycontent
            nb_json = json.dumps(notebook_dict)
        elif command == 'demo':
            # running demo
            notebook_dict['cells'][0]['source'] = demo_source
            nb_json = json.dumps(notebook_dict)
            nb_json = nb_json.replace('"null"', 'null').replace('test.nc', crd).replace('prmtop', command)
        elif _is_density_data(command):
            # check if density data
            notebook_dict['cells'][0]['source'] = density_source.replace('filename',
                    command)
            nb_json = json.dumps(notebook_dict)
        else:
            nb_json = json.dumps(notebook_dict)
            nb_json = nb_json.replace('"null"', 'null').replace('test.nc',
                    crd).replace('prmtop', command)
            assert os.path.exists(command), '{} does not exists'.format(command)

        nb_json = nb_json.replace('"null"', 'null')

        with open(notebook_name, 'w') as fh:
            fh.write(nb_json)
            create_new_nb = True
    
    
    dirname = os.path.dirname(os.path.abspath(notebook_name))
    if not args.remote:
        cm = '{jupyter} notebook {notebook_name} {browser}'.format(jupyter=args.jexe,
                                                                   notebook_name=notebook_name,
                                                                   browser=browser)
    else:
        port = get_remote_port(args.port, notebook_name)
        cm = '{jupyter} notebook --no-browser --port {port} ' \
              '--notebook-dir {dirname}'.format(jupyter=args.jexe,
                                                port=port,
                                                dirname=dirname)
        print('NOTE: make sure to open {0} in your local machine\n'.format(notebook_name))

    if args.disable_autorun or command is None:
        try:
            disable_extension(jupyter=args.jexe)
        except CalledProcessError:
            pass
    else:
        if not args.test:
            install_nbextension(jupyter=args.jexe)

    try:
        if args.test:
            # for testing cli
            pass
        else:
            subprocess.check_call(cm.split())
    except KeyboardInterrupt:
        if args.clean and create_new_nb:
            print("deleting {}".format(notebook_name))
            os.remove(notebook_name)
        if not args.disable_autorun:
            disable_extension(jupyter=args.jexe)

if __name__ == '__main__':
    main(cmd_arg=sys.argv[1:])
