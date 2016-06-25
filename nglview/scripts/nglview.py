#!/usr/bin/env python
from __future__ import absolute_import
import os, sys, argparse, json
import subprocess
from .cmd_example import CMD_EXAMPLE 

bin_path = sys.prefix + '/bin/'

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

In your local machine, run:

    {client_cm}

NOTE: you might want to replace {hostname} by full hostname with domain name

Then open your web browser, copy and paste:
    http://localhost:{port}
"""

def get_remote_port(port=None):
    import os, socket
    from nglview.scripts.app import NGLViewApp
    port = port if port is not None else NGLViewApp().get_port()

    username = os.getlogin()
    hostname = socket.gethostname()
    client_cm = "ssh -NL localhost:{port}:localhost:{port} {username}@{hostname}".format(username=username,
            hostname=hostname,
            port=port)
    print(remote_msg.format(client_cm=client_cm, port=port, hostname=hostname))
    return port

def main(notebook_dict=notebook_dict):
    PY3 = sys.version_info[0] == 3
    pyv_full_string = ','.join(str(i) for i in sys.version_info)
    pyv_short_string = str(sys.version_info[0])
    default_jexe = bin_path + 'jupyter'

    parser = argparse.ArgumentParser(description='NGLView: An IPython/Jupyter widget to '
                                     'interactively view molecular structures and trajectories.',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=CMD_EXAMPLE)
    parser.add_argument('command',
            help='command could be a topology filename (.pdb, .mol2, .parm7, ...) or \n'
                          'could be a python script (.py), a notebook (.ipynb)') 
    parser.add_argument('-c', '--crd', help='coordinate filename')
    parser.add_argument('--browser', help='web browser')
    parser.add_argument('-j', '--jexe', default=default_jexe, help='jupyter path')
    parser.add_argument('--notebook-name', default='tmpnb_ngl.ipynb', help='notebook name')
    parser.add_argument('--port', type=int, help='port number')
    parser.add_argument('--remote', action='store_true', help='create remote notebook')
    parser.add_argument('--clean-cache', action='store_true', help='delete temp file after closing notebook')
    args = parser.parse_args()

    command = parm = args.command

    crd = args.crd
    if crd is None:
        crd = parm

    browser = '--browser ' + args.browser if args.browser else ''

    create_new_nb = False

    if parm.endswith('.ipynb'):
        notebook_name = parm
    else:
        notebook_name = args.notebook_name
        if not PY3:
            kernelspec = notebook_dict['metadata']['kernelspec']
            kernelspec['display_name'] = 'Python 2'
            kernelspec['name'] = 'python2'
            notebook_dict['metadata']['kernelspec'] = kernelspec

            notebook_dict['metadata']['language_info']['codemirror_mode']['version'] = pyv_short_string
            notebook_dict['metadata']['version'] = pyv_full_string
        if parm.endswith('.py'):
            pycontent = open(parm).read().strip()
            notebook_dict['cells'][0]['source'] = pycontent
            nb_json = json.dumps(notebook_dict)
        else:
            nb_json = json.dumps(notebook_dict)
            nb_json = nb_json.replace('"null"', 'null').replace('test.nc', crd).replace('prmtop', parm)
        nb_json = nb_json.replace('"null"', 'null')

        with open(notebook_name, 'w') as fh:
            fh.write(nb_json)
            create_new_nb = True
    
    
    if not args.remote:
        cm = '{jupyter} notebook {notebook_name} {browser}'.format(jupyter=args.jexe,
                                                                   notebook_name=notebook_name,
                                                               browser=browser)
    else:
        port = get_remote_port(args.port)
        cm = '{jupyter} notebook --no-browser --port {port}'.format(jupyter=args.jexe,
                                                                    port=port)
        print('NOTE: make sure to open {0} in your local machine\n'.format(notebook_name))

    install_nbextension(jupyter=args.jexe)

    try:
        subprocess.check_call(cm.split())
    except KeyboardInterrupt:
        if args.clean_cache and create_new_nb:
            print("deleting {}".format(notebook_name))
            os.remove(notebook_name)
        disable_extension(jupyter=args.jexe)

if __name__ == '__main__':
    main()
