#!/usr/bin/env python
import os, sys, argparse, json
import subprocess

remote_msg = """
SSH port forwarding help
    In your remote machine
    ----------------------
    
        export port=8890 # if port=8890 is not available, pick another one
        jupyter notebook --port=$port --no-browser
    
    In your local machine
    ---------------------
    
        export port=8890 #  same as given port in your remote machine
        ssh -N -f -L localhost:$port:localhost:$port {username}@{hostname}
    
    Then open your favorite web browser, paste
    
        localhost:8890
    
        # Note: change 8890 to the port number you specified
Troubleshooting:
    If you get 'bind: Address already in use', please issue another port number
"""

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

def help_remote(remote_msg=remote_msg):
    import os, socket
    username = os.getlogin()
    hostname = socket.gethostname()

    print(remote_msg.format(username=username, hostname=hostname))

def install_nbextension(jupyter):
    path = os.path.dirname(__file__)
    nglview_main = os.path.join(path, 'nglview_main.js')

    cm_install = '{jupyter} nbextension install {nglview_main}'.format(jupyter=jupyter,
            nglview_main=nglview_main)
    cm_activate = '{jupyter} nbextension enable nglview_main'.format(jupyter=jupyter) 

    subprocess.check_call(cm_install.split())
    subprocess.check_call(cm_activate.split())


def main(notebook_dict=notebook_dict):
    PY3 = sys.version_info[0] == 3
    parser = argparse.ArgumentParser(description='NGLView')
    # parser.add_argument('-p', '--parm', help='Topology filename', required=True)
    parser.add_argument('parm', help='Topology filename (could be PDB, CIF, ... files)') 
    parser.add_argument('-c', '--crd', help='Coordinate filename')
    parser.add_argument('--browser', help='web browser, optional')
    parser.add_argument('-j', '--jexe', default='jupyter', help='jupyter command, optional')
    parser.add_argument('--notebook-name', default='tmpnb_ngl.ipynb', help='notebook name, optional')
    args = parser.parse_args()

    parm = args.parm

    crd = args.crd
    if crd is None:
        crd = parm

    browser = '--browser ' + args.browser if args.browser else ''

    if parm.lower() == 'remote':
        help_remote()
    else:
        if parm.endswith('.ipynb'):
            notebook_name = parm
        else:
            notebook_name = args.notebook_name
            if parm.endswith('.py'):
                pycontent = open(parm).read().strip()
                notebook_dict['cells'][0]['source'] = pycontent
                nb_json = json.dumps(notebook_dict)
            else:
                nb_json = json.dumps(notebook_dict)
                nb_json = nb_json.replace('"null"', 'null').replace('test.nc', crd).replace('prmtop', parm)
            if not PY3:
               nb_json.replace('python3', 'python')

            nb_json = nb_json.replace('"null"', 'null')

            with open(notebook_name, 'w') as fh:
                fh.write(nb_json)
        
        
        cm = '{jupyter} notebook {notebook_name} {browser}'.format(jupyter=args.jexe,
                                                                   notebook_name=notebook_name,
                                                                   browser=browser)
        print(cm)
        install_nbextension(jupyter=args.jexe)
        # subprocess.check_call(cm.split())
        # not using subprocess to avoid nasty message after closing kernel
        os.system(cm)

if __name__ == '__main__':
    main()
