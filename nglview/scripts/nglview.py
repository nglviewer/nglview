#!/usr/bin/env python
import os, sys, argparse

notebook_content = r"""
{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {
    "collapsed": true
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
   "pygments_lexer": "ipython3",
   "version": "3.5.1"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 0
}
""".strip()

def main(notebook_content=notebook_content):
    parser = argparse.ArgumentParser(description='NGLView')
    # parser.add_argument('-p', '--parm', help='Topology filename', required=True)
    parser.add_argument('parm', help='Topology filename (could be PDB, CIF, ... files)') 
    parser.add_argument('-c', '--crd', help='Coordinate filename')
    parser.add_argument('-j', '--jexe', default='jupyter', help='jupyter command, optional')
    args = parser.parse_args()

    parm = args.parm
    crd = args.crd
    if crd is None:
        crd = parm

    notebook_name = 'tmpnb_ngl.ipynb'
    notebook_content = notebook_content.replace('test.nc', crd).replace('prmtop', parm)

    with open(notebook_name, 'w') as fh:
        fh.write(notebook_content)
    
    
    cm = '{jupyter} notebook {notebook_name}'.format(jupyter=args.jexe, notebook_name=notebook_name)
    print(cm)
    os.system(cm)

if __name__ == '__main__':
    main()
