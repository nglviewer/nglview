#!/usr/bin/env python

import subprocess
from glob import glob
from random import shuffle

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--basic', action='store_true')
args = parser.parse_args()

if args.basic:
    root = 'nglview/tests/notebooks/api/'
    notebook_names  = [
                 # 'test_auto_detect_pytraj_mdtraj_mdanalysis_parmed.ipynb',
                 'test_view_trajectory.ipynb',
                 'test_no_gui_demo.ipynb',
                 'test_add_structure_then_trajectory.ipynb',
                 'test_automatically_added_attributes_0.ipynb',
                 'test_fix_player_if_adding_single_struture_first.ipynb',
                ]

    notebooks = [root + notebook_name for notebook_name in notebook_names]
else:
    notebooks = ['nglview/tests/notebooks/dummy.ipynb',]
    
    notebooks += (glob('nglview/tests/notebooks/*ipynb') +
                glob('nglview/tests/notebooks/api/*ipynb'))

# shuffle(notebooks)
def get_cell_length(nb):
    n_cells = 0

    with open(nb) as fh:
        for line in fh.readlines():
            if 'cell_type' in line:
                n_cells += 1
    return n_cells

notebooks_with_cell_lengths = [(nb, 2*get_cell_length(nb)) for nb in notebooks]

head = """
module.exports = {

"""


body_template = """
    "%s": function (browser) {
        browser.openNotebook("%s");
        browser.restartKernel(2000);
        for ( var i = 0; i < %s; i++) {
           browser.executeCell(i)
                  .pause(2500)
                  .cellHasError(i);
        }
        browser.end();
    },
"""

tail = """
}
"""

if __name__ == '__main__':

    all_notebooks = '\n'.join(body_template % (notebook, notebook, n_cells)
                              for (notebook, n_cells) in notebooks_with_cell_lengths)

    fn = 'nglview/tests/js/test.js'
    with open(fn, 'w') as fh:
        fh.write(head + all_notebooks + tail)
