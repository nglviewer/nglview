#!/usr/bin/env python

import os
import subprocess
from glob import glob
from random import shuffle

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--basic', action='store_true')
parser.add_argument('--api', action='store_true')
parser.add_argument('--single', action='store_true')
parser.add_argument('--travis', action='store_true')
parser.add_argument('-n', '--nb', nargs='?')
args = parser.parse_args()

api_root_dir = 'nglview/tests/notebooks/api/'

notebook_names_for_travis  = [
             'test_no_gui_demo.ipynb',
             'test_add_structure_then_trajectory.ipynb',
             'test_automatically_added_attributes_0.ipynb',
            ]

if args.travis:
    notebook_names = notebook_names_for_travis
    notebooks = [api_root_dir + notebook_name for notebook_name in notebook_names]
elif args.api:
    notebooks = glob(api_root_dir + '/*.ipynb')
elif args.single:
    notebooks = [args.nb]
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
                  .pause(3000)
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

    for nb_abspath, n_cells in notebooks_with_cell_lengths:
        basename = os.path.basename(nb_abspath)
        js_fn = os.path.join('nglview/tests/js',
                os.path.splitext(basename)[0] + '.js')
        with open(js_fn, 'w') as fh:
            body = body_template % (nb_abspath, nb_abspath, n_cells)
            fh.write(head + body + tail)
