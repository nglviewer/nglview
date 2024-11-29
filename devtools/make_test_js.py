#!/usr/bin/env python

import argparse
import os
import subprocess
from glob import glob
from random import shuffle

parser = argparse.ArgumentParser()
parser.add_argument('--basic', action='store_true')
parser.add_argument('--api', action='store_true')
parser.add_argument('--single', action='store_true')
parser.add_argument('-n', '--nb', nargs='?')
args = parser.parse_args()

api_root_dir = 'tests/notebooks/api/'

if args.api:
    notebooks = glob(os.path.join(api_root_dir, 'test*.ipynb'))
elif args.single:
    notebooks = [args.nb]
else:
    notebooks = []
    notebooks += (glob('tests/notebooks/test*ipynb') +
                  glob('tests/notebooks/api/test*ipynb'))

# shuffle(notebooks)
def get_cell_length(nb):
    n_cells = 0
    with open(nb) as fh:
        for line in fh:
            if 'cell_type' in line:
                n_cells += 1
    return n_cells

notebooks_with_cell_lengths = [(nb, 2 * get_cell_length(nb)) for nb in notebooks]

head = """
module.exports = {
"""

body_template = """
    "%s": function (browser) {
        browser.openNotebook("%s");
        browser.restartKernel(2000);
        for (var i = 0; i < %s; i++) {
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
    from pathlib import Path

    all_notebooks = '\n'.join(body_template % (notebook, notebook, n_cells)
                              for (notebook, n_cells) in notebooks_with_cell_lengths)

    js_dir = Path('tests/js')
    js_dir.mkdir(parents=True, exist_ok=True)

    fn = js_dir / 'test.js'
    fn.write_text(head + all_notebooks + tail)
