#!/usr/bin/env python

import argparse
import os
from glob import glob
from pathlib import Path

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--basic', action='store_true')
    parser.add_argument('--api', action='store_true')
    parser.add_argument('--single', action='store_true')
    parser.add_argument('-n', '--nb', nargs='?')
    return parser.parse_args()

def get_notebooks(args):
    api_root_dir = 'tests/notebooks/api/'
    if args.api:
        return glob(os.path.join(api_root_dir, 'test*.ipynb'))
    elif args.single:
        return [args.nb]
    else:
        return glob('tests/notebooks/test*ipynb') + glob('tests/notebooks/api/test*ipynb')

def get_cell_length(nb):
    n_cells = 0
    with open(nb) as fh:
        for line in fh:
            if 'cell_type' in line:
                n_cells += 1
    return n_cells

def generate_js(notebooks_with_cell_lengths):
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
    all_notebooks = '\n'.join(
        body_template % (notebook, notebook, n_cells)
        for (notebook, n_cells) in notebooks_with_cell_lengths)

    js_dir = Path('tests/js')
    js_dir.mkdir(parents=True, exist_ok=True)

    fn = js_dir / 'test.js'
    fn.write_text(head + all_notebooks + tail)

if __name__ == '__main__':
    args = parse_arguments()
    notebooks = get_notebooks(args)
    notebooks_with_cell_lengths = [(nb, 2 * get_cell_length(nb)) for nb in notebooks]
    generate_js(notebooks_with_cell_lengths)
