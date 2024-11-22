#!/usr/bin/env python

import json
from pathlib import Path
import sys
from packaging.version import parse as parse_version

version = sys.argv[1]
HERE = Path(__file__).parents[1].resolve().absolute()
pkg_file = HERE / "js/package.json"

content = []

with open(pkg_file) as fh:
    for line in fh:
        if line.strip().startswith('"version"'):
            content.append(f'  "version": "{version}",\n')
        else:
            content.append(line)

with open(pkg_file, 'w') as fh:
    for line in content:
        fh.write(line)


with open('nglview/_frontend.py', 'w') as fh:
    fh.write(f"__frontend_version__ = '{version}'")


print("DONE updating version")
