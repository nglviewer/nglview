import json
from pathlib import Path
import sys
from packaging.version import parse as parse_version

version = sys.argv[1]
HERE = (Path(__file__).parent / '../').resolve().absolute()
sys.path.insert(0, str(HERE))

import versioneer
dirty_version = versioneer.get_version()
clean_version = dirty_version.split('+')[0]

print(f"Current nglview version = {dirty_version}")
print(f"Front end version = {version}")

if parse_version(clean_version) != parse_version(version):
    print("WARNING: Current nglview version should be equal to front end version")

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
