#!/usr/bin/env python

import subprocess
import sys
from pathlib import Path

HERE = (Path(__file__).parent / '../').resolve().absolute()
sys.path.insert(0, str(HERE))

import nglview

print(nglview)

front_end_version = nglview.widget.__frontend_version__
latest_tag = sorted(subprocess.check_output('git tag', shell=True).decode().split('\n'))[-1]
latest_tag = latest_tag.replace('v', '')
print("latest_tag", latest_tag)
print("front_end_version", front_end_version)

print("\nMake sure to publish npm package")
output = subprocess.check_output(["npm", "search", "nglview-js-widgets"]).decode()
print(output)

if front_end_version != latest_tag:
    print(f"Version mismatch between front_end_version {front_end_version} and latest_tag {latest_tag}")
    sys.exit(1)
subprocess.check_call('python setup.py sdist', shell=True)
