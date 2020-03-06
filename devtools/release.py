#!/usr/bin/env python

import subprocess
import sys

import nglview

print(nglview)

front_end_version = nglview.widget.__frontend_version__
latest_tag = sorted(subprocess.check_output('git tag', shell=True).decode().split('\n'))[-1]
latest_tag = latest_tag.replace('v', '')
print("latest_tag", latest_tag)
print("front_end_version", front_end_version)

if front_end_version != latest_tag:
    print(f"Version mismatch between front_end_version {front_end_version} and latest_tag {latest_tag}")
    sys.exit(1)
subprocess.check_call('python setup.py sdist', shell=True)
print("Make sure to publish npm package")
