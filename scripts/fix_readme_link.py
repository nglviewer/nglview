#!/usr/bin/env python
import sys
import subprocess

readme = open(sys.argv[1]).read()

words = ["(nglview.gif)",
         "(examples/membrane.gif)",
         "(doc/interface_classes.md)",
         "(CHANGELOG.md)",
         "(examples/README.md)",
         ]

for word in words:
    new_word = "(" + 'https://github.com/arose/nglview/' + word.strip('(')
    readme = readme.replace(word, new_word)

with open('doc/index.md', 'w') as fh:
    fh.write(readme)

subprocess.check_call('pandoc doc/index.md -o doc/index.rst', shell=True)
