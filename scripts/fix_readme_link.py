#!/usr/bin/env python
import os
import sys
import subprocess

readme_file = sys.argv[1]

if os.path.abspath(readme_file) == os.path.abspath('./README'):
    print("wrong README file {}".format(readme_file))
    sys.exit(1)

readme = open(sys.argv[1]).read()

words = ["(nglview.gif)",
         "(examples/images/membrane.gif)",
         "(examples/README.md)",
         "(examples/README.md)",
         "(examples/images/nglview_gui.png)",
         ]

for word in words:
    if '.gif' in word:
        new_word = "(" + 'https://github.com/arose/nglview/blob/master/' + word.strip('(').strip(')') + "?raw=true)"
    else:
        new_word = "(" + 'https://github.com/arose/nglview/blob/master/' + word.strip('(')
    if './' in new_word:
        new_word = new_word.replace('./', '')
    print(new_word)
    readme = readme.replace(word, new_word)

readme = readme.replace('doc/interface_classes.md', 'interface_classes.html')
readme = readme.replace('CHANGELOG.md', 'CHANGELOG.html')

with open('doc/index.md', 'w') as fh:
    fh.write(readme)

subprocess.check_call('pandoc doc/index.md -o doc/index.rst', shell=True)
subprocess.check_call('pandoc ../nglview/doc/interface_classes.md  -o doc/interface_classes.rst', shell=True)
subprocess.check_call('pandoc ../nglview/CHANGELOG.md  -o doc/CHANGELOG.rst', shell=True)
