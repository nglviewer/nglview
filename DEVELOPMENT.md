
General
=======

* Please follow [pep8](https://www.python.org/dev/peps/pep-0008/) to format Python code.


Making a Release
================

We use [Versioneer](https://github.com/warner/python-versioneer) to automatically update the version string (of a release but also in development). This means for a release a new git tag should be created. The tag should be of the form vX.Y or vX.Y.Z and generally follow [pep440](https://www.python.org/dev/peps/pep-0440/) with a prefixed "v".

```bash
git tag  -a vX.Y -m "version X.Y"
git push
git push origin --tags
python setup.py sdist upload -r pypi  # better use twine for uploading, see below
```

To ensure a secure upload use `twine`:
```bash
# Create some distributions in the normal way:
python setup.py sdist
# Upload with twine:
twine upload dist/*
```
Install
=======
Install developer mode so you don't need to  re-install after chaning your code.

```python
pip install -e .
```

Update Javascript build
========================
```bash
# install nodejs
# conda install -c javascript nodejs

# if you want to update ngl code
# cp /dir/to/ngl/dist/ngl.js js/src/

# build
python setup.py build --npm

# then python setup.py install or pip install -e . # development

# you can install `watch` to autobuild if the js source code is changed.
npm install watch -g
watch 'python setup.py build --npm' js/src/
# tips
# - Use private browser mode to avoid cache
# - After changing js code:
#     - Kernel --> Restart and Clear Output
#     - Refresh webpage (F5)

# run quick test in terminal
nglview demo
```

Test notebook
=============

- [edit to add more notebooks] and update notebook files
```bash
python ./devtools/make_test_js.py
```

- run

```bash
source devtools/travis-ci/clone_nbtest.sh # only once
jupyter notebook --port=8889 &
# npm install -g nightwatch
nightwatch
```

More stuff
==========

[wiki](https://github.com/arose/nglview/wiki)
