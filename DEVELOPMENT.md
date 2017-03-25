
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
# or cd js && npm install

# known working versions:
# node: v6.6.0
# npm: 4.4.4


# then python setup.py install or pip install -e . # development, you can edit the source code without re-installing

# tips
# - Use private browser mode to avoid cache
# - After changing js code:
#     - Kernel --> Restart and Clear Output
#     - Refresh webpage (F5)

# run quick test in terminal
nglview demo
```

Using `NGL` locally
===================

1. Change 
`var NGL = require('ngl');` to `var NGL = require('./ngl');`
https://github.com/arose/nglview/blob/master/js/src/widget_ngl.js#L2

2. Then, [build NGL](https://github.com/arose/ngl/blob/master/DEVELOPMENT.md#building), then copy `ngl.js` (or `ngl.dev.js`) to `nglview/js/src/`

3. Rebuild js code
```
cd nglview (root folder having setup.py file)
python setup.py install --npm
```

You need to install `nodejs` (which includes `npm`).
Tips: `conda install nodejs -c conda-forge` (and so on)

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
