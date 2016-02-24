
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
python setup.py sdist upload -r pypi
```
