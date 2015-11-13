#!/usr/bin/env python

from setuptools import setup

if __name__ == '__main__':
    setup(
        name = "nglview",
        author = "Alexander S. Rose",
        description = "An IPython widget to interactively view molecular structures and trajectories.",
        version = "0.1",
        license = "MIT",
        url = "https://github.com/arose/ngl-widget",
        package_data = { "nglview.html": [ "static/*" ] },
        packages = [ "nglview", "nglview.html" ],
        # install_requires = [ "simpletraj" ],
    )
