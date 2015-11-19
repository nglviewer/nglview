#!/usr/bin/env python

from setuptools import setup

if __name__ == '__main__':
    setup(
        name = "nglview",
        author = "Alexander S. Rose",
        description = "An IPython widget to interactively view molecular structures and trajectories.",
        version = "0.2dev",
        license = "MIT",
        url = "https://github.com/arose/ngl-widget",
        zip_safe = False,
        package_data = { "nglview.html": [ "static/*" ] },
        packages = [ "nglview", "nglview.html" ],
        extras_require = {
            "simpletraj": [ "simpletraj" ],
            "mdtraj": [ "mdtraj" ]
        }
    )
