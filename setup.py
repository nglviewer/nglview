#!/usr/bin/env python

from setuptools import setup
import versioneer

CLASSIFIERS = [
    "Development Status :: 3 - Alpha",
    "Intended Audience :: Science/Research",
    "Intended Audience :: Developers",
    "License :: OSI Approved :: MIT License",
    "Programming Language :: JavaScript",
    "Programming Language :: Python",
    "Programming Language :: Python :: 2",
    "Programming Language :: Python :: 3",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "Topic :: Scientific/Engineering :: Chemistry",
    "Topic :: Scientific/Engineering :: Visualization",
    "Operating System :: Microsoft :: Windows",
    "Operating System :: POSIX",
    "Operating System :: Unix",
    "Operating System :: MacOS",
]


if __name__ == '__main__':
    setup(
        name="nglview",
        author="Alexander S. Rose",
        author_email="alexander.rose@weirdbyte.de",
        description="IPython widget to interactively view molecular\
                     structures and trajectories.",
        version=versioneer.get_version(),
        cmdclass=versioneer.get_cmdclass(),
        classifiers=CLASSIFIERS,
        license="MIT",
        url="https://github.com/arose/nglview",
        zip_safe=False,
        package_data={
            "nglview.datafiles": ["*"],
            "nglview.js": ["*"],
            "nglview.scripts": ["*"]
        },
        packages=["nglview", "nglview.datafiles", "nglview.js", "nglview.scripts"],
        install_requires=["jupyter", "traitlets>=4.2.1", "ipywidgets<5.0"],
        tests_require=["nose"],
        test_suite="nose.collector",
        extras_require={
            "simpletraj": ["simpletraj"],
            "mdtraj": ["mdtraj"],
            "pytraj": ["pytraj"],
            "MDAnalysis": ["MDAnalysis"],
            "ParmEd": ["parmed"],
            "rdkit": ["rdkit"],
        },
      entry_points={'console_scripts':
          ['nglview = nglview.scripts.nglview:main',]
           },
    )
