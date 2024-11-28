"""Use `pip install` rather than `python setup.py install`
to correctly distribute the notebook and lab extensions.
"""

import os
from distutils import log
from pathlib import Path

from setuptools import find_packages, setup
from jupyter_packaging import wrap_installers, get_data_files

here = os.path.dirname(os.path.abspath(__file__))
node_root = os.path.join(here, 'js')
is_repo = os.path.exists(os.path.join(here, '.git'))

log.set_verbosity(log.DEBUG)
log.info('setup.py entered')
log.info('$PATH=%s' % os.environ['PATH'])


def update_package_data(distribution):
    """update package_data to catch changes during setup"""
    build_py = distribution.get_command_obj('build_py')
    build_py.finalize_options()


HERE = Path(__file__).parent.resolve()
# The name of the project
lab_path = (HERE / "nglview" / "staticlab")
nb_path = (HERE / "nglview" / "static")
assert (nb_path/"index.js").exists(), "index.js not found in %s" % nb_path
assert (lab_path/"static/package.json").exists(), "package.json not found in %s" % lab_path
labext_name = "nglview-js-widgets"
package_data_spec = {
    labext_name: ["*"],
}

data_files_spec = [
    ("share/jupyter/labextensions/%s" % labext_name, str(lab_path), "**"),
    ("share/jupyter/labextensions/%s" % labext_name, str(HERE), "install.json"),
    ("share/jupyter/nbextensions/%s" % labext_name, str(nb_path), "**"),
    ("etc/jupyter/nbconfig/notebook.d", str(HERE), "nglview-js-widgets.json"),
]

def pre_develop():
    pass

def pre_dist():
    pass

cmdclass = wrap_installers(pre_develop=pre_develop, pre_dist=pre_dist)
data_files = get_data_files(data_files_spec)

setup_args = {
    'name': 'nglview',
    "use_scm_version": True,
    "setup_requires": ['setuptools_scm'],
    'description': 'IPython widget to interactively view molecular structures and trajectories.',
    'long_description': open('README.md').read(),
    'long_description_content_type': 'text/markdown',
    'license': "MIT",
    'package_data': {
         "nglview.datafiles": ["*"],
         "nglview.theme": ["*"],
         "nglview.static": ["*"],
         "nglview.staticlab": ["*"],
     },
    'data_files': data_files,
    'install_requires': [
        'ipywidgets>=8',
        'notebook>=7',
        'jupyterlab>=3',
        'jupyterlab_widgets',
        'numpy',
    ],
    'extras_require': {
        "simpletraj": ["simpletraj"],
        "mdtraj": ["mdtraj"],
        "pytraj": ["pytraj"],
        "MDAnalysis": ["MDAnalysis"],
        "ParmEd": ["parmed"],
        "rdkit": ["rdkit"],
        "ase": ["ase"],
        "htmd": ["htmd"],
        "qcelemental": ["qcelemental"],
    },
    'packages': find_packages(),
    'zip_safe': False,
    'cmdclass': cmdclass,
    'author': 'Alexander S. Rose, Hai Nguyen',
    'author_email': 'alexander.rose@weirdbyte.de',
    'url': 'https://github.com/arose/nglview',
    'keywords': [
        'ipython',
        'jupyter',
        'widgets',
    ],
    'python_requires': '>=3.7',
    'classifiers': [
        'Framework :: IPython',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: JavaScript',
        'Programming Language :: Python :: 3 :: Only',
        'Topic :: Multimedia :: Graphics',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Visualization',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Operating System :: MacOS',
    ],
}

setup(**setup_args)
