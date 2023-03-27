"""Use `pip install` rather than `python setup.py install`
to correctly distribute the notebook and lab extensions.
"""

import os
import subprocess
import sys
from distutils import log
from subprocess import check_call
from pathlib import Path

from setuptools import Command, find_packages, setup
from setuptools.command.egg_info import egg_info

from jupyter_packaging import (
    create_cmdclass,
)

import versioneer
from versioneer import get_cmdclass

sdist = get_cmdclass()['sdist']
build_py = get_cmdclass()['build_py']

here = os.path.dirname(os.path.abspath(__file__))
node_root = os.path.join(here, 'js')
is_repo = os.path.exists(os.path.join(here, '.git'))

npm_path = os.pathsep.join([
    os.path.join(node_root, 'node_modules', '.bin'),
                os.environ.get('PATH', os.defpath),
])

log.set_verbosity(log.DEBUG)
log.info('setup.py entered')
log.info('$PATH=%s' % os.environ['PATH'])


try:
    sys.argv.remove('--npm')
    rebuild_nglview_js = True
except ValueError:
    rebuild_nglview_js = False


def js_prerelease(command, strict=False):
    """decorator for building minified js/css prior to another command"""
    class DecoratedCommand(command):
        def run(self):
            jsdeps = self.distribution.get_command_obj('jsdeps')
            if not is_repo and all(os.path.exists(t) for t in jsdeps.targets):
                # sdist, nothing to do
                command.run(self)
                return

            try:
                self.distribution.run_command('jsdeps')
            except:
                missing = [t for t in jsdeps.targets if not os.path.exists(t)]
                if strict or missing:
                    log.warn('rebuilding js and css failed')
                    if missing:
                        log.error('missing files: %s' % missing)
                    raise e
                else:
                    log.warn('rebuilding js and css failed (not a problem)')
                    log.warn(str(e))
            command.run(self)
            update_package_data(self.distribution)
    return DecoratedCommand

def update_package_data(distribution):
    """update package_data to catch changes during setup"""
    build_py = distribution.get_command_obj('build_py')
    # distribution.package_data = find_package_data()
    # re-init build_py options which load package_data
    build_py.finalize_options()


class NPM(Command):
    # FIXME: remove this class
    description = 'install package.json dependencies using npm'

    user_options = []

    node_modules = os.path.join(node_root, 'node_modules')

    targets = [
        os.path.join(here, 'nglview', 'static', 'extension.js'),
        os.path.join(here, 'nglview', 'static', 'index.js')
    ]

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def has_npm(self):
        try:
            check_call(['npm', '--version'])
            return True
        except:
            return False

    def should_run_npm_install(self):
        package_json = os.path.join(node_root, 'package.json')
        node_modules_exists = os.path.exists(self.node_modules)
        return rebuild_nglview_js

    def run(self):
        has_npm = self.has_npm()
        env = os.environ.copy()
        env['PATH'] = npm_path

        if self.should_run_npm_install():
            if not has_npm:
                log.error("`npm` unavailable.  If you're running this command using sudo, make sure `npm` is available to sudo")
                sys.exit()

            log.info("Installing build dependencies with npm.  This may take a while...")
            check_call(['npm', 'install'], cwd=node_root, stdout=sys.stdout, stderr=sys.stderr)
            os.utime(self.node_modules, None)

        for t in self.targets:
            if not os.path.exists(t):
                msg = 'Missing file: %s' % t
                if not has_npm:
                    msg += '\nnpm is required to build a development version of widgetsnbextension'
                raise ValueError(msg)

        # update package data in case this created new files
        update_package_data(self.distribution)

HERE = Path(__file__).parent.resolve()
# The name of the project
name = "nglview-js-widgets"
lab_path = (HERE / "nglview"/ "staticlab")
package_data_spec = {
    name: ["*"],
}

labext_name = "nglview-js-widgets"
data_files_spec = [
    ("share/jupyter/labextensions/%s" % labext_name, str(lab_path), "**"),
    ("share/jupyter/labextensions/%s" % labext_name, str(HERE), "install.json"),
]

cmdclass = create_cmdclass("jsdeps",
    package_data_spec=package_data_spec,
    data_files_spec=data_files_spec
)
cmdclass['jsdeps'] = NPM
cmdclass['version'] = get_cmdclass()['version']
cmdclass['build_py'] = js_prerelease(build_py)
cmdclass['sdist'] = js_prerelease(sdist, strict=True)
cmdclass['egg_info'] = js_prerelease(egg_info)

setup_args = {
    'name': 'nglview',
    'version': versioneer.get_version(),
    'description': 'IPython widget to interactively view molecular structures and trajectories.',
    'include_package_data': True,
    'license': "MIT",
    'package_data': {
         "nglview.datafiles": ["*"],
         "nglview.scripts": ["*"],
         "nglview.theme": ["*"],
         "nglview.static": ["*"],
         "nglview.staticlab": ["*"],
     },
    'entry_points': {'console_scripts':
          ['nglview = nglview.scripts.nglview:main',]
    },
    'data_files': [
        ('share/jupyter/nbextensions/nglview-js-widgets', [
         'nglview/static/extension.js',
         'nglview/static/index.js',
         'nglview/static/index.js.map',
        ]),
        ('etc/jupyter/nbconfig/notebook.d' , ['nglview-js-widgets.json'])
    ],
    'tests_require': [
        'pytest'
    ],
    'install_requires': [
        'ipywidgets>=7',
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
    'packages': set(find_packages() + 
                ['nglview',
                 'nglview.static',
                 'nglview.staticlab',
                 'nglview.theme',
                 'nglview.datafiles',
                 'nglview.utils',
                 'nglview.tests',
                 'nglview.sandbox',
                 'nglview.contrib',
                 'nglview.scripts']),
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
    'python_requires': '>=3.6',
    'classifiers': [
        'Development Status :: 4 - Beta',
        'Framework :: IPython',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: JavaScript',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
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
