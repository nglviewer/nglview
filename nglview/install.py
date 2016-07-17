import warnings
from notebook.nbextensions import install_nbextension
from notebook.services.config import ConfigManager

import numpy as np

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from pkg_resources import resource_filename

def install(user=True, symlink=False):
    """Install the widget nbextension.

    Parameters
    ----------
    user: bool
        Install for current user instead of system-wide.
    symlink: bool
        Symlink instead of copy (for development).
    """
    staticdir = resource_filename('nglview', 'js')
    install_nbextension(staticdir, destination='nglview',
                        user=user, symlink=symlink,
                        verbose=0)

    cm = ConfigManager()
    cm.update('notebook', {
        "load_extensions": {
            "widgets/notebook/js/extension": True,
        }
    })

def enable_nglview_js():
    # place holder for ipywidget >= 5.1
    pass
