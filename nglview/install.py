import argparse
from os.path import dirname, abspath, join
from notebook.nbextensions import install_nbextension, enable_nbextension


def install(user=True, symlink=False, overwrite=True, **kwargs):
    """Install nglview nbextension.

    Parameters
    ----------
    user: bool, default True
        Install for current user instead of system-wide.
    symlink: bool, default False
        Symlink instead of copy (for development).
    overwrite: bool, default True
        Overwrite previously-installed files for this extension
    **kwargs: keyword arguments
        Other keyword arguments passed to the install_nbextension command
    """
    directory = join(dirname(abspath(__file__)), 'static')
    install_nbextension(directory, destination='nglview',
                        symlink=symlink, user=user, overwrite=overwrite,
                        **kwargs)


def enable_nglview_js(user=True):
    enable_nbextension('nglview', '', user=user)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="nglview")
    parser.add_argument("-u", "--user",
                        help="Install as current user instead of system-wide",
                        action="store_true")
    parser.add_argument("-s", "--symlink",
                        help="Symlink instead of copying files",
                        action="store_true")
    parser.add_argument("-f", "--force",
                        help="Overwrite any previously-installed files for this extension",
                        action="store_true")
    args = parser.parse_args()
    install(user=args.user, symlink=args.symlink, overwrite=args.force)
    enable_nglview_js()
