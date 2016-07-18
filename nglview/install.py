import argparse
from os.path import dirname, abspath, join
from notebook.nbextensions import install_nbextension, enable_nbextension

def install(user=True, symlink=False, overwrite=True, **kwargs):
    """Install the bqplot nbextension.
    
    Parameters
    ----------
    user: bool
        Install for current user instead of system-wide.
    symlink: bool
        Symlink instead of copy (for development).
    overwrite: bool
        Overwrite previously-installed files for this extension
    **kwargs: keyword arguments
        Other keyword arguments passed to the install_nbextension command
    """
    directory = join(dirname(abspath(__file__)), 'static')
    install_nbextension(directory, destination='nglview',
                        symlink=symlink, user=user, overwrite=overwrite,
                        **kwargs)

def enable_nglview_js():
    enable_nbextension('nglview-js', '')

    
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
