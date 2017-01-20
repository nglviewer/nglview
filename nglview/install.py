import argparse
from os.path import join
from notebook import nbextensions
from glob import glob


def install_nglview_js_widgets(user=True, symlink=False, overwrite=True, debug=False, **kwargs):
    """Install nglview-js-widgets nbextension.

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
    nglivew_js_dirs = nbextensions.install_nbextension_python('nglview',
            user=user, symlink=symlink, overwrite=overwrite, **kwargs)
    if debug:
        print(nglivew_js_dirs)
        print([glob(join(my_dir, '*')) for my_dir in nglivew_js_dirs])

def enable_nglview_js_widgets(user=True):
    nbextensions.enable_nbextension_python('nglview', user=user)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="nglview-js-widgets")
    parser.add_argument("-u", "--user",
                        help="Install as current user instead of system-wide",
                        action="store_true")
    parser.add_argument("-s", "--symlink",
                        help="Symlink instead of copying files",
                        action="store_true")
    parser.add_argument("-f", "--force",
                        help="Overwrite any previously-installed files for this extension",
                        action="store_true")
    parser.add_argument("-d", "--debug",
                        help="print nglivew-js-widgets",
                        action="store_true")
    args = parser.parse_args()
    install_nglview_js_widgets(user=args.user,
            symlink=args.symlink,
            overwrite=args.force,
            debug=args.debug)
    enable_nglview_js_widgets()
