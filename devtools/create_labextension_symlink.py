import os
import sys
import argparse
from pathlib import Path


def get_extension_paths():
    base_dir = Path(sys.prefix)
    nbextension_dest = base_dir / "share" / "jupyter" / "nbextensions" / "nglview"
    labextension_dest = base_dir / "share" / "jupyter" / "labextensions" / "nglview"
    return nbextension_dest, labextension_dest


def create_symlink(source, destination):
    try:
        if os.path.islink(destination):
            os.remove(destination)
        os.symlink(source, destination)
        print(f"Created symlink: {destination} -> {source}")
    except OSError as e:
        print(f"Failed to create symlink: {e}", file=sys.stderr)


def remove_symlink(destination):
    try:
        if os.path.islink(destination):
            os.remove(destination)
            print(f"Removed symlink: {destination}")
        else:
            print(f"No symlink to remove at: {destination}")
    except OSError as e:
        print(f"Failed to remove symlink: {e}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description="Create or remove symlinks for nbextension and labextension.")
    parser.add_argument('-d', '--delete', action='store_true', help="Remove the symlinks instead of creating them.")
    args = parser.parse_args()

    nbextension_dest, labextension_dest = get_extension_paths()

    nbextension_source = Path('nglview') / 'nbextension'
    labextension_source = Path('nglview') / 'labextension'

    if args.delete:
        remove_symlink(nbextension_dest)
        remove_symlink(labextension_dest)
    else:
        create_symlink(nbextension_source, nbextension_dest)
        create_symlink(labextension_source, labextension_dest)


if __name__ == "__main__":
    main()