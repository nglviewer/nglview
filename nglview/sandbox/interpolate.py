"""This module is meant for playground.
"""

import numpy as np


def smooth(coords_3d,
           method='filtfilt',
           inplace=True,
           atom_indices=None,
           **kwargs):
    """
    Parameters
    ----------
    coords_3d : 3d numpy's array
    method : str, {'filtfilt', 'savgol_filter'}
    atom_indices : None or 1d integer array-like
        If None, using all atoms.
    kwargs : Additional arguments for corresponding method.

    Examples
    --------
    >>> import nglview as nv
    >>> from nglview.sandbox.interpolate import smooth
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2()[:]
    >>> traj.xyz = smooth(traj.xyz)
    >>> view = nv.show_pytraj(traj)
    >>> view.clear()
    >>> view.add_licorice('protein')
    >>> view
    """
    from scipy import signal
    if inplace:
        xyz = coords_3d
    else:
        xyz = coords_3d.copy()
    n_atoms = coords_3d[0].shape[0]
    xyz0 = xyz[0].copy()
    atom_indices = atom_indices or range(n_atoms)
    for idx in atom_indices:
        for j in range(3):
            xyz_idx_j = xyz[:, idx, j]
            if method == 'filtfilt':
                b, a = signal.butter(3, 0.1)
                smooth_data = signal.filtfilt(b, a, xyz_idx_j, **kwargs)
            elif method == 'savgol_filter':
                smooth_data = signal.savgol_filter(xyz_idx_j, **kwargs)
            else:
                raise ValueError(f"%s method is not supported" % method)
            xyz[:, idx, j] = smooth_data
    xyz[0] = xyz0  # to correctly render 1st frame.
    return xyz
