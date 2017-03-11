__all__ = ['linear']


def lerp(a, b, t):
    return (b - a) * t + a


def linear(index, t, traj, step=1):
    """
    
    Parameters
    ----------
    index : int
    t : float, between 0. and 1.
    traj : nglview.Trajectory or subclass
    step : int, default=1
    """

    # need to copy coordinates to avoid early memory free
    # in pytraj
    c = traj.get_coordinates(index).copy()
    cp = traj.get_coordinates(min(index + step, traj.n_frames - 1)).copy()

    coords = lerp(cp, c, t)
    return coords
