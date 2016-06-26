def lerp(a, b, t):
    return (b - a) * t + a

def linear(index, t, traj, step=1):
    """

    Parameters
    ----------
    index : int
    step : int
    t : float, between 0. and 1.
    traj : nglview.Trajectory or subclass
    """
    # need to copy coordinates to avoid early memory free
    # in pytraj
    c = traj.get_coordinates(index).copy()
    cp = traj.get_coordinates(index + step).copy()

    coords = lerp(cp, c, t)
    return coords
