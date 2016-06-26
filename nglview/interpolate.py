def lerp(a, b, t):
    return (b - a) * t + a

def _spline(p0, p1, p2, p3, t, tension):
    v0 = ( p2 - p0 ) * tension
    v1 = ( p3 - p1 ) * tension
    t2 = t * t
    t3 = t * t2
    return (( 2 * p1 - 2 * p2 + v0 + v1 ) * t3 +
           ( -3 * p1 + 3 * p2 - 2 * v0 - v1 ) * t2 +
           v0 * t + p1)

doc ="""

    Parameters
    ----------
    index : int
    t : float, between 0. and 1.
    traj : nglview.Trajectory or subclass
    step : int, default=1
"""

def linear(index, t, traj, step=1):
    # need to copy coordinates to avoid early memory free
    # in pytraj
    c = traj.get_coordinates(index).copy()
    cp = traj.get_coordinates(min(index + step, traj.n_frames-1)).copy()

    coords = lerp(cp, c, t)
    return coords

def spline(index, t, traj, step=1):
    i = index
    ip = min(index + step, traj.n_frames - 1)
    ipp = min(index + 2 * step, traj.n_frames - 1)
    ippp = min(index + 3 * step, traj.n_frames - 1)

    c = traj.get_coordinates(i).copy()
    cp = traj.get_coordinates(ip).copy()
    cpp = traj.get_coordinates(ipp).copy()
    cppp = traj.get_coordinates(ippp).copy()

    coords = _spline(cppp, cpp, cp, c, t, tension=1)
    return coords

linear.__doc__ = doc
spline.__doc__ = doc
