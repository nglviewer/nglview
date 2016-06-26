# translated from ngl code
def _spline(p0, p1, p2, p3, t, tension):

    v0 = (p2 - p0) * tension
    v1 = (p3 - p1) * tension
    t2 = t * t
    t3 = t * t2

    return ((2 * p1 - 2 * p2 + v0 + v1) * t3 +
            (-3 * p1 + 3 * p2 - 2 * v0 - v1) * t2 +
            v0 * t + p1)

def _lerp(a, b, t):
    return a + (b - a) * t

def lerp(index, t, traj, step=1):
    """

    Parameters
    ----------
    index : int
    step : int
    t : float?
    traj : nglview.Trajectory or subclass
    """
    c = traj.get_coordinates(index)
    cp = traj.get_coordinates(index + step)

    coords = _lerp(cp, c, t)
    return coords

def spline(index, t, traj, step=1):

    i, ip, ipp, ippp = index, index + step, index + 2 * step, index + 3 * step
    c = traj.get_coordinates(i)
    cp = traj.get_coordinates(ip)
    cpp = traj.get_coordinates(ipp)
    cppp = traj.get_coordinates(ippp)

    coords = _spline(cppp, cpp, cp, c, t, 1)
    return coords
