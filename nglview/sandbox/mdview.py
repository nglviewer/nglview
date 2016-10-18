import os
import time
import pytraj as pt
from threading import Thread

class AmberMD(object):
    # TODO: doc
    """

    Examples
    --------
    >>> from nglview.sandbox.mdview import AmberMD
    >>> ambermd = AmberMD(top='./peptide.top', restart='./md.r', reference='min.rst7')
    >>> view = ambermd.initialize()
    >>> view
    >>> # another cell
    >>> ambermd.update(every=1, timeout=3000)
    >>> # do other stuff
    """
    def __init__(self, top=None, restart=None, reference=None):
        self.top = top
        assert os.path.exists(restart), '{} must exists'.format(restart)
        assert os.path.exists(reference), '{} must exists'.format(reference)
        self.restart = restart
        self.reference_traj = pt.load(reference, top=self.top)
        self.thread = None

    def initialize(self, gui=False):
        self.view = self.reference_traj.visualize(gui=gui)
        return self.view

    def update(self, every=1, autoimage=False, timeout=3000):

        def _update():
            start = time.time()
            while time.time() - start <= timeout:
                time.sleep(every)
                if autoimage:
                    traj.autoimage()
                traj = pt.load(self.restart, top=self.top)
                mask = '@C,N,O,CA,P'
                pt.superpose(traj, mask=mask, ref=self.reference_traj)
                self.view.coordinates_dict = {0: traj[0].xyz}
        # non-blocking so we can use other Jupyter's cells
        self.thread = Thread(target=_update)
        self.thread.daemon = True
        self.thread.start()
