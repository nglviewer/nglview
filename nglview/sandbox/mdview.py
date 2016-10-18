import os
import time
import pytraj as pt
from threading import Thread
import abc, six

@six.add_metaclass(abc.ABCMeta)
class BaseMD(object):
    @abc.abstractmethod
    def initialize(self):
        pass

    @abc.abstractmethod
    def update(self):
        pass

class AmberMD(BaseMD):
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

    def initialize(self, gui=False):
        self.view = self.reference_traj.visualize(gui=gui)
        return self.view

    def update(self, every=1, timeout=3000, callback=None):
        """

        Parameters
        ----------
        every : int, default 1 (s)
            update coordinates "every" second
        timeout : int, default 3000 (s)
            stop updating coordinate after "timeout"
        callback : func, optional
            If given, trajectory will be processed (autoimage, ...)
            Must follow func(traj)
        """
        def _update():
            start = time.time()
            while time.time() - start <= timeout:
                time.sleep(every)
                traj = pt.load(self.restart, top=self.top)
                if callback is not None:
                    callback(traj)
                else:
                    mask = '@C,N,O,CA,P'
                    pt.superpose(traj, mask=mask, ref=self.reference_traj)
                self._update_coordinates(traj[0].xyz)
        # non-blocking so we can use other Jupyter's cells
        self.thread = Thread(target=_update)
        self.thread.daemon = True
        self.thread.start()

    def _update_coordinates(self, xyz):
        self.view.coordinates_dict = {0: xyz}
