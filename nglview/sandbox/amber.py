import os
import threading
import time

import pytraj as pt

from .base import BaseMD


class AmberMD(BaseMD):
    # TODO: doc
    '''
    Unstable API

    Examples
    --------
    >>> from nglview.sandbox.amber import AmberMD
    >>> amber_view = AmberMD(top='./peptide.top', restart='./md.r', reference='min.rst7')
    >>> view = amber_view.initialize()
    >>> view
    >>> # another cell
    >>> amber_view.update(every=1, timeout=3000)
    >>> # do other stuff
    '''

    def __init__(self, top=None, restart=None, reference=None):
        self.top = top
        assert os.path.exists(restart), f'{restart} must exists'
        assert os.path.exists(reference), f'{reference} must exists'
        self.restart = restart
        self.reference_traj = pt.load(reference, top=self.top)
        self.thread = None
        self.event = None

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
        # always reset
        self.event = threading.Event()

        def _update(event):
            start = time.time()
            while time.time() - start <= timeout and not event.is_set():
                time.sleep(every)
                traj = pt.load(self.restart, top=self.top)
                if callback is not None:
                    callback(traj)
                else:
                    mask = '@C,N,O,CA,P'
                    pt.superpose(traj, mask=mask, ref=self.reference_traj)
                self._update_coordinates(traj[0].xyz)

        # non-blocking so we can use other Jupyter's cells
        self.thread = threading.Thread(target=_update, args=(self.event, ))
        self.thread.daemon = True
        self.thread.start()

    def stop(self):
        """Stop update"""
        if self.event is not None:
            self.event.set()

    def _update_coordinates(self, xyz):
        self.view.coordinates_dict = {0: xyz}
