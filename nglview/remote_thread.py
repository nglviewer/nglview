import threading
import time


class RemoteCallThread(threading.Thread):
    def __init__(self, view, timeout=0.1, registered_funcs=['loadFile']):
        '''

        Parameters
        ----------
        view : NGLWidget
        timeout : float (second)
            sleep every `timeout`
        registered_funcs : List[str]
            List of funtion names to wait for.
        '''
        self.q = []
        self.view = view
        self.timeout = timeout
        super().__init__()
        self.daemon = True
        self.registered_funcs = registered_funcs

    def run(self):
        '''
        This `run` method will be called if `start` the thread
        How does this work?

        First, try to pop all callbacks and execute them, if loadFile
        then wait until getting 'ok' signal from NGL. Calling those callbacks
        will block execution from other threads. This is why we let this thread
        run forever in background. This thread is usually sleeping all the time
        and let other threads do the work. It "wake up" only if its 'q' is added more
        callbacks.

        This class is needed if use call 
        add_trajectory, clear, add_representation, ... in the same notebook cell
        '''
        while True:
            try:
                callback = self.q.pop(0)
                callback(self.view)
                if callback._method_name in self.registered_funcs:
                    self.view._wait_until_finished()
            except IndexError:
                time.sleep(self.timeout)
