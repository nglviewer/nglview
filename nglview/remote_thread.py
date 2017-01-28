import time
import threading

class RemoteCallThread(threading.Thread):
    def __init__(self, view, timeout=0.1):
        self.q = []
        self.view = view
        self.timeout = timeout
        super(RemoteCallThread, self).__init__()
        self.daemon = True

    def run(self):
        # this `run` method will be called if `start` the thread
        # How does this work?
        # First, try to pop all callbacks and execute them, if loadFile
        # then wait until getting 'ok' signal from NGL. Calling those callbacks
        # will block execution from other threads. This is why we let this thread
        # run forever in background. This thread is usually sleeping all the time
        # and let other threads do the work. It "wake up" only if its 'q' is added more
        # callbacks

        # This class is needed if use call 
        # add_trajectory, clear, add_representation, ... in the same notebook cell
        while True:
            try:
                callback = self.q.pop(0)
                callback(self.view)
                if callback._method_name in ['loadFile']:
                    self.view._wait_until_finished()
            except IndexError:
                time.sleep(self.timeout)

