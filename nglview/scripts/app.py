import socket
from itertools import chain

from notebook.notebookapp import NotebookApp, random_ports
from tornado import httpserver


class NGLViewApp(NotebookApp):
    def get_port(self, port=None):

        # user-given port
        uport = [
            port,
        ] if port is not None else [
            self.port,
        ]

        def func():
            return

        self.http_server = httpserver.HTTPServer(func)
        success = None

        for port in chain(uport, random_ports(self.port,
                                              self.port_retries + 1)):
            try:
                self.http_server.listen(port, self.ip)
            except (OSError, socket.error):
                pass
            else:
                self.port = port
                success = True
                break
        if not success:
            self.log.critical(
                'ERROR: the notebook server could not be started because '
                'no available port could be found.')
            self.exit(1)
        self.http_server.stop()
        return self.port
