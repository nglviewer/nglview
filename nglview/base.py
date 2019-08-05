from ipywidgets import DOMWidget
from traitlets import Bool, List


def _singleton(cls):
    # https://www.python.org/dev/peps/pep-0318/#examples
    instances = {}
    def getinstance():
        if cls not in instances:
            instances[cls] = cls()
        return instances[cls]
    return getinstance


class BaseWidget(DOMWidget):
    _msg_q = []
    _msg_ar = List().tag(sync=True)
    _ready = Bool(False).tag(sync=True)

    def _js(self, code):
        self._call("executeCode", code)

    def _call(self, method_name, *args, **kwargs):
        msg = {"type": "callMethod",
               "methodName": method_name,
                "args": args, "kwargs": kwargs}
        if not self._ready:
            # fire later
            self._msg_q.append(msg)
        else:
            self.send(msg)
        msg_ar = self._msg_ar[:]
        msg_ar.append(msg)
        self._msg_ar = msg_ar # trigger sync
