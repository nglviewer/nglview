from ipywidgets import DOMWidget
from traitlets import Bool


class BaseWidget(DOMWidget):
    _msg_q = []
    _msg_ar = []
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
        self._msg_ar.append(msg)
