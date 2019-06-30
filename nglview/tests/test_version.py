import json
import os

import nglview as nv


def test_make_sure_versions_matched():
    nglview_dir = os.path.join(os.path.dirname(__file__), '..')

    # frontend and backend
    nglview_js_widgets = os.path.join(nglview_dir, 'static', 'index.js')
    with open(nglview_js_widgets) as fh:
        for line in fh:
            if '"name":"nglview-js-widgets"' in line:
                break
    line = line.replace('module.exports =', '')
    assert json.loads(line)['version'] == nv.widget.__frontend_version__
