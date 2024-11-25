import json
from pathlib import Path

import nglview as nv


def test_make_sure_versions_matched():
    nglview_dir = Path(__file__).resolve().parent.parent / 'nglview'

    # frontend and backend
    nglview_js_widgets = nglview_dir / 'static' / 'index.js'
    with nglview_js_widgets.open() as fh:
        for line in fh:
            if '"name":"nglview-js-widgets"' in line:
                break
    line = line.replace('module.exports =', '')
    assert json.loads(line)['version'] == nv.widget.__frontend_version__
