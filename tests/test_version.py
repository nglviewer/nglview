import json
from pathlib import Path

import nglview as nv


def test_make_sure_versions_matched():
    nglview_dir = Path(__file__).resolve().parent.parent / 'nglview'

    # frontend and backend
    json_str = (nglview_dir / "staticlab/package.json").read_text()
    assert json.loads(json_str)['version'] == nv.widget.__frontend_version__
