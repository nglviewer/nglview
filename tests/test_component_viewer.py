import nglview as nv


def test_component_viewer():
    view = nv.NGLWidget()
    c = view.add_component(nv.datafiles.PDB)
    c._call("setPosition", [1, 2, 3])
    c.add_cartoon()
    view.remove_component(c)
    assert c._view is None
