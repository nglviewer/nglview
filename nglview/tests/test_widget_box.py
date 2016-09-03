import nglview as nv

def test_widget_box():
    # empty
    box = nv.widget_box.BoxNGL()
    box._update_size()
    view = nv.demo()
    box = nv.widget_box.BoxNGL([view])
    box._update_size()

    box._is_beautified = True
    box._beautify()
    box._is_beautified = False
    box._beautify()
