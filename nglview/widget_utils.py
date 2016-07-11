
def get_widget_by_name(box, myid):
    for widget in box.children:
        if hasattr(widget, '_ngl_name') and widget._ngl_name == myid:
            return widget
