from IPython.display import display


def lab_display(view):
    from sidecar import Sidecar
    with Sidecar(title='nglview'):
        display(view)
