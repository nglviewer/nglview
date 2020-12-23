.. ipywidgets-display::
    from ipywidgets import Text
    Text("hello")

.. ipywidgets-display::

    import nglview as nv
    from IPython.display import display
    import pytraj
    traj = pytraj.datafiles.load_trpcage()
    view = nv.NGLWidget()
    view.add_trajectory(traj)
    display(view)
    view._set_serialization(frame_range=(0, 10))
