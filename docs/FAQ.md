### Can not import nglview although successfully installed it?

> You can try

> ```bash
> python -m ipykernel install --user
> ```

> Then in your Jupyter notebook, choose the right `kernel`. If you are using `python 2`, make sure to choose `Python 2` kernel.

### widget not shown?

- Could not cross validate the widget frontend and backend versions (or similiar)

Double check if you are having two ipywidgets version (e.g: one installed via pip and one installed via conda)

- Class NGLModel not found in module nglview-js-widgets

You are likely using older JavaScript distribution of nglview. Check if it is 
`./.local/share/jupyter/nbextensions/nglview-js-widgets/`, if Yes, delete it.

### Can I have two MDA.Atomgroups in the same view?

https://github.com/arose/nglview/issues/434
