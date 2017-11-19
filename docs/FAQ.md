1. Can not import nglview although successfully installed it?

> You can try

> ```bash
> python -m ipykernel install --user
> ```

> Then in your Jupyter notebook, choose the right `kernel`. If you are using `python 2`, make sure to choose `Python 2` kernel.

2. widget not shown?

- Could not cross validate the widget frontend and backend versions (or similiar)

Double check if you are having two ipywidgets version (e.g: one installed via pip and one installed via conda)

- Class NGLModel not found in module nglview-js-widgets

You are likely using older JavaScript distribution of nglview. Check if it is 
`./.local/share/jupyter/nbextensions/nglview-js-widgets/`, if Yes, delete it.

- Debug experience from users
    - https://github.com/SBRG/ssbio/wiki/Troubleshooting#nglviewer-fresh-install-tips

3. Can I have two MDA.Atomgroups in the same view?

https://github.com/arose/nglview/issues/434

4. Can nglview handle large trajectories?

Absolutely yes. In general, trajectory data in NGLview are read (e.g. loaded from file) by external libraries. Some of the libraries, including pytraj and mdanalysis, have out-of-core readers for trajectory files, that is they don’t require loading the whole trajectory into memory for accessing and process the coordinates. With respect to the data loading aspect, this features enables viewing very large trajectory files with NGLview. The corresponding command in pytraj is `iterload` (see http://amber-md.github.io/pytraj/latest/_api/pytraj.html#pytraj.iterload). Trajectories in MDAnalysis are by default read out-of-core when using the `Universe` object to load the coordinates (see https://www.mdanalysis.org/docs/documentation_pages/core/universe.html).
