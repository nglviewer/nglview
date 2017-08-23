### Can not import nglview although successfully installed it?

> You can try

> ```bash
> python -m ipykernel install --user
> ```

> Then in your Jupyter notebook, choose the right `kernel`. If you are using `python 2`, make sure to choose `Python 2` kernel.

### widget not shown?

> try

> ```bash
> conda install traitlets=4.2.1 ipywidgets==4.1.1 notebook=4.1.0
> ```

> - try `ipywidgets >= 5.0` and want to go back to `ipywidgets 4.1.1`?

> supposed you created `ipywidgets5` conda env to test, then if you want to go back to root env
> ```bash
> source deactivate ipywidgets5

> # need to re-enable ipykernel
> python -m ipykernel install --user
> ```

### Can I have two MDA.Atomgroups in the same view?

https://github.com/arose/nglview/issues/434
