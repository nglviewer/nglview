- [Can not import nglview although successfully installed it?](#can-not-import-nglview-although-successfully-installed-it)
- [widget not shown?](#widget-not-shown)
- [Can nglview handle large trajectories?](#can-nglview-handle-large-trajectories)
- [Default views using HBox are not centered](#default-views-using-hbox-are-not-centered)
- [How to make nglview view object write PNG file?](#how-to-make-nglview-view-object-write-png-file)
- [How to make embedded nglview widget into specific size?](#how-to-make-embedded-nglview-widget-into-specific-size)
- [Is it possible to make a shape assume a different position in each frame?](#is-it-possible-to-make-a-shape-assume-a-different-position-in-each-frame)

# Can not import nglview although successfully installed it?

You can try

```bash
python -m ipykernel install
```

Then in your Jupyter notebook, choose the right `kernel`. If you are using `python 2`, make sure to choose `Python 2` kernel.

# widget not shown?
- If you're using the latest nglview, you can try below first
```
jupyter-nbextension enable nglview --py --sys-prefix
```

- Could not cross validate the widget frontend and backend versions (or similiar)

Double check if you are having two ipywidgets version (e.g: one installed via pip and one installed via conda)

- Class NGLModel not found in module nglview-js-widgets (Message can be observed in the web developer console view of your favorite browser)

You are likely using older JavaScript distribution of nglview. Check if it is 
`$HOME/.local/share/jupyter/nbextensions/nglview-js-widgets/`, if Yes, delete it.

If you're using macos, might need to delete the folder `$HOME/Library/Jupyter/nbextensions/nglview-js-widgets/`

Why? This directory has a higher preference over sys-prefix so notebook will load Javascripts files from here first.

- Extensive debug experience from users
    - https://github.com/SBRG/ssbio/wiki/Troubleshooting#nglviewer-fresh-install-tips

# Can nglview render image from command line without the notebook?

No, nglview can not do that. nglview (using NGL) needs browser for rendering.

# Can nglview handle large trajectories?

Absolutely yes. In general, trajectory data in NGLview are read (e.g. loaded from file) by external libraries. Some of the libraries, including pytraj and mdanalysis, have out-of-core readers for trajectory files, that is they don’t require loading the whole trajectory into memory for accessing and process the coordinates. With respect to the data loading aspect, this features enables viewing very large trajectory files with NGLview. The corresponding command in pytraj is `iterload` (see http://amber-md.github.io/pytraj/latest/_api/pytraj.html#pytraj.iterload). Trajectories in MDAnalysis are by default read out-of-core when using the `Universe` object to load the coordinates (see https://www.mdanalysis.org/docs/documentation_pages/core/universe.html).

# Default views using HBox are not centered

Please use `view.layout.width = 'auto'`

# How to make nglview view object write PNG file?
In other words: What is a correct way to get image data to your script (not trigger image download to your browser)?

You have to run your script in another Thread to wait for the image data (via time.sleep). Using different thread to avoid blocking the notebook code (it is important to put the code in two separate cells):
```python
import time
import nglview as nv

view = nv.demo()
view
```
```python
def generate_images(v=view):
    v.clear()
    v.add_cartoon(color='red')
    im0 = v.render_image()
    v.clear()
    v.add_cartoon(color='blue')
    im1 = v.render_image()
    for im in [im0, im1]:
        while not im.value:
            time.sleep(0.1)
    for n, im in zip('ab', [im0, im1]):
        with open(f'figure_{n}.png', 'wb') as fh:
            fh.write(im.value)

import threading
thread = threading.Thread(
    target=generate_images,
)
thread.daemon = True
thread.start()
```

# How to make embedded nglview widget into specific size?

You need to wrap it in the `ipywidgets.Box` object of appropriate size (in two separate notebook cells);

```python
import nglview as nv
from ipywidgets import Box

v = nv.demo()
box = Box([v])
box.layout.width = '600px'
box.layout.height = '600px'
box
```
```python
nv.write_html('embed.html', box)
```

#  Is it possible to make a shape assume a different position in each frame?
https://github.com/nglviewer/nglview/discussions/1002#discussioncomment-2118844
