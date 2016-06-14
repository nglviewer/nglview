Changelog
=========

Version 0.5
-----------

* ADD: `delay`
* ADD: `show_url`
* ADD: `hide` and `show_only` methods
* ADD: `show_rdkit`
* ADD: `clear` for clearing representations
* ADD: `add_unitcell` representation
* ADD: attribute `view.trajectory_0`, `view.trajectory_1`, ...
* ADD: drag and drop file to widget
* ADD: `player`
* ADD: `camera` (`perspective` and `orthographic`)
* ADD: `background`
* ADD: `orientation` for syncing two widget cameras
* ADD: viewing many trajectories, `add_trajectory`, `add_component`
* ADD: `download_image`
* ADD: `render_image`
* ADD: `center_view`
* ADD: atom selection by array
* ADD: shortcut for add_representation (add_cartoon, add_rope, ...)
* ADD: `ParmEdTrajectory` adaptor and `show_parmed`
* ENH: `add_trajectory` can accept `pytraj.Trajectory`, `mdtraj.Trajectory`, `MDAnalysis.Universe`, `parmed.Structure`
* ENH: add_representations can accept both snake and CamelKeyword
* ENH: smoother trajectory playing
* ENH: smoother rendering if adding new representation
* MIGRATION: change html/static to js/
* MIGRATION: change sending base64 to binary coordinates
* MIGRATION: change `view.trajectory` to `view._trajlist`
* MIGRATION: change `get_frame_count` method to `n_frames` property
* MIGRATION: remove `get_coordinates_dict`
* FIX: symlink error

Version 0.4
-----------

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.46373.svg)](http://dx.doi.org/10.5281/zenodo.46373)

* ADD: Convenience methods to show widget from various sources
* ADD: `PyTrajTrajectory` adaptor
* ADD: `MDAnalysisTrajectory` adaptor
* ADD: `NGLWidget.add_representation()` method
* ADD: append a "WebGL not supported message" to widget if so
* ADD: `parameters` widget property, passed to NGL stage object
* ADD: `params` property for `Structure`, dict passed to NGL
* CODE: be less noisy when importing nglview
* DOC: more usage examples, API description
* DOC: added CHANGELOG file
* BUILD: added example files in the package


Version 0.3
-----------

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.44700.svg)](http://dx.doi.org/10.5281/zenodo.44700)

* MIGRATION: `Trajectory` classes need `get_frame_count` method
* MIGRATION: removed `set_frame` method use new `frame` property
* ADD: simple trajectory player
* ADD: widget resizing support
* ADD: picking support (gui info; `picked` property)
* CODE: check for file existence in `FileStructure` and `SimpletrajTrajectory`


Version 0.2
-----------

* MIGRATION: changed `get_string` to `get_structure_string`
* MIGRATION: changed `get_coordinates` to `get_coordinates_list`
* DOC: usage, interface classes
* ADD: MDTrajTrajectory adaptor
* CODE: added interface classes
* CODE: suggested packages; mdtraj, simpletraj


Version 0.1
-----------

* initial version, no release
