Changelog
=========

Dev
---

-  MIGRATION: only support ``ipywidgets >= 5.2``
-  ADD: ``ngl v0.8``

Version 0.5.2
-------------

-  ADD: ``nglview`` command line for quickly opening notebook.
-  ADD: ``add_text`` to display
-  FIX: view.background assignment.

Version 0.5.1
-------------

-  FIX: correctly handle player when adding a structure first
-  FIX: cleaned-up outdated notebooks

Version 0.5
-----------

|DOI|

-  ADD: basic use of binary messages for faster coordinate data transfer
-  ADD: updated to newer ``NGL`` version, (e.g. orthographic camera,
   anti-aliasing, much faster rendering, ...)
-  ADD: ``view.delay`` property to adjust trajectory playback
-  ADD: ``show_url`` function for loading & viewing from an url
-  ADD: ``view.hide()`` and ``view.show_only()`` convenience methods
-  ADD: ``show_rdkit`` function to view rdkit Mol objects
-  ADD: ``view.clear()`` for clearing all representations
-  ADD: ``view.add_unitcell()``, support for the new NGL unitcell
   representation
-  ADD: automatic creation of attributes: ``view.trajectory_0``,
   ``view.trajectory_1``, ``view.component_0``, ...
-  ADD: specify ``component`` index in ``view.add_representations``,
   ``view.clear_representations``
-  ADD: drag and drop files into widget
-  ADD: extra TrajectoryPlayer widget available via ``view.player``
-  ADD: ``view.camera`` property (``perspective`` and ``orthographic``,
   new in NGL)
-  ADD: ``view.background`` property to change the viewer background
   (color name or hex value)
-  ADD: ``view.orientation`` property to synchronize cameras of two
   viewer widgets
-  ADD: view multiple trajectories via ``view.add_trajectory()``,
   ``view.add_component()``
-  ADD: ``view.download_image()`` to render high quality image and
   trigger download in browser
-  ADD: ``view.render_image()`` to render image for display in the
   notebook
-  ADD: ``view.center_view()`` to center
-  ADD: select atoms by their indices, given as a List
-  ADD: shortcuts for view.add\_representation (view.add\_cartoon,
   view.add\_rope, ...)
-  ADD: ``ParmEdTrajectory`` adaptor and ``show_parmed`` function to
   view ParmEd objects
-  ENH: ``view.add_trajectory`` now accept ``pytraj.Trajectory``,
   ``mdtraj.Trajectory``, ``MDAnalysis.Universe``, ``parmed.Structure``
-  ENH: view.add\_representations now accepts keywords in both
   snake\_case and CamelCase
-  ENH: smoother trajectory playing
-  ENH: smoother rendering when adding new representation
-  MIGRATION: change html/static to js/
-  MIGRATION: change sending base64 to binary coordinates
-  MIGRATION: change ``view.trajectory`` to ``view.trajlist``
-  MIGRATION: change ``get_frame_count`` method to ``n_frames`` property
-  MIGRATION: remove ``get_coordinates_dict``
-  FIX: symlink error

Version 0.4
-----------

|DOI|

-  ADD: Convenience methods to show widget from various sources
-  ADD: ``PyTrajTrajectory`` adaptor
-  ADD: ``MDAnalysisTrajectory`` adaptor
-  ADD: ``NGLWidget.add_representation()`` method
-  ADD: append a "WebGL not supported message" to widget if so
-  ADD: ``parameters`` widget property, passed to NGL stage object
-  ADD: ``params`` property for ``Structure``, dict passed to NGL
-  CODE: be less noisy when importing nglview
-  DOC: more usage examples, API description
-  DOC: added CHANGELOG file
-  BUILD: added example files in the package

Version 0.3
-----------

|DOI|

-  MIGRATION: ``Trajectory`` classes need ``get_frame_count`` method
-  MIGRATION: removed ``set_frame`` method use new ``frame`` property
-  ADD: simple trajectory player
-  ADD: widget resizing support
-  ADD: picking support (gui info; ``picked`` property)
-  CODE: check for file existence in ``FileStructure`` and
   ``SimpletrajTrajectory``

Version 0.2
-----------

-  MIGRATION: changed ``get_string`` to ``get_structure_string``
-  MIGRATION: changed ``get_coordinates`` to ``get_coordinates_list``
-  DOC: usage, interface classes
-  ADD: MDTrajTrajectory adaptor
-  CODE: added interface classes
-  CODE: suggested packages; mdtraj, simpletraj

Version 0.1
-----------

-  initial version, no release

.. |DOI| image:: https://zenodo.org/badge/doi/10.5281/zenodo.55409.svg
   :target: http://dx.doi.org/10.5281/zenodo.55409
.. |DOI| image:: https://zenodo.org/badge/doi/10.5281/zenodo.46373.svg
   :target: http://dx.doi.org/10.5281/zenodo.46373
.. |DOI| image:: https://zenodo.org/badge/doi/10.5281/zenodo.44700.svg
   :target: http://dx.doi.org/10.5281/zenodo.44700
