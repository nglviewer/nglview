#!/usr/bin/env python

notebooks = [
        "nglview/tests/notebooks/test_delay.ipynb",
        "nglview/tests/notebooks/test_background.ipynb",
        "nglview/tests/notebooks/test_camera.ipynb",
        "nglview/tests/notebooks/test_load_binary_different_folder_ccp4.ipynb",
        "nglview/tests/notebooks/trajlist_pytraj.ipynb",
        "nglview/tests/notebooks/trajlist_mdtraj.ipynb",
        "nglview/tests/notebooks/trajlist_mdanalysis.ipynb",
        "nglview/tests/notebooks/trajlist_parmed.ipynb",
        "nglview/tests/notebooks/trajlist_simpletraj.ipynb",
        "nglview/tests/notebooks/remove_representations_by_name_shortcut.ipynb",
        "nglview/tests/notebooks/remove_representations_by_name.ipynb",
        "nglview/tests/notebooks/test_load_url.ipynb",
        "nglview/tests/notebooks/test_link_player.ipynb",
        "nglview/tests/notebooks/api/binary_vs_base64.ipynb",
        "nglview/tests/notebooks/duck.ipynb",
        "nglview/tests/notebooks/api/render_image.ipynb",
        "nglview/tests/notebooks/api/view_trajectory.ipynb"
        ]

head = """
module.exports = {

"""


body_template = """
    "%s": function (browser) {
        browser.openNotebook("%s");
        browser.restartKernel(2000);
        for ( var i = 0; i < %s; i++) {
           browser.executeCell(i)
                  .pause(2000)
                  .cellHasError(i);
        }
        browser.end();
    },
"""

tail = """
}
"""

if __name__ == '__main__':

    max_cells = 42
    notebook = 'nglview/tests/notebooks/test_auto_detect_pytraj_mdtraj_mdanalysis_parmed.ipynb'
    comprehensive_nb = body_template % (notebook, notebook, max_cells)

    max_cells = 20
    others  = '\n'.join(body_template % (notebook, notebook, max_cells)
                              for notebook in notebooks)

    all_notebooks = others + '\n' + comprehensive_nb
    fn = 'nglview/tests/js/test.js'
    with open(fn, 'w') as fh:
        fh.write(head + all_notebooks + tail)

