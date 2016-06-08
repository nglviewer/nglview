#!/usr/bin/env python

notebooks = [
        "remove_representations_by_name_shortcut.ipynb",
        "remove_representations_by_name.ipynb",
        "test_load_url.ipynb",
        "test_link_player.ipynb",
        "api/binary_vs_base64.ipynb",
        "duck.ipynb",
        "api/render_image.ipynb",
        "api/view_trajectory.ipynb"
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
    notebook = 'test_auto_detect_pytraj_mdtraj_mdanalysis_parmed.ipynb'
    comprehensive_nb = body_template % (notebook, notebook, max_cells)

    max_cells = 20
    others  = '\n'.join(body_template % (notebook, notebook, max_cells)
                              for notebook in notebooks)

    all_notebooks = comprehensive_nb + '\n' + others
    fn = 'nglview/tests/js/test.js'
    with open(fn, 'w') as fh:
        fh.write(head + all_notebooks + tail)

