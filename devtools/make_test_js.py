#!/usr/bin/env python

notebooks = [
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
        for ( var i = 0; i < 20; i++) {
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
    all_notebooks = '\n'.join(body_template % (notebook, notebook)
                              for notebook in notebooks)
    fn = 'nglview/tests/js/test.js'
    with open(fn, 'w') as fh:
        fh.write(head + all_notebooks + tail)
