
module.exports = {


    "nglview/tests/notebooks/trajlist_pytraj.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/trajlist_pytraj.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 20; i++) {
           browser.executeCell(i)
                  .pause(2000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/trajlist_mdtraj.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/trajlist_mdtraj.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 20; i++) {
           browser.executeCell(i)
                  .pause(2000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/trajlist_mdanalysis.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/trajlist_mdanalysis.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 20; i++) {
           browser.executeCell(i)
                  .pause(2000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/trajlist_parmed.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/trajlist_parmed.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 20; i++) {
           browser.executeCell(i)
                  .pause(2000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/trajlist_simpletraj.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/trajlist_simpletraj.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 20; i++) {
           browser.executeCell(i)
                  .pause(2000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/remove_representations_by_name_shortcut.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/remove_representations_by_name_shortcut.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 20; i++) {
           browser.executeCell(i)
                  .pause(2000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/remove_representations_by_name.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/remove_representations_by_name.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 20; i++) {
           browser.executeCell(i)
                  .pause(2000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/test_load_url.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/test_load_url.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 20; i++) {
           browser.executeCell(i)
                  .pause(2000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/test_link_player.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/test_link_player.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 20; i++) {
           browser.executeCell(i)
                  .pause(2000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/api/binary_vs_base64.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/binary_vs_base64.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 20; i++) {
           browser.executeCell(i)
                  .pause(2000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/duck.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/duck.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 20; i++) {
           browser.executeCell(i)
                  .pause(2000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/api/render_image.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/render_image.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 20; i++) {
           browser.executeCell(i)
                  .pause(2000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/api/view_trajectory.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/view_trajectory.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 20; i++) {
           browser.executeCell(i)
                  .pause(2000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/test_auto_detect_pytraj_mdtraj_mdanalysis_parmed.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/test_auto_detect_pytraj_mdtraj_mdanalysis_parmed.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 42; i++) {
           browser.executeCell(i)
                  .pause(2000)
                  .cellHasError(i);
        }
        browser.end();
    },

}
