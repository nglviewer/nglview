
module.exports = {


    "nglview/tests/notebooks/api/render_image.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/render_image.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 12; i++) {
           browser.executeCell(i)
                  .pause(3000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/api/test_add_shape.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/test_add_shape.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 6; i++) {
           browser.executeCell(i)
                  .pause(3000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/api/test_add_structure_then_trajectory.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/test_add_structure_then_trajectory.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 6; i++) {
           browser.executeCell(i)
                  .pause(3000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/api/test_callbacks.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/test_callbacks.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 22; i++) {
           browser.executeCell(i)
                  .pause(3000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/api/test_component_names.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/test_component_names.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 12; i++) {
           browser.executeCell(i)
                  .pause(3000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/api/test_movie_making.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/test_movie_making.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 8; i++) {
           browser.executeCell(i)
                  .pause(3000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/api/test_no_gui_demo.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/test_no_gui_demo.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 4; i++) {
           browser.executeCell(i)
                  .pause(3000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/api/test_representation_for_small_peptide.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/test_representation_for_small_peptide.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 8; i++) {
           browser.executeCell(i)
                  .pause(3000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/api/test_representations.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/test_representations.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 20; i++) {
           browser.executeCell(i)
                  .pause(3000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/api/test_set_represetation_via_gui.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/test_set_represetation_via_gui.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 22; i++) {
           browser.executeCell(i)
                  .pause(3000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/api/test_simple_gui_by_clicking.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/test_simple_gui_by_clicking.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 12; i++) {
           browser.executeCell(i)
                  .pause(3000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/api/test_sync_n_components.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/test_sync_n_components.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 12; i++) {
           browser.executeCell(i)
                  .pause(3000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/api/test_update_representation_shortcut.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/test_update_representation_shortcut.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 14; i++) {
           browser.executeCell(i)
                  .pause(3000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/api/test_view_trajectory.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/test_view_trajectory.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 40; i++) {
           browser.executeCell(i)
                  .pause(3000)
                  .cellHasError(i);
        }
        browser.end();
    },

}
