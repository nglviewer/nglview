
module.exports = {


    "nglview/tests/notebooks/api/test_auto_detect_pytraj_mdtraj_mdanalysis_parmed.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/test_auto_detect_pytraj_mdtraj_mdanalysis_parmed.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 82; i++) {
           browser.executeCell(i)
                  .pause(2000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/api/test_view_trajectory.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/test_view_trajectory.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 38; i++) {
           browser.executeCell(i)
                  .pause(2000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/api/test_no_gui_demo.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/test_no_gui_demo.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 4; i++) {
           browser.executeCell(i)
                  .pause(2000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/api/test_add_structure_then_trajectory.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/test_add_structure_then_trajectory.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 8; i++) {
           browser.executeCell(i)
                  .pause(2000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/api/test_automatically_added_attributes_0.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/test_automatically_added_attributes_0.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 42; i++) {
           browser.executeCell(i)
                  .pause(2000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/api/test_fix_player_if_adding_single_struture_first.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/test_fix_player_if_adding_single_struture_first.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 24; i++) {
           browser.executeCell(i)
                  .pause(2000)
                  .cellHasError(i);
        }
        browser.end();
    },

}
