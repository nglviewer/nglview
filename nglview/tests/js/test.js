
module.exports = {


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


    "nglview/tests/notebooks/api/test_add_structure_then_trajectory.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/test_add_structure_then_trajectory.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 8; i++) {
           browser.executeCell(i)
                  .pause(3000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/api/test_automatically_added_attributes_0.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/test_automatically_added_attributes_0.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 32; i++) {
           browser.executeCell(i)
                  .pause(3000)
                  .cellHasError(i);
        }
        browser.end();
    },

}
