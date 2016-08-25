
module.exports = {


    "nglview/tests/notebooks/api/binary_vs_base64.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/binary_vs_base64.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 34; i++) {
           browser.executeCell(i)
                  .pause(3000)
                  .cellHasError(i);
        }
        browser.end();
    },


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


    "nglview/tests/notebooks/api/test_2_widget_views.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/test_2_widget_views.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 8; i++) {
           browser.executeCell(i)
                  .pause(3000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/api/test_2_widget_views_fist_time_loaded.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/test_2_widget_views_fist_time_loaded.ipynb");
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
        for ( var i = 0; i < 10; i++) {
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


    "nglview/tests/notebooks/api/test_auto_detect_pytraj_mdtraj_mdanalysis_parmed.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/test_auto_detect_pytraj_mdtraj_mdanalysis_parmed.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 82; i++) {
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


    "nglview/tests/notebooks/api/test_component_dropdown_options.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/test_component_dropdown_options.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 12; i++) {
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


    "nglview/tests/notebooks/api/test_detach.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/test_detach.ipynb");
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


    "nglview/tests/notebooks/api/test_removing_all_comopnents_and_clear_all_info.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/test_removing_all_comopnents_and_clear_all_info.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 12; i++) {
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
