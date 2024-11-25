
module.exports = {


    "nglview/tests/notebooks/api/test_custom_color.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/test_custom_color.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 16; i++) {
           browser.executeCell(i)
                  .pause(3000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/api/test_app_layout.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/test_app_layout.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 10; i++) {
           browser.executeCell(i)
                  .pause(3000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/api/test_initial_representations.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/test_initial_representations.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 8; i++) {
           browser.executeCell(i)
                  .pause(3000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/api/test_grid_box_3_views_with_gui.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/test_grid_box_3_views_with_gui.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 12; i++) {
           browser.executeCell(i)
                  .pause(3000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/api/test_render_image.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/test_render_image.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 12; i++) {
           browser.executeCell(i)
                  .pause(3000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/api/test_grid_view_sync_camera.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/test_grid_view_sync_camera.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 14; i++) {
           browser.executeCell(i)
                  .pause(3000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/api/test_sidebar_representation.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/test_sidebar_representation.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 4; i++) {
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


    "nglview/tests/notebooks/api/test_display_same_view_several_times.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/test_display_same_view_several_times.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 8; i++) {
           browser.executeCell(i)
                  .pause(3000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/api/test_embed.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/test_embed.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 24; i++) {
           browser.executeCell(i)
                  .pause(3000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/api/test_sync_view.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/test_sync_view.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 44; i++) {
           browser.executeCell(i)
                  .pause(3000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/api/test_add_shape.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/test_add_shape.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 14; i++) {
           browser.executeCell(i)
                  .pause(3000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/api/test_movie_making.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/test_movie_making.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 10; i++) {
           browser.executeCell(i)
                  .pause(3000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/api/test_representations.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/test_representations.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 50; i++) {
           browser.executeCell(i)
                  .pause(3000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "nglview/tests/notebooks/api/test_set_coordinates_for_structure_or_trajectory.ipynb": function (browser) {
        browser.openNotebook("nglview/tests/notebooks/api/test_set_coordinates_for_structure_or_trajectory.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 18; i++) {
           browser.executeCell(i)
                  .pause(3000)
                  .cellHasError(i);
        }
        browser.end();
    },

}
