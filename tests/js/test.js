
module.exports = {

    "tests/notebooks/api/test_custom_color.ipynb": function (browser) {
        browser.openNotebook("tests/notebooks/api/test_custom_color.ipynb");
        browser.restartKernel(2000);
        for (var i = 0; i < 18; i++) {
            browser.executeCell(i)
                  .pause(3000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "tests/notebooks/api/test_gui_theme.ipynb": function (browser) {
        browser.openNotebook("tests/notebooks/api/test_gui_theme.ipynb");
        browser.restartKernel(2000);
        for (var i = 0; i < 22; i++) {
            browser.executeCell(i)
                  .pause(3000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "tests/notebooks/api/test_app_layout.ipynb": function (browser) {
        browser.openNotebook("tests/notebooks/api/test_app_layout.ipynb");
        browser.restartKernel(2000);
        for (var i = 0; i < 10; i++) {
            browser.executeCell(i)
                  .pause(3000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "tests/notebooks/api/test_render_image.ipynb": function (browser) {
        browser.openNotebook("tests/notebooks/api/test_render_image.ipynb");
        browser.restartKernel(2000);
        for (var i = 0; i < 12; i++) {
            browser.executeCell(i)
                  .pause(3000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "tests/notebooks/api/test_grid_view_sync_camera.ipynb": function (browser) {
        browser.openNotebook("tests/notebooks/api/test_grid_view_sync_camera.ipynb");
        browser.restartKernel(2000);
        for (var i = 0; i < 14; i++) {
            browser.executeCell(i)
                  .pause(3000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "tests/notebooks/api/test_view_trajectory.ipynb": function (browser) {
        browser.openNotebook("tests/notebooks/api/test_view_trajectory.ipynb");
        browser.restartKernel(2000);
        for (var i = 0; i < 40; i++) {
            browser.executeCell(i)
                  .pause(3000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "tests/notebooks/api/test_embed.ipynb": function (browser) {
        browser.openNotebook("tests/notebooks/api/test_embed.ipynb");
        browser.restartKernel(2000);
        for (var i = 0; i < 16; i++) {
            browser.executeCell(i)
                  .pause(3000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "tests/notebooks/api/test_sync_view.ipynb": function (browser) {
        browser.openNotebook("tests/notebooks/api/test_sync_view.ipynb");
        browser.restartKernel(2000);
        for (var i = 0; i < 44; i++) {
            browser.executeCell(i)
                  .pause(3000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "tests/notebooks/api/test_add_shape.ipynb": function (browser) {
        browser.openNotebook("tests/notebooks/api/test_add_shape.ipynb");
        browser.restartKernel(2000);
        for (var i = 0; i < 14; i++) {
            browser.executeCell(i)
                  .pause(3000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "tests/notebooks/api/test_set_coordinates.ipynb": function (browser) {
        browser.openNotebook("tests/notebooks/api/test_set_coordinates.ipynb");
        browser.restartKernel(2000);
        for (var i = 0; i < 18; i++) {
            browser.executeCell(i)
                  .pause(3000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "tests/notebooks/api/test_movie_making.ipynb": function (browser) {
        browser.openNotebook("tests/notebooks/api/test_movie_making.ipynb");
        browser.restartKernel(2000);
        for (var i = 0; i < 28; i++) {
            browser.executeCell(i)
                  .pause(3000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "tests/notebooks/api/test_representations.ipynb": function (browser) {
        browser.openNotebook("tests/notebooks/api/test_representations.ipynb");
        browser.restartKernel(2000);
        for (var i = 0; i < 50; i++) {
            browser.executeCell(i)
                  .pause(3000)
                  .cellHasError(i);
        }
        browser.end();
    },

}
