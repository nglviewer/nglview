
module.exports = {


    "duck.ipynb": function (browser) {
        browser.openNotebook("duck.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 20; i++) {
           browser.executeCell(i)
                  .pause(2000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "api/render_image.ipynb": function (browser) {
        browser.openNotebook("api/render_image.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 20; i++) {
           browser.executeCell(i)
                  .pause(2000)
                  .cellHasError(i);
        }
        browser.end();
    },


    "api/view_trajectory.ipynb": function (browser) {
        browser.openNotebook("api/view_trajectory.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 20; i++) {
           browser.executeCell(i)
                  .pause(2000)
                  .cellHasError(i);
        }
        browser.end();
    },

}
