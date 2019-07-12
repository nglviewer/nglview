
module.exports = {


    "./nglview/tests/notebooks/api/test_embed.ipynb": function (browser) {
        browser.openNotebook("./nglview/tests/notebooks/api/test_embed.ipynb");
        browser.restartKernel(2000);
        for ( var i = 0; i < 6; i++) {
           browser.executeCell(i)
                  .pause(3000)
                  .cellHasError(i);
        }
        browser.end();
    },

}
