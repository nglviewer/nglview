module.exports = {
    "Test repr_alias": function (browser) {
        browser.openNotebook("repr_alias.ipynb");
        
        browser.restartKernel(2000);
        for ( var i = 0; i < 23 ; i++) {
           browser.executeCell(i)
                  .pause(1000)
                  .cellHasError(i);
        }
        
        browser.end();
    },

    "Test initialize representation": function (browser) {
        browser.openNotebook("init_representations.ipynb");
        
        browser.restartKernel(2000);
        for ( var i = 0; i < 23 ; i++) {
           browser.executeCell(i)
                  .pause(1000)
                  .cellHasError(i);
        }
        
        browser.end();
    }
}
