module.exports = {
    "view_trajectory": function (browser) {
        browser.openNotebook("api/view_trajectory.ipynb");
        
        browser.restartKernel(2000);
        for ( var i = 0; i < 20; i++) {
           browser.executeCell(i)
                  .pause(5000)
                  .cellHasError(i);
        }
        browser.end();
    }
    
}
