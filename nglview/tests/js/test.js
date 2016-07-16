
module.exports = {


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

}
