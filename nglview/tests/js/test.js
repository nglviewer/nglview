var testsuits = {};

// TODO: auto-discovery
// var nbfiles = ["repr_alias.ipynb", "init_representations.ipynb"]; 
var nbfiles = ["caching.ipynb", "representations.ipynb", 
               "closest_waters.ipynb"];

for (var i = 0; i < nbfiles.length; i++) {
    var fn = nbfiles[i];

    testsuits[fn] = function (browser) {
        browser.openNotebook(fn);
        
        browser.restartKernel(2000);
        for (var i = 0; i < 23 ; i++) {
           browser.executeCell(i)
                  .pause(1000)
                  .cellHasError(i);
        }
        
        browser.end();
    }
}

module.exports = testsuits;
