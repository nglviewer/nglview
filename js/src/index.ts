module.exports = {};

var loadedModules = [
    require("./widget_ngl.js"),
]

for (var i in loadedModules) {
    if (loadedModules && Object.prototype.hasOwnProperty.call(loadedModules, i)) {
        var loadedModule = loadedModules[i];
        for (var target_name in loadedModule) {
            if (loadedModule && Object.prototype.hasOwnProperty.call(loadedModule, target_name)) {
                module.exports[target_name] = loadedModule[target_name];
            }
        }
    }
}

module.exports['version'] = require('../package.json').version;
