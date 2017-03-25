define(function() { return /******/ (function(modules) { // webpackBootstrap
/******/ 	// The module cache
/******/ 	var installedModules = {};
/******/
/******/ 	// The require function
/******/ 	function __webpack_require__(moduleId) {
/******/
/******/ 		// Check if module is in cache
/******/ 		if(installedModules[moduleId])
/******/ 			return installedModules[moduleId].exports;
/******/
/******/ 		// Create a new module (and put it into the cache)
/******/ 		var module = installedModules[moduleId] = {
/******/ 			exports: {},
/******/ 			id: moduleId,
/******/ 			loaded: false
/******/ 		};
/******/
/******/ 		// Execute the module function
/******/ 		modules[moduleId].call(module.exports, module, module.exports, __webpack_require__);
/******/
/******/ 		// Flag the module as loaded
/******/ 		module.loaded = true;
/******/
/******/ 		// Return the exports of the module
/******/ 		return module.exports;
/******/ 	}
/******/
/******/
/******/ 	// expose the modules object (__webpack_modules__)
/******/ 	__webpack_require__.m = modules;
/******/
/******/ 	// expose the module cache
/******/ 	__webpack_require__.c = installedModules;
/******/
/******/ 	// __webpack_public_path__
/******/ 	__webpack_require__.p = "";
/******/
/******/ 	// Load entry module and return exports
/******/ 	return __webpack_require__(0);
/******/ })
/************************************************************************/
/******/ ([
/* 0 */
/***/ function(module, exports, __webpack_require__) {

	module.exports = {};
	
	var loadedModules = [
	    __webpack_require__(!(function webpackMissingModule() { var e = new Error("Cannot find module \"./widget_ngl.js\""); e.code = 'MODULE_NOT_FOUND'; throw e; }())),
	    //require("./ngl.js"),
	]
	
	for (var i in loadedModules) {
	    if (loadedModules.hasOwnProperty(i)) {
	        var loadedModule = loadedModules[i];
	        for (var target_name in loadedModule) {
	            if (loadedModule.hasOwnProperty(target_name)) {
	                module.exports[target_name] = loadedModule[target_name];
	            }
	        }
	    }
	}
	
	module.exports['version'] = __webpack_require__(2).version;


/***/ },
/* 1 */,
/* 2 */
/***/ function(module, exports) {

	module.exports = {
		"name": "nglview-js-widgets",
		"version": "0.5.4-dev.1",
		"description": "nglview-js-widgets",
		"author": "Hai Nguyen <hainm.comp@gmail.com>, Alexander Rose <alexander.rose@weirdbyte.de>",
		"license": "MIT",
		"main": "src/index.js",
		"repository": {
			"type": "git",
			"url": "https://github.com/arose/nglview.git"
		},
		"bugs": {
			"url": "https://github.com/arose/nglview/issues"
		},
		"files": [
			"dist"
		],
		"keywords": [
			"molecular graphics",
			"molecular structure",
			"jupyter",
			"widgets",
			"ipython",
			"ipywidgets",
			"science"
		],
		"scripts": {
			"lint": "eslint src test",
			"prepublish": "webpack",
			"test": "mocha"
		},
		"devDependencies": {
			"@jupyterlab/extension-builder": "^0.8.1",
			"babel-eslint": "^7.0.0",
			"babel-register": "^6.11.6",
			"eslint": "^3.2.2",
			"eslint-config-google": "^0.7.1",
			"json-loader": "^0.5.4",
			"ngl": "0.10.0-dev.8",
			"webpack": "^1.12.14"
		},
		"dependencies": {
			"jquery": "^2.1.4",
			"jquery-ui": "^1.12.1",
			"jupyter-js-widgets": "^2.0.6",
			"@jupyterlab/nbwidgets": "^0.6.0",
			"underscore": "^1.8.3"
		}
	};

/***/ }
/******/ ])});;
//# sourceMappingURL=index.js.map