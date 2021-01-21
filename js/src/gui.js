/**
 * @file  Gui
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */
var NGL = require('ngl');
var UI = require('./ui/ui.js').UI;
var signals = require("./lib/signals.min.js");


HTMLElement.prototype.getBoundingClientRect = (function () {
  // workaround for ie11 behavior with disconnected dom nodes

  var _getBoundingClientRect = HTMLElement.prototype.getBoundingClientRect

  return function getBoundingClientRect () {
    try {
      return _getBoundingClientRect.apply(this, arguments)
    } catch (e) {
      return {
        top: 0,
        left: 0,
        width: this.width,
        height: this.height
      }
    }
  }
}())

NGL.Widget = function () {

}

NGL.Widget.prototype = {
  constructor: NGL.Widget
}

NGL.createParameterInput = function (p, v) {
  if (!p) return

  var value = v === undefined ? p.value : v
  var input

  if (p.type === 'number') {
    input = new UI.Number(0)
      .setRange(p.min, p.max)
      .setPrecision(p.precision)
      .setValue(parseFloat(value))
  } else if (p.type === 'integer') {
    input = new UI.Integer(parseInt(value))
      .setRange(p.min, p.max)
  } else if (p.type === 'range') {
    input = new UI.Range(p.min, p.max, value, p.step)
      .setValue(parseFloat(value))
  } else if (p.type === 'boolean') {
    input = new UI.Checkbox(value)
  } else if (p.type === 'text') {
    input = new UI.Input(value)
  } else if (p.type === 'select') {
    input = new UI.Select()
      .setWidth('')
      .setOptions(p.options)
      .setValue(value)
  } else if (p.type === 'color') {
    input = new UI.ColorPopupMenu(p.label)
      .setValue(value)
  } else if (p.type === 'vector3') {
    input = new UI.Vector3(value)
      .setPrecision(p.precision)
  } else if (p.type === 'hidden') {

    // nothing to display

  } else {
    console.warn(
      'NGL.createParameterInput: unknown parameter type ' +
      "'" + p.type + "'"
    )
  }

  return input
}

/// /////////////
// Preferences

NGL.Preferences = function (id, defaultParams) {
  this.signals = {
    keyChanged: new signals.Signal()
  }

  this.id = id || 'ngl-gui'
  var dp = Object.assign({}, defaultParams)

  this.storage = {
    impostor: true,
    quality: 'auto',
    sampleLevel: 0,
    theme: 'dark',
    backgroundColor: 'black',
    overview: false,
    rotateSpeed: 2.0,
    zoomSpeed: 1.2,
    panSpeed: 0.8,
    clipNear: 0,
    clipFar: 100,
    clipDist: 10,
    fogNear: 50,
    fogFar: 100,
    cameraFov: 40,
    cameraType: 'perspective',
    lightColor: 0xdddddd,
    lightIntensity: 1.0,
    ambientColor: 0xdddddd,
    ambientIntensity: 0.2,
    hoverTimeout: 0
  }

  // overwrite default values with params
  for (var key in this.storage) {
    if (dp[ key ] !== undefined) {
      this.storage[ key ] = dp[ key ]
    }
  }
}

NGL.Preferences.prototype = {

  constructor: NGL.Preferences,

  getKey: function (key) {
    return this.storage[ key ]
  },

  setKey: function (key, value) {
    this.storage[ key ] = value
    this.signals.keyChanged.dispatch(key, value)
  },
}

// Stage

class StageWidget{
  constructor(view){
      // view: NGLView of NGLModel
      var el
      var stage = view.stage
      var viewport = new UI.Panel()

      if (view.options.parent){
          el = view.el.parentElement
      }else{
          el = view.el
      }
      viewport.setPosition("absolute")
      viewport.dom = stage.viewer.container
      this.el = el
      this.widgetList = []

      // Turn off in Jupyter notebook so user can run the next cell.
      // ensure initial focus on viewer canvas for key-stroke listening
      // stage.viewer.renderer.domElement.focus()

      var preferences = new NGL.Preferences('ngl-stage-widget', stage.getParameters())

      var pp = {}
      for (var name in preferences.storage) {
        pp[ name ] = preferences.getKey(name)
      }
      stage.setParameters(pp)

      preferences.signals.keyChanged.add(function (key, value) {
        var sp = {}
        sp[ key ] = value
        stage.setParameters(sp)
        // FIXME: remove?
        // if (key === 'theme') {
        //   setTheme(value)
        // }
      }, this)

      //

      var cssLinkElement = document.createElement('link')
      cssLinkElement.rel = 'stylesheet'
      cssLinkElement.id = 'theme'

      function setTheme (value) {
        var cssPath, bgColor
        if (value === 'light') {
          cssPath = NGL.cssDirectory + 'light.css'
          bgColor = 'white'
        } else {
          cssPath = NGL.cssDirectory + 'dark.css'
          bgColor = 'black'
        }
        cssLinkElement.href = cssPath
        stage.setParameters({ backgroundColor: bgColor })
      }

      // FIXME: remove?
      // setTheme(preferences.getKey('theme'))
      el.appendChild(cssLinkElement)

      var toolbar = new NGL.ToolbarWidget(stage).setId('toolbar_ngl')
      el.appendChild(toolbar.dom)

      var menubar = new NGL.MenubarWidget(stage, preferences).setId('menubar_ngl')
      el.appendChild(menubar.dom)

      var sidebar = new NGL.SidebarWidget(stage).setId('sidebar_ngl')
      el.appendChild(sidebar.dom)


      this.widgetList.push(toolbar)
      this.widgetList.push(menubar)
      this.widgetList.push(sidebar)

      //

      // el.body.style.touchAction = 'none'
      el.style.touchAction = 'none'

      //

      stage.handleResize()
      // FIXME hack for ie11
      setTimeout(function () { stage.handleResize() }, 500)

      //

      var doResizeLeft = false
      var movedResizeLeft = false
      var minResizeLeft = false
      var handleResizeInNotebook

      var handleResizeLeft = function (clientX) {
        if (clientX >= 50 && clientX <= window.innerWidth - 10) {
          sidebar.setWidth(window.innerWidth - clientX + 'px')
          viewport.setWidth(clientX + 'px')
          toolbar.setWidth(clientX + 'px')
          stage.handleResize()
        }
        var sidebarWidth = sidebar.dom.getBoundingClientRect().width
        if (clientX === undefined) {
          var mainWidth = window.innerWidth - sidebarWidth
          viewport.setWidth(mainWidth + 'px')
          toolbar.setWidth(mainWidth + 'px')
          stage.handleResize()
        }
        if (sidebarWidth <= 10) {
          minResizeLeft = true
        } else {
          minResizeLeft = false
        }
        handleResizeInNotebook();
      }
      handleResizeLeft = NGL.throttle(
        handleResizeLeft, 50, { leading: true, trailing: true }
      )

      var handleResizeInNotebook = function(){
          // FIXME
          var sw = sidebar.dom.getBoundingClientRect().width
          var ew = el.getBoundingClientRect().width
          if (ew == sw) {
              sw = ew / 4
          }

          var w = ew - sw + 'px'

          stage.viewer.container.style.width = w 
          stage.handleResize()
      }

      var resizeLeft = new UI.Panel()
        .setClass('ResizeLeft')
        .onMouseDown(function () {
          doResizeLeft = true
          movedResizeLeft = false
        })
        .onClick(function () {
          if (minResizeLeft) {
            handleResizeLeft(window.innerWidth - 300)
          } else if (!doResizeLeft && !movedResizeLeft) {
            handleResizeLeft(window.innerWidth - 10)
          }
        })

      sidebar.add(resizeLeft)
      this.viewport = viewport
      this.toolbar = toolbar
      this.menubar = menubar
      this.sidebar = sidebar

      handleResizeInNotebook()

      stage.signals.fullscreenChanged.add(function (isFullscreen) {
        const box = stage.viewer.container.parentElement.getBoundingClientRect()
        stage.setSize(box.width+"px", isFullscreen?box.height+'px':'300px')
        view.handleResize()
      })

      window.addEventListener(
        'mousemove', function (event) {
          if (doResizeLeft) {
            document.body.style.cursor = 'col-resize'
            movedResizeLeft = true
            handleResizeLeft(event.clientX)
          }
        }, false
      )

      window.addEventListener(
        'mouseup', function (event) {
          doResizeLeft = false
          document.body.style.cursor = ''
        }, false
      )

      window.addEventListener(
        'resize', function (event) {
          handleResizeLeft()
        }, false
      )

      //

      document.addEventListener('dragover', function (e) {
        e.stopPropagation()
        e.preventDefault()
        e.dataTransfer.dropEffect = 'none'
      }, false)

      document.addEventListener('drop', function (e) {
        e.stopPropagation()
        e.preventDefault()
      }, false)
  }

  dispose(){
       for (var i in this.widgetList){
           this.widgetList[i].dispose()
       }
  }
}

// Viewport

NGL.ViewportWidget = function (stage) {
  var viewer = stage.viewer
  var renderer = viewer.renderer

  var container = new UI.Panel()
  container.dom = viewer.container
  container.setPosition('absolute')

  var fileTypesOpen = NGL.flatten([
    NGL.ParserRegistry.getStructureExtensions(),
    NGL.ParserRegistry.getVolumeExtensions(),
    NGL.ParserRegistry.getSurfaceExtensions(),
    NGL.DecompressorRegistry.names
  ])

  // event handlers

  container.dom.addEventListener('dragover', function (e) {
    e.stopPropagation()
    e.preventDefault()
    e.dataTransfer.dropEffect = 'copy'
  }, false)

  container.dom.addEventListener('drop', function (e) {
    e.stopPropagation()
    e.preventDefault()

    var fn = function (file, callback) {
      var ext = file.name.split('.').pop().toLowerCase()
      if (NGL.ScriptExtensions.includes(ext)) {
        stage.loadScript(file).then(callback)
      } else if (fileTypesOpen.includes(ext)) {
        stage.loadFile(file, { defaultRepresentation: true }).then(callback)
      } else {
        console.error('unknown filetype: ' + ext)
        callback()
      }
    }
    var queue = new NGL.Queue(fn, e.dataTransfer.files)
  }, false)

  return container
}

// Toolbar

NGL.ToolbarWidget = function (stage) {
  var container = new UI.Panel()

  var messageText = new UI.Text()
  var messagePanel = new UI.Panel()
    .setDisplay('inline')
    .setFloat('left')
    .add(messageText)

  var statsText = new UI.Text()
  var statsPanel = new UI.Panel()
    .setDisplay('inline')
    .setFloat('right')
    .add(statsText)

  stage.signals.clicked.add(function (pickingProxy) {
    messageText.setValue(pickingProxy ? pickingProxy.getLabel() : 'nothing')
  })

  stage.viewer.stats.signals.updated.add(function () {
    if (NGL.Debug) {
      statsText.setValue(
        stage.viewer.stats.lastDuration.toFixed(2) + ' ms | ' +
        stage.viewer.stats.lastFps + ' fps'
      )
    } else {
      statsText.setValue('')
    }
  })

  container.add(messagePanel, statsPanel)

  return container
}

// Menubar

NGL.MenubarWidget = function (stage, preferences) {
  var container = new UI.Panel()

  container.add(new NGL.MenubarFileWidget(stage))
  container.add(new NGL.MenubarViewWidget(stage, preferences))
  if (NGL.examplesListUrl && NGL.examplesScriptUrl) {
    container.add(new NGL.MenubarExamplesWidget(stage))
  }
  container.add(new NGL.MenubarHelpWidget(stage, preferences))

  container.add(
    new UI.Panel().setClass('menu').setFloat('right').add(
      new UI.Text('NGL Viewer ' + NGL.Version).setClass('title')
    )
  )

  return container
}

NGL.MenubarFileWidget = function (stage) {
  var fileTypesOpen = NGL.flatten([
    NGL.ParserRegistry.getStructureExtensions(),
    NGL.ParserRegistry.getVolumeExtensions(),
    NGL.ParserRegistry.getSurfaceExtensions(),
    NGL.DecompressorRegistry.names,
    NGL.ScriptExtensions
  ])

  function fileInputOnChange (e) {
    var fn = function (file, callback) {
      var ext = file.name.split('.').pop().toLowerCase()
      if (NGL.ScriptExtensions.includes(ext)) {
        stage.loadScript(file).then(callback)
      } else if (fileTypesOpen.includes(ext)) {
        stage.loadFile(file, { defaultRepresentation: true }).then(callback)
      } else {
        console.error('unknown filetype: ' + ext)
        callback()
      }
    }
    var queue = new NGL.Queue(fn, e.target.files)
  }

  var fileInput = document.createElement('input')
  fileInput.type = 'file'
  fileInput.multiple = true
  fileInput.style.display = 'none'
  fileInput.accept = '.' + fileTypesOpen.join(',.')
  fileInput.addEventListener('change', fileInputOnChange, false)

  // export image

  var exportImageWidget = new NGL.ExportImageWidget(stage)
    .setDisplay('none')
    .attach(stage.viewer.container.parentElement)

  // event handlers

  function onOpenOptionClick () {
    fileInput.click()
  }

  function onImportOptionClick () {
    var dirWidget
    function onListingClick (info) {
      var ext = info.path.split('.').pop().toLowerCase()
      if (NGL.ScriptExtensions.includes(ext)) {
        stage.loadScript(NGL.ListingDatasource.getUrl(info.path))
        dirWidget.dispose()
      } else if (fileTypesOpen.includes(ext)) {
        stage.loadFile(NGL.ListingDatasource.getUrl(info.path), {
          defaultRepresentation: true
        })
        dirWidget.dispose()
      } else {
        console.error('unknown filetype: ' + ext)
      }
    }

    dirWidget = new NGL.DirectoryListingWidget(
      NGL.ListingDatasource, stage, 'Import file',
      fileTypesOpen, onListingClick
    )

    dirWidget
      .setOpacity('0.9')
      .setLeft('50px')
      .setTop('80px')
      .attach()
  }

  function onExportImageOptionClick () {
    exportImageWidget
      .setOpacity('0.9')
      .setLeft('50px')
      .setTop('80px')
      .setDisplay('block')
  }

  function onScreenshotOptionClick () {
    stage.makeImage({
      factor: 1,
      antialias: true,
      trim: false,
      transparent: false
    }).then(function (blob) {
      NGL.download(blob, 'screenshot.png')
    })
  }

  function onPdbInputKeyDown (e) {
    if (e.keyCode === 13) {
      stage.loadFile('rcsb://' + e.target.value.trim(), {
        defaultRepresentation: true
      })
      e.target.value = ''
    }
  }

  function onFirstModelOnlyChange (e) {
    stage.defaultFileParams.firstModelOnly = e.target.checked
  }

  function onCAlphaOnlyChange (e) {
    stage.defaultFileParams.cAlphaOnly = e.target.checked
  }

  // configure menu contents

  var createOption = UI.MenubarHelper.createOption
  var createInput = UI.MenubarHelper.createInput
  var createCheckbox = UI.MenubarHelper.createCheckbox
  var createDivider = UI.MenubarHelper.createDivider

  var menuConfig = [
    createOption('Open...', onOpenOptionClick),
    createInput('PDB', onPdbInputKeyDown),
    createCheckbox('firstModelOnly', false, onFirstModelOnlyChange),
    createCheckbox('cAlphaOnly', false, onCAlphaOnlyChange),
    createDivider(),
    createOption('Screenshot', onScreenshotOptionClick, 'camera'),
    createOption('Export image...', onExportImageOptionClick)
  ]

  if (NGL.ListingDatasource) {
    menuConfig.splice(
      1, 0, createOption('Import...', onImportOptionClick)
    )
  }

  var optionsPanel = UI.MenubarHelper.createOptionsPanel(menuConfig)
  optionsPanel.dom.appendChild(fileInput)

  return UI.MenubarHelper.createMenuContainer('File', optionsPanel)
}

NGL.MenubarViewWidget = function (stage, preferences) {
  // event handlers

  function onLightThemeOptionClick () {
    preferences.setKey('theme', 'light')
  }

  function onDarkThemeOptionClick () {
    preferences.setKey('theme', 'dark')
  }

  function onPerspectiveCameraOptionClick () {
    stage.setParameters({ cameraType: 'perspective' })
  }

  function onOrthographicCameraOptionClick () {
    stage.setParameters({ cameraType: 'orthographic' })
  }

  function onStereoCameraOptionClick () {
    stage.setParameters({ cameraType: 'stereo' })
  }

  function onFullScreenOptionClick () {
    stage.toggleFullscreen(stage.viewer.container.parentElement);
  }

  function onCenterOptionClick () {
    stage.autoView(1000)
  }

  function onToggleSpinClick () {
    stage.toggleSpin()
  }

  function onToggleRockClick () {
    stage.toggleRock()
  }

  function onGetOrientationClick () {
    window.prompt(
      'Get orientation',
      JSON.stringify(
        stage.viewerControls.getOrientation().toArray(),
        function (k, v) {
          return v.toFixed ? Number(v.toFixed(2)) : v
        }
      )
    )
  }

  function onSetOrientationClick () {
    stage.viewerControls.orient(
      JSON.parse(window.prompt('Set orientation'))
    )
  }

  var that = this;
  stage.signals.fullscreenChanged.add(function (isFullscreen) {
    // const box = stage.viewer.container.parentElement.getBoundingClientRect()
    // stage.setSize(box.width+"px", box.height+"px")
    // stage.handleResize()
    var icon = menuConfig[ 4 ].children[ 0 ]
    if (isFullscreen) {
      icon.switchClass('compress', 'expand')
    } else {
      icon.switchClass('expand', 'compress')
    }
  })

  // configure menu contents

  var createOption = UI.MenubarHelper.createOption
  var createDivider = UI.MenubarHelper.createDivider

  var menuConfig = [
    // createOption('Light theme', onLightThemeOptionClick),
    // createOption('Dark theme', onDarkThemeOptionClick),
    // createDivider(),
    createOption('Perspective', onPerspectiveCameraOptionClick),
    createOption('Orthographic', onOrthographicCameraOptionClick),
    createOption('Stereo', onStereoCameraOptionClick),
    createDivider(),
    createOption('Full screen', onFullScreenOptionClick, 'expand'),
    createOption('Center', onCenterOptionClick, 'bullseye'),
    createDivider(),
    createOption('Toggle spin', onToggleSpinClick),
    createOption('Toggle rock', onToggleRockClick),
    // createDivider(),
    // createOption('Get orientation', onGetOrientationClick),
    // createOption('Set orientation', onSetOrientationClick)
  ]

  var optionsPanel = UI.MenubarHelper.createOptionsPanel(menuConfig)

  return UI.MenubarHelper.createMenuContainer('View', optionsPanel)
}

NGL.MenubarExamplesWidget = function (stage) {
  // configure menu contents

  var createOption = UI.MenubarHelper.createOption
  var optionsPanel = UI.MenubarHelper.createOptionsPanel([])
  optionsPanel.setWidth('300px')

  var xhr = new XMLHttpRequest()
  xhr.open('GET', NGL.examplesListUrl)
  xhr.responseType = 'json'
  xhr.onload = function (e) {
    var response = this.response
    if (typeof response === 'string') {
      // for ie11
      response = JSON.parse(response)
    }
    response.sort().forEach(function (name) {
      var option = createOption(name, function () {
        stage.loadScript(NGL.examplesScriptUrl + name + '.js')
      })
      optionsPanel.add(option)
    })
  }
  xhr.send()

  return UI.MenubarHelper.createMenuContainer('Examples', optionsPanel)
}

NGL.MenubarHelpWidget = function (stage, preferences) {
  // event handlers

  function onOverviewOptionClick () {
    overviewWidget
      .setOpacity('0.9')
      .setDisplay('block')
      .setWidgetPosition(50, 80)
  }

  function onDocOptionClick () {
    window.open(NGL.documentationUrl, '_blank')
  }

  function onDebugOnClick () {
    NGL.setDebug(true)
    stage.viewer.updateHelper()
    stage.viewer.requestRender()
  }

  function onDebugOffClick () {
    NGL.setDebug(false)
    stage.viewer.updateHelper()
    stage.viewer.requestRender()
  }

  function onPreferencesOptionClick () {
    preferencesWidget
      .setOpacity('0.9')
      .setDisplay('block')
      .setWidgetPosition(50, 80)
  }

  // export image

  var preferencesWidget = new NGL.PreferencesWidget(stage, preferences)
    .setDisplay('none')
    .attach(stage.viewer.container.parentElement)

  // overview

  var overviewWidget = new NGL.OverviewWidget(stage, preferences)
    .setDisplay('none')
    .attach(stage.viewer.container.parentElement)

  if (preferences.getKey('overview')) {
    onOverviewOptionClick()
  }

  // configure menu contents

  var createOption = UI.MenubarHelper.createOption
  var createDivider = UI.MenubarHelper.createDivider

  var menuConfig = [
    createOption('Overview', onOverviewOptionClick),
    // createOption('Documentation', onDocOptionClick),
    // createDivider(),
    // createOption('Debug on', onDebugOnClick),
    // createOption('Debug off', onDebugOffClick),
    // createDivider(),
    createOption('Preferences', onPreferencesOptionClick, 'sliders')
  ]

  var optionsPanel = UI.MenubarHelper.createOptionsPanel(menuConfig)

  return UI.MenubarHelper.createMenuContainer('Help', optionsPanel)
}

// Overview

NGL.OverviewWidget = function (stage, preferences) {
  var container = new UI.OverlayPanel()

  var xOffset = 0
  var yOffset = 0

  var prevX = 0
  var prevY = 0

  function onMousemove (e) {
    if (prevX === 0) {
      prevX = e.clientX
      prevY = e.clientY
    }
    xOffset -= prevX - e.clientX
    yOffset -= prevY - e.clientY
    prevX = e.clientX
    prevY = e.clientY
    container.dom.style.top = yOffset + 'px'
    container.dom.style.left = xOffset + 'px'
  }

  function setWidgetPosition (left, top) {
    xOffset = left
    yOffset = top
    prevX = 0
    prevY = 0
    container.dom.style.top = yOffset + 'px'
    container.dom.style.left = xOffset + 'px'
  }
  container.setWidgetPosition = setWidgetPosition

  var headingPanel = new UI.Panel()
    .setBorderBottom('1px solid #555')
    .setHeight('25px')
    .setCursor('move')
    .onMouseDown(function (e) {
      if (e.which === 1) {
        document.addEventListener('mousemove', onMousemove)
      }
      document.addEventListener('mouseup', function (e) {
        document.removeEventListener('mousemove', onMousemove)
      })
    })

  var listingPanel = new UI.Panel()
    .setMarginTop('10px')
    .setMinHeight('100px')
    .setMaxHeight('500px')
    .setMaxWidth('600px')
    .setOverflow('auto')

  headingPanel.add(
    new UI.Text('NGL Viewer').setFontStyle('italic'),
    new UI.Html('&nbsp;&mdash;&nbsp;Overview')
  )
  headingPanel.add(
    new UI.Icon('times')
      .setCursor('pointer')
      .setMarginLeft('20px')
      .setFloat('right')
      .onClick(function () {
        container.setDisplay('none')
      })
  )

  container.add(headingPanel)
  container.add(listingPanel)

  //

  function addIcon (name, text) {
    var panel = new UI.Panel()

    var icon = new UI.Icon(name)
      .setWidth('20px')
      .setFloat('left')

    var label = new UI.Text(text)
      .setDisplay('inline')
      .setMarginLeft('5px')

    panel
      .setMarginLeft('20px')
      .add(icon, label)
    listingPanel.add(panel)
  }

  listingPanel
    .add(new UI.Panel().add(new UI.Html("To load a new structure use the <i>File</i> menu in the top left via drag'n'drop.")))
    .add(new UI.Break())

  listingPanel
    .add(new UI.Panel().add(new UI.Text('A number of clickable icons provide common actions. Most icons can be clicked on, just try it or hover the mouse pointer over it to see a tooltip.')))
    .add(new UI.Break())

  addIcon('eye', 'Controls the visibility of a component.')
  addIcon('trash-o', 'Deletes a component. Note that a second click is required to confirm the action.')
  addIcon('bullseye', 'Centers a component.')
  addIcon('bars', 'Opens a menu with further options.')
  addIcon('square', 'Opens a menu with coloring options.')
  addIcon('filter', 'Indicates atom-selection input fields.')

  listingPanel
    .add(new UI.Text('Mouse controls'))
    .add(new UI.Html(
      '<ul>' +
          '<li>Left button hold and move to rotate camera around center.</li>' +
          '<li>Left button click to pick atom.</li>' +
          '<li>Middle button hold and move to zoom camera in and out.</li>' +
          '<li>Middle button click to center camera on atom.</li>' +
          '<li>Right button hold and move to translate camera in the screen plane.</li>' +
      '</ul>'
    ))

  listingPanel
    .add(new UI.Panel().add(new UI.Html(
      'For more information please visit the ' +
      "<a href='" + NGL.documentationUrl + "' target='_blank'>documentation pages</a>."
    )))

  var overview = preferences.getKey('overview')
  var showOverviewCheckbox = new UI.Checkbox(overview)
    .onClick(function () {
      preferences.setKey(
        'overview',
        showOverviewCheckbox.getValue()
      )
    })

  listingPanel
    .add(new UI.HorizontalRule()
      .setBorderTop('1px solid #555')
      .setMarginTop('15px')
    )
    .add(new UI.Panel().add(
      showOverviewCheckbox,
      new UI.Text(
        'Show on startup. Always available from Menu > Help > Overview.'
      ).setMarginLeft('5px')
    ))

  return container
}

// Preferences

NGL.PreferencesWidget = function (stage, preferences) {
  var container = new UI.OverlayPanel()

  var xOffset = 0
  var yOffset = 0

  var prevX = 0
  var prevY = 0

  function onMousemove (e) {
    if (prevX === 0) {
      prevX = e.clientX
      prevY = e.clientY
    }
    xOffset -= prevX - e.clientX
    yOffset -= prevY - e.clientY
    prevX = e.clientX
    prevY = e.clientY
    container.dom.style.top = yOffset + 'px'
    container.dom.style.left = xOffset + 'px'
  }

  function setWidgetPosition (left, top) {
    xOffset = left
    yOffset = top
    prevX = 0
    prevY = 0
    container.dom.style.top = yOffset + 'px'
    container.dom.style.left = xOffset + 'px'
  }
  container.setWidgetPosition = setWidgetPosition

  var headingPanel = new UI.Panel()
    .setBorderBottom('1px solid #555')
    .setHeight('25px')
    .setCursor('move')
    .onMouseDown(function (e) {
      if (e.which === 1) {
        document.addEventListener('mousemove', onMousemove)
      }
      document.addEventListener('mouseup', function (e) {
        document.removeEventListener('mousemove', onMousemove)
      })
    })

  var listingPanel = new UI.Panel()
    .setMarginTop('10px')
    .setMinHeight('100px')
    .setMaxHeight('500px')
    .setOverflow('auto')

  headingPanel.add(new UI.Text('Preferences'))
  headingPanel.add(
    new UI.Icon('times')
      .setCursor('pointer')
      .setMarginLeft('20px')
      .setFloat('right')
      .onClick(function () {
        container.setDisplay('none')
      })
  )

  container.add(headingPanel)
  container.add(listingPanel)

  //

  Object.keys(NGL.UIStageParameters).forEach(function (name) {
    var p = NGL.UIStageParameters[ name ]
    if (p.label === undefined) p.label = name
    var input = NGL.createParameterInput(p, stage.parameters[ name ])

    if (!input) return

    preferences.signals.keyChanged.add(function (key, value) {
      if (key === name) input.setValue(value)
    })

    function setParam () {
      var sp = {}
      sp[ name ] = input.getValue()
      preferences.setKey(name, sp[ name ])
    }

    var ua = navigator.userAgent
    if (p.type === 'range' && !/Trident/.test(ua) && !/MSIE/.test(ua)) {
      input.onInput(setParam)
    } else {
      input.onChange(setParam)
    }

    listingPanel
      .add(new UI.Text(name).setWidth('120px'))
      .add(input)
      .add(new UI.Break())
  })

  return container
}

// Export image

NGL.ExportImageWidget = function (stage) {
  var container = new UI.OverlayPanel()

  var headingPanel = new UI.Panel()
    .setBorderBottom('1px solid #555')
    .setHeight('25px')

  var listingPanel = new UI.Panel()
    .setMarginTop('10px')
    .setMinHeight('100px')
    .setMaxHeight('500px')
    .setOverflow('auto')

  headingPanel.add(new UI.Text('Image export'))
  headingPanel.add(
    new UI.Icon('times')
      .setCursor('pointer')
      .setMarginLeft('20px')
      .setFloat('right')
      .onClick(function () {
        container.setDisplay('none')
      })
  )

  container.add(headingPanel)
  container.add(listingPanel)

  var factorSelect = new UI.Select()
    .setOptions({
      '1': '1x',
      '2': '2x',
      '3': '3x',
      '4': '4x',
      '5': '5x',
      '6': '6x',
      '7': '7x',
      '8': '8x',
      '9': '9x',
      '10': '10x'
    })
    .setValue('4')

  var antialiasCheckbox = new UI.Checkbox()
    .setValue(true)

  var trimCheckbox = new UI.Checkbox()
    .setValue(false)

  var transparentCheckbox = new UI.Checkbox()
    .setValue(false)

  var progress = new UI.Progress()
    .setDisplay('none')

  var exportButton = new UI.Button('export')
    .onClick(function () {
      exportButton.setDisplay('none')
      progress.setDisplay('inline-block')
      function onProgress (i, n, finished) {
        if (i === 1) {
          progress.setMax(n)
        }
        if (i >= n) {
          progress.setIndeterminate()
        } else {
          progress.setValue(i)
        }
        if (finished) {
          progress.setDisplay('none')
          exportButton.setDisplay('inline-block')
        }
      }

      setTimeout(function () {
        stage.makeImage({
          factor: parseInt(factorSelect.getValue()),
          antialias: antialiasCheckbox.getValue(),
          trim: trimCheckbox.getValue(),
          transparent: transparentCheckbox.getValue(),
          onProgress: onProgress
        }).then(function (blob) {
          NGL.download(blob, 'screenshot.png')
        })
      }, 50)
    })

  function addEntry (label, entry) {
    listingPanel
      .add(new UI.Text(label).setWidth('80px'))
      .add(entry || new UI.Panel())
      .add(new UI.Break())
  }

  addEntry('scale', factorSelect)
  addEntry('antialias', antialiasCheckbox)
  addEntry('trim', trimCheckbox)
  addEntry('transparent', transparentCheckbox)

  listingPanel.add(
    new UI.Break(),
    exportButton, progress
  )

  return container
}

// Sidebar

NGL.SidebarWidget = function (stage, view=undefined) {
  var signals = stage.signals
  var container = new UI.Panel()
  container.view = view

  var widgetContainer = new UI.Panel()
    .setClass('Content')

  var compList = []
  var widgetList = []

  function handleComponent(component) {
    var widget

    switch (component.type) {
      case 'structure':
        widget = new NGL.StructureComponentWidget(component, stage)
        break

      case 'surface':
        widget = new NGL.SurfaceComponentWidget(component, stage)
        break

      case 'volume':
        widget = new NGL.VolumeComponentWidget(component, stage)
        break

      case 'shape':
        widget = new NGL.ShapeComponentWidget(component, stage)
        break

      default:
        console.warn('NGL.SidebarWidget: component type unknown', component)
        return
    }

    widgetContainer.add(widget)

    compList.push(component)
    widgetList.push(widget)

  }

  // In case user adds components directly from notebook
  stage.compList.forEach(function(comp){
      handleComponent(comp)
  })

  signals.componentAdded.add(handleComponent)

  signals.componentRemoved.add(function (component) {
    var idx = compList.indexOf(component)

    if (idx !== -1) {
      widgetList[ idx ].dispose()

      compList.splice(idx, 1)
      widgetList.splice(idx, 1)
    }
  })

  // actions

  var expandAll = new UI.Icon('plus-square')
    .setTitle('expand all')
    .setCursor('pointer')
    .onClick(function () {
      widgetList.forEach(function (widget) {
        widget.expand()
      })
    })

  var collapseAll = new UI.Icon('minus-square')
    .setTitle('collapse all')
    .setCursor('pointer')
    .setMarginLeft('10px')
    .onClick(function () {
      widgetList.forEach(function (widget) {
        widget.collapse()
      })
    })

  var centerAll = new UI.Icon('bullseye')
    .setTitle('center all')
    .setCursor('pointer')
    .setMarginLeft('10px')
    .onClick(function () {
      stage.autoView(1000)
    })

  var disposeAll = new UI.DisposeIcon()
    .setMarginLeft('10px')
    .setDisposeFunction(function () {
      stage.removeAllComponents()
    })

  var syncRepr = new UI.Icon('superpowers')
    .setTitle('sync representation')
    .setCursor('pointer')
    .setMarginLeft('10px')
    .onClick(function(){
        if (view){
            view.syncReprWithMe()
        }
    })

  var settingsMenu = new UI.PopupMenu('cogs', 'Settings', 'window')
    .setIconTitle('settings')
    .setMarginLeft('10px')
  settingsMenu.entryLabelWidth = '120px'

  // Busy indicator

  var busy = new UI.Panel()
    .setDisplay('inline')
    .add(
      new UI.Icon('spinner')
        .addClass('spin')
        .setMarginLeft('45px')
    )

  stage.tasks.signals.countChanged.add(function (delta, count) {
    if (count > 0) {
      actions.add(busy)
    } else {
      try {
        actions.remove(busy)
      } catch (e) {
        // already removed
      }
    }
  })

  var paramNames = [
    'clipNear', 'clipFar', 'clipDist', 'fogNear', 'fogFar',
    'lightColor', 'lightIntensity', 'ambientColor', 'ambientIntensity'
  ]

  paramNames.forEach(function (name) {
    var p = NGL.UIStageParameters[ name ]
    if (p.label === undefined) p.label = name
    var input = NGL.createParameterInput(p, stage.parameters[ name ])

    if (!input) return

    stage.signals.parametersChanged.add(function (params) {
      input.setValue(params[ name ])
    })

    function setParam () {
      var sp = {}
      sp[ name ] = input.getValue()
      stage.setParameters(sp)
    }

    var ua = navigator.userAgent
    if (p.type === 'range' && !/Trident/.test(ua) && !/MSIE/.test(ua)) {
      input.onInput(setParam)
    } else {
      input.onChange(setParam)
    }

    settingsMenu.addEntry(name, input)
  })

  //

  var actions = new UI.Panel()
    .setClass('Panel Sticky')
    .add(
      expandAll,
      collapseAll,
      centerAll,
      disposeAll,
      syncRepr,
      settingsMenu
    )

  container.add(
    actions,
    widgetContainer
  )

  return container
}

// Component

NGL.StructureComponentWidget = function (component, stage) {
  var signals = component.signals
  var container = new UI.CollapsibleIconPanel('minus-square', 'plus-square')

  var reprContainer = new UI.Panel()
  var trajContainer = new UI.Panel()

  function handleRepr(repr){
    reprContainer.add(
      new NGL.RepresentationElementWidget(repr, stage, component)
    )
  }

  component.reprList.forEach(function(repr){
      handleRepr(repr)
  })

  signals.representationAdded.add(handleRepr)

  signals.trajectoryAdded.add(function (traj) {
    trajContainer.add(new NGL.TrajectoryElementWidget(traj, stage))
  })

  signals.defaultAssemblyChanged.add(function () {
    assembly.setValue(component.parameters.defaultAssembly)
  })

  // Selection

  container.add(
    new UI.SelectionPanel(component.selection)
      .setMarginLeft('20px')
      .setInputWidth('214px')
  )

  // Export PDB

  var pdb = new UI.Button('export').onClick(function () {
    var pdbWriter = new NGL.PdbWriter(component.structure)
    pdbWriter.download('structure')
    componentPanel.setMenuDisplay('none')
  })

  // Add representation

  var repr = new UI.Select()
    .setColor('#444')
    .setOptions((function () {
      var reprOptions = { '': '[ add ]' }
      NGL.RepresentationRegistry.names.forEach(function (key) {
        reprOptions[ key ] = key
      })
      return reprOptions
    })())
    .onChange(function () {
      component.addRepresentation(repr.getValue())
      repr.setValue('')
      componentPanel.setMenuDisplay('none')
    })

  // Assembly

  var assembly = new UI.Select()
    .setColor('#444')
    .setOptions((function () {
      var biomolDict = component.structure.biomolDict
      var assemblyOptions = {
        '': (component.structure.unitcell ? 'AU' : 'FULL')
      }
      Object.keys(biomolDict).forEach(function (k) {
        assemblyOptions[ k ] = k
      })
      return assemblyOptions
    })())
    .setValue(component.parameters.defaultAssembly)
    .onChange(function () {
      component.setDefaultAssembly(assembly.getValue())
      componentPanel.setMenuDisplay('none')
    })

  // Open trajectory

  var trajExt = []
  NGL.ParserRegistry.getTrajectoryExtensions().forEach(function (ext) {
    trajExt.push('.' + ext, '.' + ext + '.gz')
  })

  function framesInputOnChange (e) {
    var fn = function (file, callback) {
      NGL.autoLoad(file).then(function (frames) {
        component.addTrajectory(frames)
        callback()
      })
    }
    var queue = new NGL.Queue(fn, e.target.files)
    e.target.value = ''
  }

  var framesInput = document.createElement('input')
  framesInput.type = 'file'
  framesInput.multiple = true
  framesInput.style.display = 'none'
  framesInput.accept = trajExt.join(',')
  framesInput.addEventListener('change', framesInputOnChange, false)

  var traj = new UI.Button('open').onClick(function () {
    framesInput.click()
    componentPanel.setMenuDisplay('none')
  })

  // Import remote trajectory

  var remoteTraj = new UI.Button('import').onClick(function () {
    componentPanel.setMenuDisplay('none')

    // TODO list of extensions should be provided by trajectory datasource
    var remoteTrajExt = [
      'xtc', 'trr', 'netcdf', 'dcd', 'ncdf', 'nc', 'gro', 'pdb',
      'lammpstrj', 'xyz', 'mdcrd', 'gz', 'binpos', 'h5', 'dtr',
      'arc', 'tng', 'trj', 'trz'
    ]
    var dirWidget

    function onListingClick (info) {
      var ext = info.path.split('.').pop().toLowerCase()
      if (remoteTrajExt.indexOf(ext) !== -1) {
        component.addTrajectory(info.path + '?struc=' + component.structure.path)
        dirWidget.dispose()
      } else {
        NGL.log('unknown trajectory type: ' + ext)
      }
    }

    dirWidget = new NGL.DirectoryListingWidget(
      NGL.ListingDatasource, stage, 'Import trajectory',
      remoteTrajExt, onListingClick
    )

    dirWidget
      .setOpacity('0.9')
      .setLeft('50px')
      .setTop('80px')
      .attach()
  })

  // Superpose

  function setSuperposeOptions () {
    var superposeOptions = { '': '[ structure ]' }
    stage.eachComponent(function (o, i) {
      if (o !== component) {
        superposeOptions[ i ] = o.name
      }
    }, NGL.StructureComponent)
    superpose.setOptions(superposeOptions)
  }

  stage.signals.componentAdded.add(setSuperposeOptions)
  stage.signals.componentRemoved.add(setSuperposeOptions)

  var superpose = new UI.Select()
    .setColor('#444')
    .onChange(function () {
      component.superpose(
        stage.compList[ superpose.getValue() ],
        true
      )
      component.autoView(1000)
      superpose.setValue('')
      componentPanel.setMenuDisplay('none')
    })

  setSuperposeOptions()

  // Principal axes

  var alignAxes = new UI.Button('align').onClick(function () {
    var pa = component.structure.getPrincipalAxes()
    var q = pa.getRotationQuaternion()
    q.multiply(component.quaternion.clone().inverse())
    stage.animationControls.rotate(q)
    stage.animationControls.move(component.getCenter())
  })

  // Measurements removal

  var removeMeasurements = new UI.Button('remove').onClick(function () {
    component.removeAllMeasurements()
  })

  // Annotations visibility

  var showAnnotations = new UI.Button('show').onClick(function () {
    component.annotationList.forEach(function (annotation) {
      annotation.setVisibility(true)
    })
  })

  var hideAnnotations = new UI.Button('hide').onClick(function () {
    component.annotationList.forEach(function (annotation) {
      annotation.setVisibility(false)
    })
  })

  var annotationButtons = new UI.Panel()
    .setDisplay('inline-block')
    .add(showAnnotations, hideAnnotations)

  // Open validation

  function validationInputOnChange (e) {
    var fn = function (file, callback) {
      NGL.autoLoad(file, { ext: 'validation' }).then(function (validation) {
        component.structure.validation = validation
        callback()
      })
    }
    var queue = new NGL.Queue(fn, e.target.files)
  }

  var validationInput = document.createElement('input')
  validationInput.type = 'file'
  validationInput.style.display = 'none'
  validationInput.accept = '.xml'
  validationInput.addEventListener('change', validationInputOnChange, false)

  var vali = new UI.Button('open').onClick(function () {
    validationInput.click()
    componentPanel.setMenuDisplay('none')
  })

  // Position

  var position = new UI.Vector3()
    .onChange(function () {
      component.setPosition(position.getValue())
    })

  // Rotation

  var q = new NGL.Quaternion()
  var e = new NGL.Euler()
  var rotation = new UI.Vector3()
    .setRange(-6.28, 6.28)
    .onChange(function () {
      e.setFromVector3(rotation.getValue())
      q.setFromEuler(e)
      component.setRotation(q)
    })

  // Scale

  var scale = new UI.Number(1)
    .setRange(0.01, 100)
    .onChange(function () {
      component.setScale(scale.getValue())
    })

  // Matrix

  signals.matrixChanged.add(function () {
    position.setValue(component.position)
    rotation.setValue(e.setFromQuaternion(component.quaternion))
    scale.setValue(component.scale.x)
  })

  // Component panel

  var componentPanel = new UI.ComponentPanel(component)
    .setDisplay('inline-block')
    .setMargin('0px')
    .addMenuEntry('PDB file', pdb)
    .addMenuEntry('Representation', repr)
    .addMenuEntry('Assembly', assembly)
    .addMenuEntry('Superpose', superpose)
    .addMenuEntry(
      'File', new UI.Text(component.structure.path)
        .setMaxWidth('100px')
        .setOverflow('auto')
      // .setWordWrap( "break-word" )
    )
    .addMenuEntry('Trajectory', traj)
    .addMenuEntry('Principal axes', alignAxes)
    .addMenuEntry('Measurements', removeMeasurements)
    .addMenuEntry('Annotations', annotationButtons)
    .addMenuEntry('Validation', vali)
    .addMenuEntry('Position', position)
    .addMenuEntry('Rotation', rotation)
    .addMenuEntry('Scale', scale)

  if (NGL.ListingDatasource && NGL.TrajectoryDatasource) {
    componentPanel.addMenuEntry('Remote trajectory', remoteTraj)
  }

  // Fill container

  container
    .addStatic(componentPanel)
    .add(trajContainer)
    .add(reprContainer)

  return container
}

NGL.SurfaceComponentWidget = function (component, stage) {
  var signals = component.signals
  var container = new UI.CollapsibleIconPanel('minus-square', 'plus-square')

  var reprContainer = new UI.Panel()

  signals.representationAdded.add(function (repr) {
    reprContainer.add(
      new NGL.RepresentationElementWidget(repr, stage)
    )
  })

  // Add representation

  var repr = new UI.Select()
    .setColor('#444')
    .setOptions((function () {
      var reprOptions = {
        '': '[ add ]',
        'surface': 'surface',
        'dot': 'dot'
      }
      return reprOptions
    })())
    .onChange(function () {
      component.addRepresentation(repr.getValue())
      repr.setValue('')
      componentPanel.setMenuDisplay('none')
    })

  // Position

  var position = new UI.Vector3()
    .onChange(function () {
      component.setPosition(position.getValue())
    })

  // Rotation

  var q = new NGL.Quaternion()
  var e = new NGL.Euler()
  var rotation = new UI.Vector3()
    .setRange(-6.28, 6.28)
    .onChange(function () {
      e.setFromVector3(rotation.getValue())
      q.setFromEuler(e)
      component.setRotation(q)
    })

  // Scale

  var scale = new UI.Number(1)
    .setRange(0.01, 100)
    .onChange(function () {
      component.setScale(scale.getValue())
    })

  // Matrix

  signals.matrixChanged.add(function () {
    position.setValue(component.position)
    rotation.setValue(e.setFromQuaternion(component.quaternion))
    scale.setValue(component.scale.x)
  })

  // Component panel

  var componentPanel = new UI.ComponentPanel(component)
    .setDisplay('inline-block')
    .setMargin('0px')
    .addMenuEntry('Representation', repr)
    .addMenuEntry(
      'File', new UI.Text(component.surface.path)
        .setMaxWidth('100px')
        .setWordWrap('break-word'))
    .addMenuEntry('Position', position)
    .addMenuEntry('Rotation', rotation)
    .addMenuEntry('Scale', scale)

  // Fill container

  container
    .addStatic(componentPanel)
    .add(reprContainer)

  return container
}

NGL.VolumeComponentWidget = function (component, stage) {
  var signals = component.signals
  var container = new UI.CollapsibleIconPanel('minus-square', 'plus-square')

  var reprContainer = new UI.Panel()

  signals.representationAdded.add(function (repr) {
    reprContainer.add(
      new NGL.RepresentationElementWidget(repr, stage)
    )
  })

  // Add representation

  var repr = new UI.Select()
    .setColor('#444')
    .setOptions((function () {
      var reprOptions = {
        '': '[ add ]',
        'surface': 'surface',
        'dot': 'dot',
        'slice': 'slice'
      }
      return reprOptions
    })())
    .onChange(function () {
      component.addRepresentation(repr.getValue())
      repr.setValue('')
      componentPanel.setMenuDisplay('none')
    })

  // Position

  var position = new UI.Vector3()
    .onChange(function () {
      component.setPosition(position.getValue())
    })

  // Rotation

  var q = new NGL.Quaternion()
  var e = new NGL.Euler()
  var rotation = new UI.Vector3()
    .setRange(-6.28, 6.28)
    .onChange(function () {
      e.setFromVector3(rotation.getValue())
      q.setFromEuler(e)
      component.setRotation(q)
    })

  // Scale

  var scale = new UI.Number(1)
    .setRange(0.01, 100)
    .onChange(function () {
      component.setScale(scale.getValue())
    })

  // Matrix

  signals.matrixChanged.add(function () {
    position.setValue(component.position)
    rotation.setValue(e.setFromQuaternion(component.quaternion))
    scale.setValue(component.scale.x)
  })

  // Component panel

  var componentPanel = new UI.ComponentPanel(component)
    .setDisplay('inline-block')
    .setMargin('0px')
    .addMenuEntry('Representation', repr)
    .addMenuEntry(
      'File', new UI.Text(component.volume.path)
        .setMaxWidth('100px')
        .setWordWrap('break-word'))
    .addMenuEntry('Position', position)
    .addMenuEntry('Rotation', rotation)
    .addMenuEntry('Scale', scale)

  // Fill container

  container
    .addStatic(componentPanel)
    .add(reprContainer)

  return container
}

NGL.ShapeComponentWidget = function (component, stage) {
  var signals = component.signals
  var container = new UI.CollapsibleIconPanel('minus-square', 'plus-square')

  var reprContainer = new UI.Panel()

  signals.representationAdded.add(function (repr) {
    reprContainer.add(
      new NGL.RepresentationElementWidget(repr, stage)
    )
  })

  // Add representation

  var repr = new UI.Select()
    .setColor('#444')
    .setOptions((function () {
      var reprOptions = {
        '': '[ add ]',
        'buffer': 'buffer'
      }
      return reprOptions
    })())
    .onChange(function () {
      component.addRepresentation(repr.getValue())
      repr.setValue('')
      componentPanel.setMenuDisplay('none')
    })

  // Position

  var position = new UI.Vector3()
    .onChange(function () {
      component.setPosition(position.getValue())
    })

    // Rotation

  var q = new NGL.Quaternion()
  var e = new NGL.Euler()
  var rotation = new UI.Vector3()
    .setRange(-6.28, 6.28)
    .onChange(function () {
      e.setFromVector3(rotation.getValue())
      q.setFromEuler(e)
      component.setRotation(q)
    })

  // Scale

  var scale = new UI.Number(1)
    .setRange(0.01, 100)
    .onChange(function () {
      component.setScale(scale.getValue())
    })

  // Matrix

  signals.matrixChanged.add(function () {
    position.setValue(component.position)
    rotation.setValue(e.setFromQuaternion(component.quaternion))
    scale.setValue(component.scale.x)
  })

  // Component panel

  var componentPanel = new UI.ComponentPanel(component)
    .setDisplay('inline-block')
    .setMargin('0px')
    .addMenuEntry('Representation', repr)
    .addMenuEntry(
      'File', new UI.Text(component.shape.path)
        .setMaxWidth('100px')
        .setWordWrap('break-word'))
    .addMenuEntry('Position', position)
    .addMenuEntry('Rotation', rotation)
    .addMenuEntry('Scale', scale)

  // Fill container

  container
    .addStatic(componentPanel)
    .add(reprContainer)

  return container
}

// Representation

NGL.RepresentationElementWidget = function (element, stage, component=undefined) {
  var signals = element.signals

  var container = new UI.CollapsibleIconPanel('minus-square', 'plus-square')
    .setMarginLeft('20px')

  signals.visibilityChanged.add(function (value) {
    toggle.setValue(value)
  })

  signals.nameChanged.add(function (value) {
    name.setValue(value)
  })

  signals.disposed.add(function () {
    menu.dispose()
    container.dispose()
  })

  // Name

  var name = new UI.EllipsisText(element.name)
    .setWidth('103px')

    // Actions

  var toggle = new UI.ToggleIcon(element.visible, 'eye', 'eye-slash')
    .setTitle('hide/show')
    .setCursor('pointer')
    .setMarginLeft('25px')
    .onClick(function () {
      element.setVisibility(!element.visible)
    })

  var center = new UI.Icon('bullseye')
    .setTitle('center')
    .setCursor('pointer')
    .setMarginLeft('10px')
    .onClick(function () {
      component.autoView(element.repr.selection.string, 1000)
    })

  var disposeIcon = new UI.DisposeIcon()
    .setMarginLeft('10px')
    .setDisposeFunction(function () {
      element.dispose()
    })

  if (component){
    container
      .addStatic(name)
      .addStatic(toggle)
      .addStatic(center)
      .addStatic(disposeIcon)
  }else{
    container
      .addStatic(name)
      .addStatic(toggle)
      .addStatic(disposeIcon)
  }

  // Selection

  if ((element.parent.type === 'structure' ||
          element.parent.type === 'trajectory') &&
        element.repr.selection && element.repr.selection.type === 'selection'
  ) {
    container.add(
      new UI.SelectionPanel(element.repr.selection)
        .setMarginLeft('20px')
        .setInputWidth('194px')
    )
  }

  // Menu

  var menu = new UI.PopupMenu('bars', 'Representation')
    .setMarginLeft('45px')
    .setEntryLabelWidth('190px')

  menu.addEntry('type', new UI.Text(element.repr.type))

  // Parameters

  var repr = element.repr
  var rp = repr.getParameters()

  Object.keys(repr.parameters).forEach(function (name) {
    if (!repr.parameters[ name ]) return
    var p = Object.assign({}, repr.parameters[ name ])
    p.value = rp[ name ]
    if (p.label === undefined) p.label = name
    var input = NGL.createParameterInput(p)

    if (!input) return

    signals.parametersChanged.add(function (params) {
      if (typeof input.setValue === 'function') {
        input.setValue(params[ name ])
      }
    })

    function setParam () {
      var po = {}
      po[ name ] = input.getValue()
      element.setParameters(po)
      stage.viewer.requestRender()
    }

    var ua = navigator.userAgent
    if (p.type === 'range' && !/Trident/.test(ua) && !/MSIE/.test(ua)) {
      input.onInput(setParam)
    } else {
      input.onChange(setParam)
    }

    menu.addEntry(name, input)
  })

  container
    .addStatic(menu)

  return container
}

// Trajectory

NGL.TrajectoryElementWidget = function (element, stage) {
  var signals = element.signals
  var traj = element.trajectory

  var container = new UI.CollapsibleIconPanel('minus-square', 'plus-square')
    .setMarginLeft('20px')

  signals.disposed.add(function () {
    menu.dispose()
    container.dispose()
  })

  var frameCount = new UI.Panel()
    .setDisplay('inline')
    .add(new UI.Icon('spinner')
      .addClass('spin')
      .setMarginRight('99px')
    )

  var frameTime = new UI.Panel()
    .setMarginLeft('5px')
    .setDisplay('inline')

  function setFrame (value) {
    frame.setValue(value)
    if (traj.deltaTime && value >= 0) {
      var t = traj.getFrameTime(value) / 1000
      time.setValue(t.toFixed(9).replace(/\.?0+$/g, '') + 'ns')
    } else {
      time.setValue('')
    }
    frameRange.setValue(value)
    frameCount.clear().add(frame.setWidth('40px'))
    frameTime.clear().add(time.setWidth('90px'))
  }

  function init (value) {
    frame.setRange(-1, value - 1)
    frameRange.setRange(-1, value - 1)

    setFrame(traj.currentFrame)

    if (element.parameters.defaultStep !== undefined) {
      step.setValue(element.parameters.defaultStep)
    } else {
      // 1000 = n / step
      step.setValue(Math.ceil((value + 1) / 100))
    }

    player.setParameters({step: step.getValue()})
    player.setParameters({end: value - 1})
  }

  signals.countChanged.add(init)
  signals.frameChanged.add(setFrame)

  // Name

  var name = new UI.EllipsisText(element.parameters.name)
    .setWidth('103px')

  signals.nameChanged.add(function (value) {
    name.setValue(value)
  })

  container.addStatic(name)
  container.addStatic(frameTime)

  // frames

  var frame = new UI.Integer(-1)
    .setWidth('40px')
    .setTextAlign('right')
    .setMarginLeft('5px')
    .setRange(-1, -1)
    .onChange(function (e) {
      traj.setFrame(frame.getValue())
      menu.setMenuDisplay('none')
    })

  var time = new UI.Text()
    .setTextAlign('right')
    .setWidth('90px')

  var step = new UI.Integer(1)
    .setWidth('50px')
    .setRange(1, 10000)
    .onChange(function () {
      player.setParameters({step: step.getValue()})
    })

  var frameRow = new UI.Panel()

  var frameRange = new UI.Range(-1, -1, -1, 1)
    .setWidth('147px')
    .setMargin('0px')
    .setPadding('0px')
    .setBorder('0px')
    .onInput(function (e) {
      var value = frameRange.getValue()

      if (value === traj.currentFrame) {
        return
      }

      if (traj.player && traj.player.isRunning) {
        traj.setPlayer()
        traj.setFrame(value)
      } else if (!traj.inProgress) {
        traj.setFrame(value)
      }
    })

  var interpolateType = new UI.Select()
    .setColor('#444')
    .setOptions({
      '': 'none',
      'linear': 'linear',
      'spline': 'spline'
    })
    .setValue(element.parameters.defaultInterpolateType)
    .onChange(function () {
      player.setParameters({interpolateType: interpolateType.getValue()})
    })

  var interpolateStep = new UI.Integer(element.parameters.defaultInterpolateStep)
    .setWidth('30px')
    .setRange(1, 50)
    .onChange(function () {
      player.setParameters({interpolateStep: interpolateStep.getValue()})
    })

  var playDirection = new UI.Select()
    .setColor('#444')
    .setOptions({
      'forward': 'forward',
      'backward': 'backward',
      'bounce': 'bounce'
    })
    .setValue(element.parameters.defaultDirection)
    .onChange(function () {
      player.setParameters({direction: playDirection.getValue()})
    })

  var playMode = new UI.Select()
    .setColor('#444')
    .setOptions({
      'loop': 'loop',
      'once': 'once'
    })
    .setValue(element.parameters.defaultMode)
    .onChange(function () {
      player.setParameters({mode: playMode.getValue()})
    })

  // player

  var timeout = new UI.Integer(element.parameters.defaultTimeout)
    .setWidth('30px')
    .setRange(10, 1000)
    .onChange(function () {
      player.setParameters({timeout: timeout.getValue()})
    })

  var player = new NGL.TrajectoryPlayer(traj, {
    step: step.getValue(),
    timeout: timeout.getValue(),
    start: 0,
    end: traj.frameCount - 1,
    interpolateType: interpolateType.getValue(),
    interpolateStep: interpolateStep.getValue(),
    direction: playDirection.getValue(),
    mode: playMode.getValue()
  })
  traj.setPlayer(player)

  var playerButton = new UI.ToggleIcon(true, 'play', 'pause')
    .setMarginRight('10px')
    .setMarginLeft('20px')
    .setCursor('pointer')
    .setWidth('12px')
    .setTitle('play')
    .onClick(function () {
      player.toggle()
    })

  player.signals.startedRunning.add(function () {
    playerButton
      .setTitle('pause')
      .setValue(false)
  })

  player.signals.haltedRunning.add(function () {
    playerButton
      .setTitle('play')
      .setValue(true)
  })

  frameRow.add(playerButton, frameRange, frameCount)

  // Selection

  container.add(
    new UI.SelectionPanel(traj.selection)
      .setMarginLeft('20px')
      .setInputWidth('194px')
  )

  // Options

  var setCenterPbc = new UI.Checkbox(traj.centerPbc)
    .onChange(function () {
      element.setParameters({
        'centerPbc': setCenterPbc.getValue()
      })
    })

  var setRemovePeriodicity = new UI.Checkbox(traj.removePeriodicity)
    .onChange(function () {
      element.setParameters({
        'removePeriodicity': setRemovePeriodicity.getValue()
      })
    })

  var setRemovePbc = new UI.Checkbox(traj.removePbc)
    .onChange(function () {
      element.setParameters({
        'removePbc': setRemovePbc.getValue()
      })
    })

  var setSuperpose = new UI.Checkbox(traj.superpose)
    .onChange(function () {
      element.setParameters({
        'superpose': setSuperpose.getValue()
      })
    })

  var setDeltaTime = new UI.Number(traj.deltaTime)
    .setWidth('55px')
    .setRange(0, 1000000)
    .onChange(function () {
      element.setParameters({
        'deltaTime': setDeltaTime.getValue()
      })
    })

  var setTimeOffset = new UI.Number(traj.timeOffset)
    .setWidth('55px')
    .setRange(0, 1000000000)
    .onChange(function () {
      element.setParameters({
        'timeOffset': setTimeOffset.getValue()
      })
    })

  signals.parametersChanged.add(function (params) {
    setCenterPbc.setValue(traj.centerPbc)
    setRemovePeriodicity.setValue(traj.removePeriodicity)
    setRemovePbc.setValue(traj.removePbc)
    setSuperpose.setValue(traj.superpose)
    setDeltaTime.setValue(traj.deltaTime)
    setTimeOffset.setValue(traj.timeOffset)
    traj.setFrame(frame.getValue())
  })

  // Dispose

  var dispose = new UI.DisposeIcon()
    .setDisposeFunction(function () {
      element.parent.removeTrajectory(element)
    })

  //

  if (traj.frameCount) {
    init(traj.frameCount)
  }

  // Menu

  var menu = new UI.PopupMenu('bars', 'Trajectory')
    .setMarginLeft('10px')
    .setEntryLabelWidth('130px')
    .addEntry('Center', setCenterPbc)
    .addEntry('Remove Periodicity', setRemovePeriodicity)
    .addEntry('Remove PBC', setRemovePbc)
    .addEntry('Superpose', setSuperpose)
    .addEntry('Step size', step)
    .addEntry('Interpolation type', interpolateType)
    .addEntry('Interpolation steps', interpolateStep)
    .addEntry('Play timeout', timeout)
    .addEntry('Play direction', playDirection)
    .addEntry('Play mode', playMode)
    .addEntry('Delta time [ps]', setDeltaTime)
    .addEntry('Time offset [ps]', setTimeOffset)
    .addEntry('File',
      new UI.Text(traj.trajPath)
        .setMaxWidth('100px')
        .setWordWrap('break-word'))
    .addEntry('Dispose', dispose)

  container
    .addStatic(menu)

  container
    .add(frameRow)

  return container
}

// Listing

NGL.DirectoryListingWidget = function (datasource, stage, heading, filter, callback) {
  // from http://stackoverflow.com/a/20463021/1435042
  function fileSizeSI (a, b, c, d, e) {
    return (b = Math, c = b.log, d = 1e3, e = c(a) / c(d) | 0, a / b.pow(d, e)).toFixed(2) +
            String.fromCharCode(160) + (e ? 'kMGTPEZY'[--e] + 'B' : 'Bytes')
  }

  function getFolderDict (path) {
    path = path || ''
    var options = { '': '' }
    var full = []
    path.split('/').forEach(function (chunk) {
      full.push(chunk)
      options[ full.join('/') ] = chunk
    })
    return options
  }

  var container = new UI.OverlayPanel()

  var headingPanel = new UI.Panel()
    .setBorderBottom('1px solid #555')
    .setHeight('30px')

  var listingPanel = new UI.Panel()
    .setMarginTop('10px')
    .setMinHeight('100px')
    .setMaxHeight('500px')
    .setPaddingRight('15px')
    .setOverflow('auto')

  var folderSelect = new UI.Select()
    .setColor('#444')
    .setMarginLeft('20px')
    .setWidth('')
    .setMaxWidth('200px')
    .setOptions(getFolderDict())
    .onChange(function () {
      datasource.getListing(folderSelect.getValue())
        .then(onListingLoaded)
    })

  heading = heading || 'Directoy listing'

  headingPanel.add(new UI.Text(heading))
  headingPanel.add(folderSelect)
  headingPanel.add(
    new UI.Icon('times')
      .setCursor('pointer')
      .setMarginLeft('20px')
      .setFloat('right')
      .onClick(function () {
        container.dispose()
      })
  )

  container.add(headingPanel)
  container.add(listingPanel)

  function onListingLoaded (listing) {
    var folder = listing.path
    var data = listing.data

    NGL.lastUsedDirectory = folder
    listingPanel.clear()

    folderSelect
      .setOptions(getFolderDict(folder))
      .setValue(folder)

    data.forEach(function (info) {
      var ext = info.path.split('.').pop().toLowerCase()
      if (filter && !info.dir && filter.indexOf(ext) === -1) {
        return
      }

      var icon, name
      if (info.dir) {
        icon = 'folder-o'
        name = info.name
      } else {
        icon = 'file-o'
        name = info.name + String.fromCharCode(160) +
                '(' + fileSizeSI(info.size) + ')'
      }

      var pathRow = new UI.Panel()
        .setDisplay('block')
        .setWhiteSpace('nowrap')
        .add(new UI.Icon(icon).setWidth('20px'))
        .add(new UI.Text(name))
        .onClick(function () {
          if (info.dir) {
            datasource.getListing(info.path)
              .then(onListingLoaded)
          } else {
            callback(info)
          }
        })

      if (info.restricted) {
        pathRow.add(new UI.Icon('lock').setMarginLeft('5px'))
      }

      listingPanel.add(pathRow)
    })
  }

  datasource.getListing(NGL.lastUsedDirectory)
    .then(onListingLoaded)

  return container
}

module.exports = {
    "StageWidget": StageWidget
}
