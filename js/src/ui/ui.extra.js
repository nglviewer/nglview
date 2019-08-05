/**
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// Html

var UI = require("./ui.js").UI;

UI.Html = function (html) {
  UI.Element.call(this)

  var dom = document.createElement('span')
  dom.className = 'Html'
  dom.style.cursor = 'default'
  dom.style.display = 'inline-block'
  dom.style.verticalAlign = 'middle'

  this.dom = dom
  this.setValue(html)

  return this
}

UI.Html.prototype = Object.create(UI.Element.prototype)

UI.Html.prototype.setValue = function (value) {
  if (value !== undefined) {
    this.dom.innerHTML = value
  }

  return this
}

// Form

UI.Form = function () {
  UI.Panel.call(this)

  var dom = document.createElement('form')
  dom.className = 'Form'
  dom.method = 'post'
  dom.action = ''
  dom.target = '_blank'
  dom.enctype = 'multipart/form-data'

  this.dom = dom

  return this
}

UI.Form.prototype = Object.create(UI.Panel.prototype)

UI.Form.prototype.setMethod = function (value) {
  this.dom.method = value

  return this
}

UI.Form.prototype.setAction = function (value) {
  this.dom.action = value

  return this
}

UI.Form.prototype.setTarget = function (value) {
  this.dom.target = value

  return this
}

UI.Form.prototype.setEnctype = function (value) {
  this.dom.enctype = value

  return this
}

// File

UI.File = function () {
  UI.Input.call(this)

  this.dom.className = 'File'
  this.dom.type = 'file'
  this.dom.multiple = false

  return this
}

UI.File.prototype = Object.create(UI.Input.prototype)

UI.File.prototype.setMultiple = function (value) {
  this.dom.multiple = value

  return this
}

UI.File.prototype.getFiles = function (value) {
  return this.dom.files
}

// Hidden

UI.Hidden = function () {
  UI.Input.call(this)

  this.dom.type = 'hidden'

  return this
}

UI.Hidden.prototype = Object.create(UI.Input.prototype)

// Toggle Button

UI.ToggleButton = function (labelA, labelB, callbackA, callbackB) {
  UI.Button.call(this, labelB)

  var flag = true

  this.onClick(function () {
    if (flag) {
      flag = false
      this.setLabel(labelA)
      callbackB()
    } else {
      flag = true
      this.setLabel(labelB)
      callbackA()
    }
  }.bind(this))

  return this
}

UI.ToggleButton.prototype = Object.create(UI.Button.prototype)

// Ellipsis Text

UI.EllipsisText = function (text) {
  UI.Text.call(this, text)

  this.setWhiteSpace('nowrap')
  this.setOverflow('hidden')
  this.setTextOverflow('ellipsis')

  return this
}

UI.EllipsisText.prototype = Object.create(UI.Text.prototype)

UI.EllipsisText.prototype.setValue = function (value) {
  if (value !== undefined) {
    this.dom.textContent = value
    this.setTitle(value)
  }

  return this
}

// Ellipsis Multiline Text

UI.EllipsisMultilineText = function (text) {
  // http://www.mobify.com/blog/multiline-ellipsis-in-pure-css/

  UI.Element.call(this)

  var dom = document.createElement('span')
  dom.className = 'EllipsisMultilineText'
  dom.style.cursor = 'default'
  dom.style.display = 'inline-block'
  dom.style.verticalAlign = 'middle'

  var content = document.createElement('p')
  dom.appendChild(content)

  this.dom = dom
  this.content = content

  this.setValue(text)

  return this
}

UI.EllipsisMultilineText.prototype = Object.create(UI.Element.prototype)

UI.EllipsisMultilineText.prototype.setValue = function (value) {
  if (value !== undefined) {
    this.content.textContent = value
    this.setTitle(value)
  }

  return this
}

// Overlay Panel

UI.OverlayPanel = function () {
  UI.Panel.call(this)

  this.dom.className = 'Panel OverlayPanel'
  this.dom.tabIndex = '-1'
  this.dom.style.outline = 'none'

  return this
}

UI.OverlayPanel.prototype = Object.create(UI.Panel.prototype)

UI.OverlayPanel.prototype.attach = function (node) {
  node = node || document.body

  node.appendChild(this.dom)

  return this
}

// Icon (requires font awesome)

UI.Icon = function (value) {
  UI.Panel.call(this)

  var dom = document.createElement('span')
  dom.className = 'Icon fa'

  this.dom = dom

  if (value) this.addClass.apply(this, arguments)

  return this
}

UI.Icon.prototype = Object.create(UI.Panel.prototype)

UI.Icon.prototype.hasClass = function (value) {
  return this.dom.classList.contains('fa-' + value)
}

UI.Icon.prototype.addClass = function (value) {
  for (var i = 0; i < arguments.length; i++) {
    this.dom.classList.add('fa-' + arguments[ i ])
  }

  return this
}

UI.Icon.prototype.setClass = function (value) {
  this.dom.className = 'Icon fa'

  for (var i = 0; i < arguments.length; i++) {
    this.dom.classList.add('fa-' + arguments[ i ])
  }

  return this
}

UI.Icon.prototype.removeClass = function (value) {
  for (var i = 0; i < arguments.length; i++) {
    this.dom.classList.remove('fa-' + arguments[ i ])
  }

  return this
}

UI.Icon.prototype.switchClass = function (newValue, oldValue) {
  this.removeClass(oldValue, newValue)
  this.addClass(newValue)

  return this
}

// Toggle Icon

UI.ToggleIcon = function (value, classTrue, classFalse) {
  UI.Icon.call(this, value ? classTrue : classFalse)

  this.value = value
  this.classTrue = classTrue
  this.classFalse = classFalse

  return this
}

UI.ToggleIcon.prototype = Object.create(UI.Icon.prototype)

UI.ToggleIcon.prototype.setValue = function (value) {
  this.value = value

  if (value) {
    this.switchClass(this.classTrue, this.classFalse)
  } else {
    this.switchClass(this.classFalse, this.classTrue)
  }

  return this
}

UI.ToggleIcon.prototype.getValue = function () {
  return this.value
}

// Dispose Icon

UI.DisposeIcon = function () {
  UI.Icon.call(this, 'trash-o')

  var flag = false
  var scope = this

  this.setTitle('delete')
  this.setCursor('pointer')

  this.onClick(function () {
    if (flag === true) {
      if (typeof scope.disposeFunction === 'function') {
        scope.disposeFunction()
      }
    } else {
      scope.setColor('rgb(178, 34, 34)')
      scope.dom.classList.add('deleteInfo')
      flag = true

      setTimeout(function () {
        scope.setColor('#888')
        scope.dom.classList.remove('deleteInfo')
        flag = false
      }, 1500)
    }
  })

  return this
}

UI.DisposeIcon.prototype = Object.create(UI.Icon.prototype)

UI.DisposeIcon.prototype.setDisposeFunction = function (fn) {
  this.disposeFunction = fn

  return this
}

// Progress

UI.Progress = function (max, value) {
  UI.Element.call(this)

  var dom = document.createElement('progress')
  dom.className = 'Progress'

  dom.max = max || 1.0
  if (value !== undefined) dom.value = value

  this.dom = dom

  return this
}

UI.Progress.prototype = Object.create(UI.Element.prototype)

UI.Progress.prototype.getValue = function () {
  return this.dom.value
}

UI.Progress.prototype.setValue = function (value) {
  this.dom.value = value

  return this
}

UI.Progress.prototype.setMax = function (value) {
  this.dom.max = value

  return this
}

UI.Progress.prototype.setIndeterminate = function () {
  this.dom.removeAttribute('value')

  return this
}

// Range

UI.Range = function (min, max, value, step) {
  UI.Element.call(this)

  var dom = document.createElement('input')
  dom.className = 'Range'
  dom.type = 'range'

  dom.min = min.toPrecision(3)
  dom.max = max.toPrecision(3)
  dom.value = value.toPrecision(3)
  dom.step = step.toPrecision(3)

  this.dom = dom

  return this
}

UI.Range.prototype = Object.create(UI.Element.prototype)

UI.Range.prototype.getValue = function () {
  return parseFloat(this.dom.value)
}

UI.Range.prototype.setRange = function (min, max) {
  this.dom.min = min
  this.dom.max = max

  return this
}

UI.Range.prototype.setValue = function (value) {
  this.dom.value = value

  return this
}

UI.Range.prototype.setStep = function (value) {
  this.dom.step = value

  return this
}

// AdaptiveTextArea

UI.AdaptiveTextArea = function () {
  // http://www.brianchu.com/blog/2013/11/02/creating-an-auto-growing-text-input/

  UI.Element.call(this)

  var container = document.createElement('div')
  container.className = 'AdaptiveTextAreaContainer'

  var textarea = document.createElement('textarea')
  textarea.className = 'AdaptiveTextArea'

  var size = document.createElement('div')
  size.className = 'AdaptiveTextAreaSize'

  container.appendChild(textarea)
  container.appendChild(size)

  textarea.addEventListener('input', function (event) {
    size.innerHTML = textarea.value + '\n'
  }, false)

  this.textarea = textarea
  this.size = size
  this.dom = container

  return this
}

UI.AdaptiveTextArea.prototype = Object.create(UI.Element.prototype)

UI.AdaptiveTextArea.prototype.getValue = function () {
  return this.textarea.value
}

UI.AdaptiveTextArea.prototype.setValue = function (value) {
  this.textarea.value = value
  this.size.innerHTML = value + '\n'

  return this
}

UI.AdaptiveTextArea.prototype.setSpellcheck = function (value) {
  this.textarea.spellcheck = value

  return this
}

UI.AdaptiveTextArea.prototype.setBackgroundColor = function (value) {
  this.textarea.style.backgroundColor = value

  return this
}

// Virtual List

UI.VirtualList = function (items, itemHeight, height, generatorFn) {
  // based on Virtual DOM List
  // https://github.com/sergi/virtual-list
  // The MIT License (MIT)
  // Copyright (C) 2013 Sergi Mansilla

  UI.Element.call(this)

  items = items || []
  itemHeight = itemHeight || 20
  height = height || 300
  generatorFn = generatorFn || function () {}

  var dom = document.createElement('div')
  dom.className = 'VirtualList'
  dom.style.height = height + 'px'

  var totalRows = items.length
  var screenItemsCount = Math.ceil(height / itemHeight)
  var cachedItemsCount = screenItemsCount * 3
  var lastRepaintY
  var maxBuffer = screenItemsCount * itemHeight
  var lastScrolled = 0
  var renderChunkCallback = function () {}

  var list = document.createElement('div')
  list.style.width = '100%'
  list.style.height = height + 'px'
  list.style[ 'overflow-y' ] = 'auto'
  list.style.position = 'relative'
  list.style.padding = 0

  var scroller = document.createElement('div')
  scroller.style.opacity = 0
  scroller.style.position = 'absolute'
  scroller.style.top = 0
  scroller.style.left = 0
  scroller.style.width = '1px'
  scroller.style.height = (itemHeight * totalRows) + 'px'

  function createRow (i) {
    var item = generatorFn(i)
    item.classList.add('VirtualListRow')
    item.style.height = itemHeight + 'px'
    item.style.top = (i * itemHeight) + 'px'
    return item
  }

  function renderChunk (from) {
    var finalItem = Math.min(totalRows, from + cachedItemsCount)
    renderChunkCallback(from, finalItem)
    // Append all the new rows in a document fragment
    // that we will later append to the parent node
    var fragment = document.createDocumentFragment()
    for (var i = from; i < finalItem; i++) {
      fragment.appendChild(createRow(i))
    }
    // Hide and mark obsolete nodes for deletion.
    for (var j = 1, l = list.childNodes.length; j < l; j++) {
      list.childNodes[ j ].style.display = 'none'
      list.childNodes[ j ].setAttribute('data-rm', '1')
    }
    list.appendChild(fragment)
  };

  // As soon as scrolling has stopped, this interval asynchronously
  // removes all the nodes that are not used anymore
  var rmNodeInterval = setInterval(function () {
    // check if list is still attached to dom
    var element = dom
    while (element !== document && element.parentNode) {
      element = element.parentNode
    }
    // if list not attached to dom, clear interval
    if (element !== document) {
      clearInterval(rmNodeInterval)
    }
    // remove tagged nodes
    if (Date.now() - lastScrolled > 100) {
      var badNodes = list.querySelectorAll('[data-rm="1"]')
      for (var i = 0, l = badNodes.length; i < l; i++) {
        list.removeChild(badNodes[ i ])
      }
    }
  }, 500)

  function onScroll (e) {
    var scrollTop = e.target.scrollTop // Triggers reflow
    if (!lastRepaintY || Math.abs(scrollTop - lastRepaintY) > maxBuffer) {
      var first = Math.floor(scrollTop / itemHeight)
      renderChunk(Math.max(0, first - screenItemsCount))
      lastRepaintY = scrollTop
    }
    lastScrolled = Date.now()
    e.preventDefault && e.preventDefault()
  }

  // API

  this.setItems = function (value) {
    items = value
    totalRows = items.length
    scroller.style.height = (itemHeight * totalRows) + 'px'
    renderChunk(0)
    return this
  }

  this.setItemHeight = function (value) {
    itemHeight = value
    screenItemsCount = Math.ceil(height / itemHeight)
    cachedItemsCount = screenItemsCount * 3
    maxBuffer = screenItemsCount * itemHeight
    scroller.style.height = (itemHeight * totalRows) + 'px'
    renderChunk(0)
    return this
  }

  this.setHeight = function (value) {
    UI.Element.prototype.setHeight.call(this, value + 'px')
    height = value
    screenItemsCount = Math.ceil(height / itemHeight)
    cachedItemsCount = screenItemsCount * 3
    maxBuffer = screenItemsCount * itemHeight
    list.style.height = height + 'px'
    scroller.style.height = height + 'px'
    renderChunk(0)
    return this
  }

  this.setGeneratorFn = function (value) {
    generatorFn = value
    renderChunk(0)
    return this
  }

  this.setRenderChunkCallback = function (value) {
    renderChunkCallback = value
  }

  this.redraw = function () {
    var first = Math.floor(list.scrollTop / itemHeight)
    renderChunk(Math.max(0, first - screenItemsCount))
    lastRepaintY = list.scrollTop
    return this
  }

  //

  list.appendChild(scroller)
  dom.appendChild(list)
  list.addEventListener('scroll', onScroll)
  renderChunk(0)

  this.dom = dom

  return this
}

UI.VirtualList.prototype = Object.create(UI.Element.prototype)

// Virtual Table

UI.VirtualTable = function (items, itemHeight, height, columns, params) {
  var p = params || {}

  UI.Panel.call(this)

  // this.setOverflow( "scroll" );

  var defaultWidth = p.defaultWidth !== undefined ? p.defaultWidth : 30
  var defaultMargin = p.defaultMargin !== undefined ? p.defaultMargin : 5
  var defaultAlign = p.defaultAlign !== undefined ? p.defaultAlign : 'left'
  var onRowSelect = p.onRowSelect

  // header

  var header = new UI.Panel()
    .setWhiteSpace('nowrap')
    .setDisplay('inline-block')
    .setOverflow('')
    .setWidth('100%')

  var fullWidth = 0

  var selected = []

  var numericalSort = function (a, b) {
    return a - b
  }

  var lexicalSort = function (a, b) {
    return a.localeCompare(b)
  }

  var sortColumn = function (idx, flag) {
    var sort
    if (typeof items[ 0 ][ idx ] === 'string') {
      sort = lexicalSort
    } else {
      sort = numericalSort
    }
    items.sort(function (a, b) {
      if (flag) {
        return sort(b[ idx ], a[ idx ])
      } else {
        return sort(a[ idx ], b[ idx ])
      }
    })
    virtualList.redraw()
    return this
  }

  var selectRow = function (event, idx) {
    selected.length = 0
    if (onRowSelect) onRowSelect(event, idx)
    if (idx !== undefined) {
      selected.push(items[ idx ][ 0 ])
    }
    virtualList.redraw()
    return this
  }

  columns.forEach(function (col) {
    var width = col.width || defaultWidth
    var margin = col.margin || defaultMargin

    var text = new UI.EllipsisText()
      .setValue(col.name)
      .setWidth(width + 'px')
      .setTextAlign(col.align || defaultAlign)
      .setMarginLeft(margin + 'px')
      .setMarginRight(margin + 'px')
      .setCursor('pointer')
      .onClick(function (e) {
        var flag = col.__sortFlag === 'ASC'
        sortColumn(col.index, flag)
        if (flag) {
          col.__sortFlag = 'DESC'
        } else {
          col.__sortFlag = 'ASC'
        }
      })

    header.add(text)

    fullWidth += width + 2 * margin
  })

  // list

  var generatorFn = function (index) {
    var panel = new UI.Panel()

    columns.forEach(function (col) {
      var value = items[ index ][ col.index ]
      if (col.format) value = col.format(value)

      var width = col.width || defaultWidth
      var margin = col.margin || defaultMargin

      var element
      if (typeof value === 'object') {
        element = value
      } else {
        element = new UI.Text()
          .setValue(value)
      }

      element
        .setWidth(width + 'px')
        .setTextAlign(col.align || defaultAlign)
        .setMarginLeft(margin + 'px')
        .setMarginRight(margin + 'px')
        .onClick(function (event) {
          if (typeof col.onClick === 'function') {
            col.onClick(event, index, value)
          }
        })
        .onMouseOver(function (event) {
          if (typeof col.onMouseOver === 'function') {
            col.onMouseOver(event, index, value)
          }
        })
        .onMouseOut(function (event) {
          if (typeof col.onMouseOut === 'function') {
            col.onMouseOut(event, index, value)
          }
        })

      panel.add(element)
    })

    panel
      .setCursor('pointer')
      .onClick(function (event) {
        selectRow(event, index)
      })

    if (selected.indexOf(items[ index ][ 0 ]) !== -1) {
      panel.dom.classList.add('highlight')
    }

    return panel.dom
  }

  var virtualList = new UI.VirtualList(
    items, itemHeight, height, generatorFn
  ).setWidth((fullWidth + 20) + 'px')

  //

  this.add(header, virtualList)

  // API

  this.header = header
  this.list = virtualList
  this.sortColumn = sortColumn
  this.selectRow = function (idx) {
    selectRow(undefined, idx)
  }

  return this
}

UI.VirtualTable.prototype = Object.create(UI.Panel.prototype)

UI.PopupMenu = function (iconClass, heading, constraintTo) {
  constraintTo = constraintTo || 'scrollParent'

  UI.Panel.call(this)

  var entryLabelWidth = '100px'

  var icon = new UI.Icon(iconClass || 'bars')

  var panel = new UI.OverlayPanel()
    .setDisplay('none')
    .attach(this.dom)

  var xOffset = 0
  var yOffset = 0

  var prevX = 0
  var prevY = 0

  function onMousemove (e) {
    if (prevX === 0) {
      prevX = e.clientX
      prevY = e.clientY
    }
    xOffset += prevX - e.clientX
    yOffset += prevY - e.clientY
    prevX = e.clientX
    prevY = e.clientY
  }

  var headingPanel = new UI.Panel()
    .setBorderBottom('1px solid #555')
    .setMarginBottom('10px')
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

  headingPanel
    .add(
      new UI.Icon('times')
        .setFloat('right')
        .setCursor('pointer')
        .onClick(function () {
          this.setMenuDisplay('none')
        }.bind(this))
    )
    .add(
      new UI.Text(heading)
    )

  panel.add(headingPanel)

  icon.setTitle('menu')
  icon.setCursor('pointer')
  icon.onClick(function (e) {
    if (panel.getDisplay() === 'block') {
      this.setMenuDisplay('none')
      return
    }

    panel.setMaxHeight((window.innerHeight / 1.2) + 'px')
    this.setMenuDisplay('block')

    xOffset = 5
    yOffset = 0

  }.bind(this))

  this.add(icon)

  this.setClass('')
    .setDisplay('inline')

  this.icon = icon
  this.panel = panel
  this.entryLabelWidth = entryLabelWidth

  return this
}

UI.PopupMenu.prototype = Object.create(UI.Panel.prototype)

UI.PopupMenu.prototype.addEntry = function (label, entry) {
  this.panel
    .add(new UI.Text(label)
      // .setWhiteSpace( "nowrap" )
      .setWidth(this.entryLabelWidth))
    .add(entry || new UI.Panel())
    .add(new UI.Break())

  return this
}

UI.PopupMenu.prototype.setEntryLabelWidth = function (value) {
  this.entryLabelWidth = value

  return this
}

UI.PopupMenu.prototype.setMenuDisplay = function (value) {
  this.panel.setDisplay(value)

  if (value !== 'none') this.panel.dom.focus()

  return this
}

UI.PopupMenu.prototype.setIconTitle = function (value) {
  this.icon.setTitle(value)

  return this
}

UI.PopupMenu.prototype.dispose = function () {
  this.panel.dispose()

  UI.Element.prototype.dispose.call(this)
}

// Collapsible Icon Panel

UI.CollapsibleIconPanel = function (iconClass1, iconClass2) {
  UI.Panel.call(this)

  this.dom.className = 'Panel CollapsiblePanel'

  if (iconClass1 === undefined) {
    // iconClass1 = iconClass1 || "plus-square";
    // iconClass2 = iconClass2 || "minus-square";

    iconClass1 = iconClass1 || 'chevron-down'
    iconClass2 = iconClass2 || 'chevron-right'
  }

  this.button = new UI.Icon(iconClass1)
    .setTitle('expand/collapse')
    .setCursor('pointer')
    .setWidth('12px')
    .setMarginRight('6px')
  this.addStatic(this.button)

  var scope = this
  this.button.dom.addEventListener('click', function (event) {
    scope.toggle()
  }, false)

  this.content = document.createElement('div')
  this.content.className = 'CollapsibleContent'
  this.dom.appendChild(this.content)

  this.isCollapsed = false

  this.iconClass1 = iconClass1
  this.iconClass2 = iconClass2

  return this
}

UI.CollapsibleIconPanel.prototype = Object.create(UI.CollapsiblePanel.prototype)

UI.CollapsibleIconPanel.prototype.setCollapsed = function (setCollapsed) {
  if (setCollapsed) {
    this.dom.classList.add('collapsed')

    if (this.iconClass2) {
      this.button.switchClass(this.iconClass2, this.iconClass1)
    } else {
      this.button.addClass('rotate-90')
    }
  } else {
    this.dom.classList.remove('collapsed')

    if (this.iconClass2) {
      this.button.switchClass(this.iconClass1, this.iconClass2)
    } else {
      this.button.removeClass('rotate-90')
    }
  }

  this.isCollapsed = setCollapsed
}

// Color picker (requires FlexiColorPicker)
// https://github.com/DavidDurman/FlexiColorPicker
// https://github.com/zvin/FlexiColorPicker

UI.ColorPicker = function () {
  var scope = this

  UI.Panel.call(this)

  // slider

  this.slideWrapper = new UI.Panel()
    .setClass('slide-wrapper')

  this.sliderIndicator = new UI.Panel()
    .setClass('slide-indicator')

  this.slider = new UI.Panel()
    .setClass('slide')
    .setWidth('25px')
    .setHeight('80px')

  this.slideWrapper.add(
    this.slider,
    this.sliderIndicator
  )

  // picker

  this.pickerWrapper = new UI.Panel()
    .setClass('picker-wrapper')

  this.pickerIndicator = new UI.Panel()
    .setClass('picker-indicator')

  this.picker = new UI.Panel()
    .setClass('picker')
    .setWidth('130px')
    .setHeight('80px')

  this.pickerWrapper.add(
    this.picker,
    this.pickerIndicator
  )

  // event

  var changeEvent = document.createEvent('Event')
  changeEvent.initEvent('change', true, true)

  // finalize

  this.add(
    this.pickerWrapper,
    this.slideWrapper
  )

  this.colorPicker = ColorPicker(

    this.slider.dom,
    this.picker.dom,

    function (hex, hsv, rgb, pickerCoordinate, sliderCoordinate) {
      if (!pickerCoordinate && sliderCoordinate && hsv.s < 0.05) {
        hsv.s = 0.5
        hsv.v = 0.7
        scope.colorPicker.setHsv(hsv)

        return
      }

      ColorPicker.positionIndicators(
        scope.sliderIndicator.dom, scope.pickerIndicator.dom,
        sliderCoordinate, pickerCoordinate
      )

      scope.hex = hex
      scope.hsv = hsv
      scope.rgb = rgb

      if (!scope._settingValue) {
        scope.dom.dispatchEvent(changeEvent)
      }
    }

  )

  this.colorPicker.fixIndicators(
    this.sliderIndicator.dom,
    this.pickerIndicator.dom
  )

  return this
}

UI.ColorPicker.prototype = Object.create(UI.Panel.prototype)

UI.ColorPicker.prototype.setValue = function (value) {
  if (value !== this.hex) {
    this._settingValue = true
    this.colorPicker.setHex(value)
    this._settingValue = false
  }

  return this
}

UI.ColorPicker.prototype.getValue = function () {
  return this.hex
}

module.exports = {
    "UI": UI
}
