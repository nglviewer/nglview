/**
 * @author mrdoob / http://mrdoob.com/
 * The MIT License, Copyright &copy; 2010-2016 three.js authors
 */

// TODO changes by Alexander S. Rose
// - more events and properties
// - ctrlKey modifier for Number and Integer
// - UI.Element.prototype.getBox()
// - UI.Element.prototype.dispose()
// - UI.Element.prototype.setTitle()
// - UI.Element.prototype.getStyle()

var UI = {}

UI.Element = function () {}

UI.Element.prototype = {

  setId: function (id) {
    this.dom.id = id

    return this
  },

  setTitle: function (title) {
    this.dom.title = title

    return this
  },

  setClass: function (name) {
    this.dom.className = name

    return this
  },

  setStyle: function (style, array) {
    for (var i = 0; i < array.length; i++) {
      this.dom.style[ style ] = array[ i ]
    }
  },

  getStyle: function (style) {
    return this.dom.style[ style ]
  },

  getBox: function () {
    return this.dom.getBoundingClientRect()
  },

  setDisabled: function (value) {
    this.dom.disabled = value

    return this
  },

  setTextContent: function (value) {
    this.dom.textContent = value

    return this
  },

  dispose: function () {
    this.dom.parentNode.removeChild(this.dom)
  }

}

// properties

var properties = [
  'position', 'left', 'top', 'right', 'bottom', 'width', 'height', 'border',
  'borderLeft', 'borderTop', 'borderRight', 'borderBottom', 'borderColor',
  'display', 'overflow', 'overflowX', 'overflowY', 'margin', 'marginLeft',
  'marginTop', 'marginRight',
  'marginBottom', 'padding', 'paddingLeft', 'paddingTop', 'paddingRight',
  'paddingBottom', 'color', 'backgroundColor', 'opacity', 'fontSize',
  'fontWeight', 'fontStyle', 'fontFamily', 'textTransform', 'cursor',
  'verticalAlign', 'clear', 'float', 'zIndex', 'minHeight', 'maxHeight',
  'minWidth', 'maxWidth', 'wordBreak', 'wordWrap', 'spellcheck',
  'lineHeight', 'whiteSpace', 'textOverflow', 'textAlign', 'pointerEvents'
]

properties.forEach(function (property) {
  var methodSuffix = property.substr(0, 1).toUpperCase() +
                        property.substr(1, property.length)

  UI.Element.prototype[ 'set' + methodSuffix ] = function () {
    this.setStyle(property, arguments)
    return this
  }

  UI.Element.prototype[ 'get' + methodSuffix ] = function () {
    return this.getStyle(property)
  }
})

// events

var events = [
  'KeyUp', 'KeyDown', 'KeyPress',
  'MouseOver', 'MouseOut', 'MouseDown', 'MouseUp', 'MouseMove',
  'Click', 'Change', 'Input', 'Scroll'
]

events.forEach(function (event) {
  var method = 'on' + event

  UI.Element.prototype[ method ] = function (callback) {
    this.dom.addEventListener(event.toLowerCase(), callback.bind(this), false)

    return this
  }
})

// Panel

UI.Panel = function () {
  UI.Element.call(this)

  var dom = document.createElement('div')
  dom.className = 'Panel'

  this.dom = dom
  this.children = []

  return this
}

UI.Panel.prototype = Object.create(UI.Element.prototype)

UI.Panel.prototype.add = function () {
  for (var i = 0; i < arguments.length; i++) {
    this.dom.appendChild(arguments[ i ].dom)
    this.children.push(arguments[ i ])
  }

  return this
}

UI.Panel.prototype.remove = function () {
  for (var i = 0; i < arguments.length; i++) {
    this.dom.removeChild(arguments[ i ].dom)

    var idx = this.children.indexOf(arguments[ i ])
    if (idx !== -1) this.children.splice(idx, 1)
  }

  return this
}

UI.Panel.prototype.clear = function () {
  while (this.dom.children.length) {
    this.dom.removeChild(this.dom.lastChild)
  }

  this.children.length = 0

  return this
}

// Collapsible Panel

UI.CollapsiblePanel = function () {
  UI.Panel.call(this)

  this.dom.className = 'Panel CollapsiblePanel'

  this.button = document.createElement('div')
  this.button.className = 'CollapsiblePanelButton'
  this.dom.appendChild(this.button)

  var scope = this
  this.button.addEventListener('click', function (event) {
    scope.toggle()
  }, false)

  this.content = document.createElement('div')
  this.content.className = 'CollapsibleContent'
  this.dom.appendChild(this.content)

  this.isCollapsed = false

  return this
}

UI.CollapsiblePanel.prototype = Object.create(UI.Panel.prototype)

UI.CollapsiblePanel.prototype.addStatic = function () {
  for (var i = 0; i < arguments.length; i++) {
    this.dom.insertBefore(arguments[ i ].dom, this.content)
  }

  return this
}

UI.CollapsiblePanel.prototype.removeStatic = UI.Panel.prototype.remove

UI.CollapsiblePanel.prototype.clearStatic = function () {
  this.dom.childNodes.forEach(function (child) {
    if (child !== this.content) {
      this.dom.removeChild(child)
    }
  })
}

UI.CollapsiblePanel.prototype.add = function () {
  for (var i = 0; i < arguments.length; i++) {
    this.content.appendChild(arguments[ i ].dom)
  }

  return this
}

UI.CollapsiblePanel.prototype.remove = function () {
  for (var i = 0; i < arguments.length; i++) {
    this.content.removeChild(arguments[ i ].dom)
  }

  return this
}

UI.CollapsiblePanel.prototype.clear = function () {
  while (this.content.children.length) {
    this.content.removeChild(this.content.lastChild)
  }
}

UI.CollapsiblePanel.prototype.toggle = function () {
  this.setCollapsed(!this.isCollapsed)
}

UI.CollapsiblePanel.prototype.collapse = function () {
  this.setCollapsed(true)
}

UI.CollapsiblePanel.prototype.expand = function () {
  this.setCollapsed(false)
}

UI.CollapsiblePanel.prototype.setCollapsed = function (setCollapsed) {
  if (setCollapsed) {
    this.dom.classList.add('collapsed')
  } else {
    this.dom.classList.remove('collapsed')
  }

  this.isCollapsed = setCollapsed
}

// Text

UI.Text = function (text) {
  UI.Element.call(this)

  var dom = document.createElement('span')
  dom.className = 'Text'
  dom.style.display = 'inline-block'
  dom.style.verticalAlign = 'middle'

  this.dom = dom
  this.setValue(text)

  return this
}

UI.Text.prototype = Object.create(UI.Element.prototype)

UI.Text.prototype.setValue = function (value) {
  if (value !== undefined) {
    this.dom.textContent = value
  }

  return this
}

UI.Text.prototype.setName = function (value) {
  this.dom.name = value

  return this
}
// Input

UI.Input = function (value) {
  UI.Element.call(this)

  var scope = this

  var dom = document.createElement('input')
  dom.className = 'Input'
  dom.style.padding = '2px'
  dom.style.border = '1px solid #ccc'

  dom.addEventListener('keydown', function (event) {
    event.stopPropagation()
  }, false)

  this.dom = dom
  this.setValue(value || '')

  return this
}

UI.Input.prototype = Object.create(UI.Element.prototype)

UI.Input.prototype.getValue = function () {
  return this.dom.value
}

UI.Input.prototype.setValue = function (value) {
  this.dom.value = value

  return this
}

UI.Input.prototype.setName = function (value) {
  this.dom.name = value

  return this
}

// TextArea

UI.TextArea = function () {
  UI.Element.call(this)

  var scope = this

  var dom = document.createElement('textarea')
  dom.className = 'TextArea'
  dom.style.padding = '2px'
  dom.style.border = '1px solid #ccc'

  dom.addEventListener('keydown', function (event) {
    event.stopPropagation()
  }, false)

  this.dom = dom

  return this
}

UI.TextArea.prototype = Object.create(UI.Element.prototype)

UI.TextArea.prototype.getValue = function () {
  return this.dom.value
}

UI.TextArea.prototype.setValue = function (value) {
  this.dom.value = value

  return this
}

// Select

UI.Select = function () {
  UI.Element.call(this)

  var scope = this

  var dom = document.createElement('select')
  dom.className = 'Select'
  dom.style.width = '64px'
  dom.style.height = '16px'
  dom.style.border = '0px'
  dom.style.padding = '0px'

  this.dom = dom

  return this
}

UI.Select.prototype = Object.create(UI.Element.prototype)

UI.Select.prototype.setMultiple = function (boolean) {
  this.dom.multiple = boolean

  return this
}

UI.Select.prototype.setOptions = function (options) {
  var selected = this.dom.value

  while (this.dom.children.length > 0) {
    this.dom.removeChild(this.dom.firstChild)
  }

  for (var key in options) {
    var option = document.createElement('option')
    option.value = key
    option.innerHTML = options[ key ]
    this.dom.appendChild(option)
  }

  this.dom.value = selected

  return this
}

UI.Select.prototype.getValue = function () {
  return this.dom.value
}

UI.Select.prototype.setValue = function (value) {
  this.dom.value = value

  return this
}

UI.Select.prototype.setName = function (value) {
  this.dom.name = value

  return this
}

// FancySelect

UI.FancySelect = function () {
  UI.Element.call(this)

  var scope = this

  var dom = document.createElement('div')
  dom.className = 'FancySelect'
  dom.tabIndex = 0 // keyup event is ignored without setting tabIndex

  // Broadcast for object selection after arrow navigation
  var changeEvent = document.createEvent('HTMLEvents')
  changeEvent.initEvent('change', true, true)

  // Prevent native scroll behavior
  dom.addEventListener('keydown', function (event) {
    switch (event.keyCode) {
      case 38: // up
      case 40: // down
        event.preventDefault()
        event.stopPropagation()
        break
    }
  }, false)

  // Keybindings to support arrow navigation
  dom.addEventListener('keyup', function (event) {
    switch (event.keyCode) {
      case 38: // up
      case 40: // down
        scope.selectedIndex += (event.keyCode == 38) ? -1 : 1

        if (scope.selectedIndex >= 0 && scope.selectedIndex < scope.options.length) {
          // Highlight selected dom elem and scroll parent if needed
          scope.setValue(scope.options[ scope.selectedIndex ].value)

          scope.dom.dispatchEvent(changeEvent)
        }

        break
    }
  }, false)

  this.dom = dom

  this.options = []
  this.selectedIndex = -1
  this.selectedValue = null

  return this
}

UI.FancySelect.prototype = Object.create(UI.Element.prototype)

UI.FancySelect.prototype.setOptions = function (options) {
  var scope = this

  var changeEvent = document.createEvent('HTMLEvents')
  changeEvent.initEvent('change', true, true)

  while (scope.dom.children.length > 0) {
    scope.dom.removeChild(scope.dom.firstChild)
  }

  scope.options = []

  for (var i = 0; i < options.length; i++) {
    var option = options[ i ]

    var div = document.createElement('div')
    div.className = 'option'
    div.innerHTML = option.html
    div.value = option.value
    scope.dom.appendChild(div)

    scope.options.push(div)

    div.addEventListener('click', function (event) {
      scope.setValue(this.value)
      scope.dom.dispatchEvent(changeEvent)
    }, false)
  }

  return scope
}

UI.FancySelect.prototype.getValue = function () {
  return this.selectedValue
}

UI.FancySelect.prototype.setValue = function (value) {
  for (var i = 0; i < this.options.length; i++) {
    var element = this.options[ i ]

    if (element.value === value) {
      element.classList.add('active')

      // scroll into view

      var y = element.offsetTop - this.dom.offsetTop
      var bottomY = y + element.offsetHeight
      var minScroll = bottomY - this.dom.offsetHeight

      if (this.dom.scrollTop > y) {
        this.dom.scrollTop = y
      } else if (this.dom.scrollTop < minScroll) {
        this.dom.scrollTop = minScroll
      }

      this.selectedIndex = i
    } else {
      element.classList.remove('active')
    }
  }

  this.selectedValue = value

  return this
}

// Checkbox

UI.Checkbox = function (boolean) {
  UI.Element.call(this)

  var scope = this

  var dom = document.createElement('input')
  dom.className = 'Checkbox'
  dom.type = 'checkbox'

  this.dom = dom
  this.setValue(boolean)

  return this
}

UI.Checkbox.prototype = Object.create(UI.Element.prototype)

UI.Checkbox.prototype.getValue = function () {
  return this.dom.checked
}

UI.Checkbox.prototype.setValue = function (value) {
  if (value !== undefined) {
    this.dom.checked = value
  }

  return this
}

UI.Checkbox.prototype.setName = function (value) {
  this.dom.name = value

  return this
}
// Color

UI.Color = function () {
  UI.Element.call(this)

  var scope = this

  var dom = document.createElement('input')
  dom.className = 'Color'
  dom.style.width = '64px'
  dom.style.height = '16px'
  dom.style.border = '0px'
  dom.style.padding = '0px'
  dom.style.backgroundColor = 'transparent'

  try {
    dom.type = 'color'
    dom.value = '#ffffff'
  } catch (exception) {}

  this.dom = dom

  return this
}

UI.Color.prototype = Object.create(UI.Element.prototype)

UI.Color.prototype.getValue = function () {
  return this.dom.value
}

UI.Color.prototype.getHexValue = function () {
  return parseInt(this.dom.value.substr(1), 16)
}

UI.Color.prototype.setValue = function (value) {
  this.dom.value = value

  return this
}

UI.Color.prototype.setHexValue = function (hex) {
  this.dom.value = '#' + ('000000' + hex.toString(16)).slice(-6)

  return this
}

// Number

UI.Number = function (number) {
  UI.Element.call(this)

  var scope = this

  var dom = document.createElement('input')
  dom.className = 'Number'
  dom.value = '0.00'

  dom.addEventListener('keydown', function (event) {
    event.stopPropagation()

    if (event.keyCode === 13) dom.blur()
  }, false)

  this.min = -Infinity
  this.max = Infinity

  this.precision = 2
  this.step = 1

  this.dom = dom
  this.setValue(number)

  var changeEvent = document.createEvent('HTMLEvents')
  changeEvent.initEvent('change', true, true)

  var distance = 0
  var onMouseDownValue = 0

  var pointer = [ 0, 0 ]
  var prevPointer = [ 0, 0 ]

  var onMouseDown = function (event) {
    event.preventDefault()

    distance = 0

    onMouseDownValue = parseFloat(dom.value)

    prevPointer = [ event.clientX, event.clientY ]

    document.addEventListener('mousemove', onMouseMove, false)
    document.addEventListener('mouseup', onMouseUp, false)
  }

  var onMouseMove = function (event) {
    var currentValue = dom.value

    pointer = [ event.clientX, event.clientY ]

    distance += (pointer[ 0 ] - prevPointer[ 0 ]) - (pointer[ 1 ] - prevPointer[ 1 ])

    var modifier = 50
    if (event.shiftKey) modifier = 5
    if (event.ctrlKey) modifier = 500

    var number = onMouseDownValue + (distance / modifier) * scope.step

    dom.value = Math.min(scope.max, Math.max(scope.min, number)).toFixed(scope.precision)

    if (currentValue !== dom.value) dom.dispatchEvent(changeEvent)

    prevPointer = [ event.clientX, event.clientY ]
  }

  var onMouseUp = function (event) {
    document.removeEventListener('mousemove', onMouseMove, false)
    document.removeEventListener('mouseup', onMouseUp, false)

    if (Math.abs(distance) < 2) {
      dom.focus()
      dom.select()
    }
  }

  var onChange = function (event) {
    var number = parseFloat(dom.value)

    dom.value = isNaN(number) === false ? number : 0
  }

  var onFocus = function (event) {
    dom.style.backgroundColor = ''
    dom.style.borderColor = '#ccc'
    dom.style.cursor = ''
  }

  var onBlur = function (event) {
    dom.style.backgroundColor = 'transparent'
    dom.style.borderColor = 'transparent'
    dom.style.cursor = 'col-resize'
  }

  dom.addEventListener('mousedown', onMouseDown, false)
  dom.addEventListener('change', onChange, false)
  dom.addEventListener('focus', onFocus, false)
  dom.addEventListener('blur', onBlur, false)

  return this
}

UI.Number.prototype = Object.create(UI.Element.prototype)

UI.Number.prototype.getValue = function () {
  return parseFloat(this.dom.value)
}

UI.Number.prototype.setValue = function (value) {
  if (isNaN(value)) {
    this.dom.value = NaN
  } else if (value !== undefined) {
    this.dom.value = value.toFixed(this.precision)
  }

  return this
}

UI.Number.prototype.setRange = function (min, max) {
  this.min = min
  this.max = max

  return this
}

UI.Number.prototype.setPrecision = function (precision) {
  this.precision = precision
  this.setValue(parseFloat(this.dom.value))

  return this
}

UI.Number.prototype.setName = function (value) {
  this.dom.name = value

  return this
}
// Integer

UI.Integer = function (number) {
  UI.Element.call(this)

  var scope = this

  var dom = document.createElement('input')
  dom.className = 'Number'
  dom.value = '0.00'

  dom.addEventListener('keydown', function (event) {
    event.stopPropagation()
  }, false)

  this.min = -Infinity
  this.max = Infinity

  this.step = 1

  this.dom = dom
  this.setValue(number)

  var changeEvent = document.createEvent('HTMLEvents')
  changeEvent.initEvent('change', true, true)

  var distance = 0
  var onMouseDownValue = 0

  var pointer = [ 0, 0 ]
  var prevPointer = [ 0, 0 ]

  var onMouseDown = function (event) {
    event.preventDefault()

    distance = 0

    onMouseDownValue = parseFloat(dom.value)

    prevPointer = [ event.clientX, event.clientY ]

    document.addEventListener('mousemove', onMouseMove, false)
    document.addEventListener('mouseup', onMouseUp, false)
  }

  var onMouseMove = function (event) {
    var currentValue = dom.value

    pointer = [ event.clientX, event.clientY ]

    distance += (pointer[ 0 ] - prevPointer[ 0 ]) - (pointer[ 1 ] - prevPointer[ 1 ])

    var modifier = 50
    if (event.shiftKey) modifier = 5
    if (event.ctrlKey) modifier = 500

    var number = onMouseDownValue + (distance / modifier) * scope.step

    dom.value = Math.min(scope.max, Math.max(scope.min, number)) | 0

    if (currentValue !== dom.value) dom.dispatchEvent(changeEvent)

    prevPointer = [ event.clientX, event.clientY ]
  }

  var onMouseUp = function (event) {
    document.removeEventListener('mousemove', onMouseMove, false)
    document.removeEventListener('mouseup', onMouseUp, false)

    if (Math.abs(distance) < 2) {
      dom.focus()
      dom.select()
    }
  }

  var onChange = function (event) {
    var number = parseInt(dom.value)

    if (isNaN(number) === false) {
      dom.value = number
    }
  }

  var onFocus = function (event) {
    dom.style.backgroundColor = ''
    dom.style.borderColor = '#ccc'
    dom.style.cursor = ''
  }

  var onBlur = function (event) {
    dom.style.backgroundColor = 'transparent'
    dom.style.borderColor = 'transparent'
    dom.style.cursor = 'col-resize'
  }

  dom.addEventListener('mousedown', onMouseDown, false)
  dom.addEventListener('change', onChange, false)
  dom.addEventListener('focus', onFocus, false)
  dom.addEventListener('blur', onBlur, false)

  return this
}

UI.Integer.prototype = Object.create(UI.Element.prototype)

UI.Integer.prototype.getValue = function () {
  return parseInt(this.dom.value)
}

UI.Integer.prototype.setValue = function (value) {
  if (value !== undefined) {
    this.dom.value = value | 0
  }

  return this
}

UI.Integer.prototype.setRange = function (min, max) {
  this.min = min
  this.max = max

  return this
}

UI.Integer.prototype.setName = function (value) {
  this.dom.name = value

  return this
}

UI.Integer.prototype.setStep = function (step) {
  this.step = step

  return this
}
// Break

UI.Break = function () {
  UI.Element.call(this)

  var dom = document.createElement('br')
  dom.className = 'Break'

  this.dom = dom

  return this
}

UI.Break.prototype = Object.create(UI.Element.prototype)

// HorizontalRule

UI.HorizontalRule = function () {
  UI.Element.call(this)

  var dom = document.createElement('hr')
  dom.className = 'HorizontalRule'

  this.dom = dom

  return this
}

UI.HorizontalRule.prototype = Object.create(UI.Element.prototype)

// Button

UI.Button = function (value) {
  UI.Element.call(this)

  var scope = this

  var dom = document.createElement('button')
  dom.className = 'Button'

  this.dom = dom
  this.dom.textContent = value

  return this
}

UI.Button.prototype = Object.create(UI.Element.prototype)

UI.Button.prototype.setLabel = function (value) {
  this.dom.textContent = value

  return this
}

// Helper

UI.MenubarHelper = {

  createMenuContainer: function (name, optionsPanel) {
    var container = new UI.Panel()
    var title = new UI.Panel()
    title.setClass('title')

    title.setTextContent(name)
    title.setMargin('0px')
    title.setPadding('8px')

    container.setClass('menu')
    container.add(title)
    container.add(optionsPanel)

    return container
  },

  createOption: function (name, callbackHandler, icon) {
    var option = new UI.Panel()
    option.setClass('option')

    if (icon) {
      option.add(new UI.Icon(icon).setWidth('20px'))
      option.add(new UI.Text(name))
    } else {
      option.setTextContent(name)
    }

    option.onClick(callbackHandler)

    return option
  },

  createOptionsPanel: function (menuConfig) {
    var options = new UI.Panel()
    options.setClass('options')

    menuConfig.forEach(function (option) {
      options.add(option)
    })

    return options
  },

  createInput: function (name, callbackHandler) {
    var panel = new UI.Panel()
      .setClass('option')

    var text = new UI.Text()
      .setWidth('70px')
      .setValue(name)

    var input = new UI.Input()
      .setWidth('40px')
      .onKeyDown(callbackHandler)

    panel.add(text)
    panel.add(input)

    return panel
  },

  createCheckbox: function (name, value, callbackHandler) {
    var panel = new UI.Panel()
      .setClass('option')

    var text = new UI.Text()
      .setWidth('70px')
      .setValue(name)

    var checkbox = new UI.Checkbox()
      .setValue(value)
      .onClick(callbackHandler)

    panel.add(checkbox)
    panel.add(text)

    return panel
  },

  createDivider: function () {
    return new UI.HorizontalRule()
  }

}

module.exports = {
    "UI": UI
}
