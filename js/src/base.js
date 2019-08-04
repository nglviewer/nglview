var __extends = (this && this.__extends) || function (d, b) {
    for (var p in b) if (b.hasOwnProperty(p)) d[p] = b[p];
    function __() { this.constructor = d; }
    d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
};
var widgets = require("@jupyter-widgets/base");
var BaseView = (function (_super) {
    __extends(BaseView, _super);
    function BaseView() {
        _super.apply(this, arguments);
    }
    BaseView.prototype.render = function () {
        var _this = this;
        this.handleMessage();
        this.displayed.then(function () {
            _this.model.set("_ready", true);
            _this.touch();
        });
    };
    BaseView.prototype.executeCode = function (code) {
        eval(code);
    };
    BaseView.prototype.handleMessage = function () {
        this.model.on("msg:custom", function (msg) {
            this.on_msg(msg);
        }.bind(this));
    };
    BaseView.prototype.on_msg = function (msg) {
        if (msg.type == 'callMethod') {
            this[msg.methodName].apply(this, msg.args, msg.kwargs);
        }
    };
    return BaseView;
})(widgets.DOMWidgetView);
exports.BaseView = BaseView;
