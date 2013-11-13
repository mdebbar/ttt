(function(exports) {

  var data = {};
  var listeners = {};
  var GLOBAL_KEY = '__';

  function attachListener(cb, key) {
    if (!listeners[key]) {
      listeners[key] = [];
    }
    listeners[key].push(cb);
  }

  var Store = {
    get: function(key, def) {
      return data.hasOwnProperty(key) ? data[key] : def;
    },
    set: function(key, val) {
      data[key] = val;
      this.notify(key, val);
    },
    listen: function(cb) {
      var keys = Array.prototype.slice.call(arguments, 1);
      if (keys.length == 0) {
        keys = [GLOBAL_KEY];
      }
      keys.forEach(attachListener.bind(null, cb));
    },
    notify: function(key, val) {
      if (listeners[key]) {
        listeners[key].forEach(function(cb) {
          cb(key, val);
        });
      }
      if (listeners[GLOBAL_KEY]) {
        listeners[GLOBAL_KEY].forEach(function(cb) {
          cb(key, val);
        });
      }
    }
  };

  exports.Store = Store;
})(this);
