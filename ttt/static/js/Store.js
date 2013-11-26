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

  function mergeInto(a, b) {
    for (var k in b) {
      if (b.hasOwnProperty(k)) {
        a[k] = b[k];
      }
    }
  }

  exports.Store = {

    get: function(key, def) {
      return data.hasOwnProperty(key) ? data[key] : def;
    },

    set: function(key, val) {
      data[key] = val;
      this.notify(key, val);
    },

    update: function(key, val, notify) {
      if (!data.hasOwnProperty(key)) {
        throw new Error('Can\'t update `' + key + '` doesn\'t exist');
      }
      var currVal = data[key];
      if (!currVal || typeof currVal !== 'object') {
        throw new Error('Can\'t update `' + key + '` type isn\'t supported', key);
      }

      if (currVal instanceof Array) {
        currVal.push.apply(currVal, val);
      }
      else {
        mergeInto(currVal, val);
      }

      if (notify || typeof notify === 'undefined') {
        this.notify(key, currVal);
      }
    },

    listen: function(/*args..., callback*/) {
      var callback = Array.prototype.pop.call(arguments);
      var keys = Array.prototype.slice.call(arguments, 0);
      
      if (keys.length == 0) {
        keys = [GLOBAL_KEY];
      }
      keys.forEach(attachListener.bind(null, callback));
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

})(this);
