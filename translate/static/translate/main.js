
(function(exports) {

  var $form;
  var $reset;
  var $submit;
  var $seq;

  function TranslateForm(_$form, subElements) {
    $form = _$form;
    $reset = subElements.reset;
    $submit = subElements.submit;
    $seq = subElements.seq;

    Store.set('seq', '');

    $reset.click(reset);
    $submit.click(submit);

    function reset() {
      $seq.val('');
      Store.set('seq', '');
      return false;
    }

    function submit() {
      Store.set('seq', $seq.val());
      return false;
    }
  }

  exports.TranslateForm = TranslateForm;
})(this);

(function(exports) {

  var data = {};

  var Store = {
    get: function(key) {
      return data[key];
    },
    set: function(key, val) {
      data[key] = val;
    }
  };

  exports.Store = Store;
})(this);
