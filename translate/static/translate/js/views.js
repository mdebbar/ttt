
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

    Store.listen(update, 'seq');

    $reset.click(reset);
    $submit.click(submit);

    function update(_, seq) {
      $seq.val(seq);
    }

    function reset() {
      Store.set('seq', null);
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

  var $container;
  var $frames;

  function FrameView(_$container) {
    $container = _$container;
    $frames = $container.find('.frame');
    Store.listen(fill, 'seq');
  }

  function fill(_, seq) {
    if (!seq) {
      $container.hide();
      return;
    }

    $container.show();
    Frames.translateFrames(seq).forEach(function(frame, ii) {
      fillOne($frames.eq(ii), frame);
    });
  }

  function fillOne($frame, frame) {
    $frame.find('.panel-title').text(frame.title);
    $frame.find('.panel-body').html(renderFrame(frame));
  }

  function renderFrame(frame) {
    if (Store.get('pref', {}).short) {
      return frame.short.map(renderAmino);
    } else {
      return frame.aminos.map(renderAmino);
    }
  }

  var within = false;
  function renderAmino(amino) {
    var klass = 'amino badge';
    if (Frames.isStart(amino)) {
      within = true;
      klass += ' start';
    } else if (Frames.isStop(amino)) {
      within = false;
      klass += ' stop';
    } else if (within) {
      klass += ' filled';
    }
    return $('<span>', {'class': klass}).text(amino);
  }

  exports.FrameView = FrameView;
})(this);


(function(exports) {

  var $container;
  var $short;
  var $autoUpdate;

  function OptionsView(_$container, options) {
    $container = _$container;
    $short = options.short;
    $autoUpdate = options.autoUpdate;
  }

  exports.OptionsView = OptionsView;
})(this);
