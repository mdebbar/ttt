
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
    Store.listen(fill, 'frames', 'options');
  }

  function fill() {
    var frames = Store.get('frames');
    if (!frames) {
      $container.hide();
      return;
    }

    $container.show();
    if (frames.length < 6) {
      $container.find('.loader').show();
    } else {
      $container.find('.loader').hide();
    }
    frames.forEach(function(frame, ii) {
      fillOne($frames.eq(ii), frame);
    });
  }

  function fillOne($frame, frame) {
    $frame.find('.panel-title').text(frame.title);
    $frame.find('.panel-body').html(renderFrame(frame));
  }

  function renderFrame(frame) {
    within = false;
    return frame.aminos.map(renderAmino).join('');
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
    return '<span class="' + klass + '">' + amino + '</span>';
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
