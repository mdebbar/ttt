
/*********************************************
 **                FASTA FILE               **
 *********************************************/

(function(exports) {
  
  function FastaFile($container) {
    
  }
  
  exports.FastaFile = FastaFile;
})(this);


/*********************************************
 **              TRANSLATE FORM             **
 *********************************************/

(function(exports) {

  var $form;
  var $reset;
  var $submit;
  var $seq;

  function TranslateForm(form, subElements) {
    $form = $(form);
    $reset = $(subElements.reset);
    $submit = $(subElements.submit);
    $seq = $(subElements.seq);

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

/*********************************************
 **                 FRAMES                  **
 *********************************************/

(function(exports) {

  var $container;
  var $frames;
  
  function FrameView(container) {
    $container = $(container);
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


/*********************************************
 **                 LOADER                  **
 *********************************************/

(function(exports) {
  
  var $loader;
  var $progress;
  
  const TOTAL_FRAMES = 6;
  
  function Loader(container, subElems) {
    $loader = $(container);
    $progress = $(subElems.progressbar);
    
    Store.listen(update, 'frames');
  }
  
  function update() {
    var frames = Store.get('frames');
    if (frames && frames.length < TOTAL_FRAMES) {
      $loader.show();
      var width = Math.floor(frames.length / TOTAL_FRAMES * 100);
      $progress.css({width: String(width) + '%'});
    } else {
      $loader.hide();
    }
  }
  
  exports.Loader = Loader;
})(this);


/*********************************************
 **                 OPTIONS                 **
 *********************************************/

(function(exports) {

  var $container;
  var $short;
  var $autoUpdate;

  function OptionsView(container, options) {
    $container = $(container);
    $short = $(options.short);
    $autoUpdate = $(options.autoUpdate);
  }

  exports.OptionsView = OptionsView;
})(this);
