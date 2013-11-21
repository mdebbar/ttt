
/*********************************************
 **                FASTA FILE               **
 *********************************************/

(function(exports) {
  
  var $label;
  var $input;
  var enabled = !!(window.File && window.FileReader);
  
  function FastaFile(label) {
    $label = $(label);
    var input = document.getElementById($label.prop('for'));
    $input = $(input);
    
    if (!enabled) {
      $label.addClass('disabled');
    }
    
    $input.change(fileSelected);
  }

  function fileSelected(event) {
    var file = event.target.files[0];

    if (!file) {
      alert("Failed to load file");
    } else if (!file.type.match('text.*')) {
      alert(file.name + " is not a text file.");
    } else {
      readFileContents(file);
    }
    // this is required so that onchange event fires on every file selection
    this.value = null;
  }
  
  function readFileContents(file) {
    var reader = new FileReader();
    reader.onload = function(event) {
      Store.set('seq', event.target.result);
    };
    reader.readAsText(file);
  }
  
  exports.FastaFile = FastaFile;
})(this);


/*********************************************
 **              TRANSLATE FORM             **
 *********************************************/

(function(exports) {

  var $reset;
  var $submit;
  var $seq;

  function TranslateForm(_, subElements) {
    $reset = $(subElements.reset);
    $submit = $(subElements.submit);
    $seq = $(subElements.seq);
    
    FastaFile(subElements.fastafilelabel);

    Store.listen(update, 'seq');

    $reset.click(reset);
    $submit.click(submit);

    function update(_, seq) {
      $seq.val(seq);
    }

    function reset(event) {
      event.preventDefault();
      Store.set('seq', null);
    }

    function submit(event) {
      event.preventDefault();
      Store.set('seq', $seq.val());
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
  
  FrameView.TOTAL_FRAMES = 6;

  function fill() {
    var frames = Store.get('frames');
    if (!frames || frames.length < FrameView.TOTAL_FRAMES) {
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
    var aminos = Store.get('options', {}).short ?
      Frames.shorten(frame.aminos) : frame.aminos;
    return aminos.map(renderAmino).join('');
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
    var isPlainText = Store.get('options', {}).plaintext;
    return isPlainText ? amino : '<span class="' + klass + '">' + amino + '</span>';
  }

  exports.FrameView = FrameView;
})(this);


/*********************************************
 **                 LOADER                  **
 *********************************************/

(function(exports) {
  
  var $loader;
  var $progress;
  
  function Loader(container, subElems) {
    $loader = $(container);
    $progress = $(subElems.progressbar);
    
    Store.listen(update, 'frames');
  }
  
  function update() {
    var frames = Store.get('frames');
    if (frames && frames.length < FrameView.TOTAL_FRAMES) {
      $loader.show();
      var width = Math.floor(frames.length / FrameView.TOTAL_FRAMES * 100);
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
  var $plaintext;

  function OptionsView(container, options) {
    $container = $(container);
    $short = $(options.short);
    $plaintext = $(options.plaintext);
    
    // initialize options
    Store.set('options', {
      short: options.short.checked,
      plaintext: options.plaintext.checked
    });
    
    handleOptionCheckbox($short, 'short');
    handleOptionCheckbox($plaintext, 'plaintext');
  }
  
  function handleOptionCheckbox($input, optionKey) {
    $input.change(function() {
      var newOptions = {};
      newOptions[optionKey] = this.checked;
      setTimeout(function() {
        Store.update('options', newOptions);
      }, 0);
    });
  }

  exports.OptionsView = OptionsView;
})(this);
