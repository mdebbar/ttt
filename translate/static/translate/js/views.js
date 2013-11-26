
/*********************************************
 **                FASTA FILE               **
 *********************************************/

(function(exports) {
  
  var $label;
  var $input;
  var $droppable;
  var enabled = !!(window.File && window.FileReader);
  
  function FastaFile(_, elements) {
    $label = $(elements.label);
    $droppable = $(elements.droppable);
    
    var input = document.getElementById($label.prop('for'));
    $input = $(input);
    
    if (!enabled) {
      $label.addClass('disabled');
      return;
    }
    
    $input.change(fileChanged);
    $droppable.on({
      dragenter: handleDragEnter,
      dragover: handleDragOver,
      dragleave: handleDragLeave,
      drop: handleDrop
    });
  }
  
  function fileChanged(event) {
    loadFile(event.target.files[0]);
    // this is required so that onchange event fires on every file selection
    this.value = null;
  }

  function loadFile(file) {
    if (!file) {
      alert("Failed to load file");
    } else if (!file.type.match('text.*')) {
      alert(file.name + " is not a text file.");
    } else {
      readFileContents(file);
    }
  }
  
  function readFileContents(file) {
    var reader = new FileReader();
    reader.onload = function(event) {
      Store.set('seq', event.target.result);
    };
    reader.readAsText(file);
  }
  
  // Drag n Drop handlers
  
  function handleDragEnter(event) {
    try {
      if(event.relatedTarget.nodeType == 3) return;
    } catch(err) {}
    if(event.target === event.relatedTarget) return;
    
    $droppable.addClass('dragging');
  }
  
  function handleDragLeave(event) {
    try {
      if(event.relatedTarget.nodeType == 3) return;
    } catch(err) {}
    if(event.target === event.relatedTarget) return;
    
    $droppable.removeClass('dragging');
  }

  function handleDragOver(event) {
    event.preventDefault();
  }
  
  function handleDrop(event) {
    event.stopPropagation();
    event.preventDefault();
    
    $droppable.removeClass('dragging');
    loadFile(event.originalEvent.dataTransfer.files[0]);
  }
  
  exports.FastaFile = FastaFile;
})(this);


/*********************************************
 **              TRANSLATE FORM             **
 *********************************************/

(function(exports) {

  var $reset;
  var $translate;
  var $blast;
  var $seq;

  function TranslateForm(_, subElements) {
    $reset = $(subElements.reset);
    $translate = $(subElements.translate);
    $blast = $(subElements.blast);
    $seq = $(subElements.seq);
    
    FastaFile(_, {
      label: subElements.fastafilelabel,
      droppable: subElements.seq
    });

    Store.listen('seq', update);

    $reset.click(reset);
    $translate.click(translate);
    $blast.click(blast);

    function update(_, seq) {
      $seq.val(seq);
    }

    function reset(event) {
      event.preventDefault();
      Store.set('seq', null);
    }

    function translate(event) {
      event.preventDefault();
      Store.set('blast_cancel', true);
      Store.set('seq', $seq.val());
    }
    
    function blast(event) {
      event.preventDefault();
      Store.set('frames', null);
      Store.set('blast_api', {
        seq: $seq.val(),
        algo: 'blastn'
      });
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
    Store.listen('frames', 'options', render);
  }
  
  FrameView.TOTAL_FRAMES = 6;
  
  function render() {
    var frames = Store.get('frames');
    if (!frames || frames.length < FrameView.TOTAL_FRAMES) {
      $container.hide();
      return;
    }
    
    $container.show();
    frames.forEach(function(frame, ii) {
      renderOne($frames.eq(ii), frame);
    });
  }

  function renderOne($frame, frame) {
    $frame.find('.panel-title').text(frame.title);
    $frame.find('.panel-body').html(renderFrame(frame));
  }

  function renderFrame(frame) {
    within = false;
    var aminos = Store.get('options', {}).short ?
      Frames.shorten(frame.aminos) :
      frame.aminos;
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
  
  function FramesLoader(container, subElems) {
    $loader = $(container);
    $progress = $(subElems.progressbar);
    
    Store.listen('frames', update);
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
  
  exports.FramesLoader = FramesLoader;
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
    
    Store.listen('frames', hideOrShow);
    
    handleOptionCheckbox($short, 'short');
    handleOptionCheckbox($plaintext, 'plaintext');
  }
  
  function hideOrShow() {
    var frames = Store.get('frames');
    if (frames && frames.length >= FrameView.TOTAL_FRAMES) {
      $container.show();
    } else {
      $container.hide();
    }
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
