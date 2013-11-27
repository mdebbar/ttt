
/*********************************************
 **                 LOADER                  **
 *********************************************/

(function(exports) {
  
  var $loader;
  
  function BlastLoader(container) {
    $loader = $(container);
    
    Store.listen('blast_api', function() {
      $loader.show();
    });

    Store.listen('blast_cancel', 'blast_fail', 'blast_response', function() {
      $loader.hide();
    });
  }
  
  exports.BlastLoader = BlastLoader;
})(this);


/*********************************************
 **               BLAST VIEW                **
 *********************************************/

(function(exports) {
  
  var $container;
  var $alignments;
  
  function BlastView(container, subElems) {
    $container = $(container);
    $alignments = $(subElems.alignments);
    
    Store.listen('blast_fail', function() {
      alert('Sorry, something went wrong! Please try again later.');
    });
    
    Store.listen('blast_cancel', 'blast_fail', function() {
      $container.hide();
    });
    
    Store.listen('blast_response', render);
    
    $alignments.on('click', '.alignment-title', function(event) {
      event.preventDefault();
      $(this).parent().toggleClass('expanded');
    });
    
    $alignments.tooltip({
      delay: 100,
      selector: '.alignment-content span',
      title: function() {
        return String($(this).index());
      }
    });
  }
  
  function render() {
    $container.show();
    
    var record = Store.get('blast_response').records[0];
    var alignments = record.alignments;
    $alignments.html(alignments.map(renderAlignmentRow));
  }
  
  function renderAlignmentRow(alignment) {
    return '<tr><td>' + renderAlignment(alignment) + '</td></tr>';
  }
  
  function renderAlignment(alignment) {
    return '<div class="alignment">' +
      renderTitle(alignment) +
      renderContent(alignment) +
      '</div>';
  }
  
  function renderTitle(alignment) {
    return '<a class="alignment-title" href="#">' + alignment.title + '</a>';
  }
  
  function renderContent(alignment) {
    var hsp = alignment.hsps[0];
    var pieces = ['<div class="alignment-content">'];
    for (var ii = 0; ii < hsp.query.length; ii++) {
      pieces.push(renderHSPPiece(hsp, ii));
    }
    pieces.push('</div>');
    return pieces.join('');
  }
  
  function renderHSPPiece(hsp, ii) {
    var bar = hsp.match[ii] === ' ' ? '&nbsp;' : hsp.match[ii];
    return [
      '<span data-toggle="tooltip">',
      hsp.query[ii],
      bar,
      hsp.subject[ii],
      '</span>'
    ].join('');
  }
  
  exports.BlastView = BlastView;
})(this);