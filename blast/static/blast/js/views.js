
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
  
  function BlastView(container) {
    $container = $(container);
    
    Store.listen('blast_fail', function() {
      alert('Sorry, something went wrong! Please try again later.');
    });
    
    Store.listen('blast_cancel', function() {
      $container.hide();
    });
    
    Store.listen('blast_response', render);
  }
  
  function render() {
    var data = Store.get('blast_response');
    // TODO: render the Blast response in UI
  }
  
  exports.BlastView = BlastView;
})(this);