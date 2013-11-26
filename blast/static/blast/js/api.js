(function(exports) {
  
  exports.API = {
    blast: function(seq, algo) {
      return $.ajax('/blast/', {
        method: 'POST',
        data: {seq: seq, algo: algo}
      });
    }
  };
  
  Store.listen('blast_api', function(_, requestData) {
    if (!requestData.seq) {
      setTimeout(function() {
        Store.set('blast_cancel', true);
      });
      return;
    }

    API.blast(requestData.seq, requestData.algo)
      .done(function(data) {
        Store.set('blast_response', data);
      })
      .fail(function() {
        Store.set('blast_fail', true);
      });
  });
  
})(this);