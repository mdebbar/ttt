(function(exports) {
  COMPLEMENTS = {
    A: 'T',
    C: 'G',
    G: 'C',
    T: 'A'
  };
  // Complement
  function cmp(codon) {
    return COMPLEMENTS[codon];
  }

  SEQ_CLEANER = /[^ATGC]/g;
  function clean(seq) {
    return seq.toUpperCase().replace(/U/g, 'T').replace(SEQ_CLEANER, '');
  }

  var MET = 'Met';
  var START = MET;
  var STOP = 'STOP';
  var TRANSLATION_MAP = {
    TTT: 'Phe',  TCT: 'Ser',  TAT: 'Tyr',  TGT: 'Cys',
    TTC: 'Phe',  TCC: 'Ser',  TAC: 'Tyr',  TGC: 'Cys',
    TTA: 'Leu',  TCA: 'Ser',  TAA: STOP ,  TGA: STOP ,
    TTG: 'Leu',  TCG: 'Ser',  TAG: STOP ,  TGG: 'Trp',

    CTT: 'Leu',  CCT: 'Pro',  CAT: 'His',  CGT: 'Arg',
    CTC: 'Leu',  CCC: 'Pro',  CAC: 'His',  CGC: 'Arg',
    CTA: 'Leu',  CCA: 'Pro',  CAA: 'Gln',  CGA: 'Arg',
    CTG: 'Leu',  CCG: 'Pro',  CAG: 'Gln',  CGG: 'Arg',

    ATT: 'Ile',  ACT: 'Thr',  AAT: 'Asn',  AGT: 'Ser',
    ATC: 'Ile',  ACC: 'Thr',  AAC: 'Asn',  AGC: 'Ser',
    ATA: 'Ile',  ACA: 'Thr',  AAA: 'Lys',  AGA: 'Arg',
    ATG: START,  ACG: 'Thr',  AAG: 'Lys',  AGG: 'Arg',

    GTT: 'Val',  GCT: 'Ala',  GAT: 'Asp',  GGT: 'Gly',
    GTC: 'Val',  GCC: 'Ala',  GAC: 'Asp',  GGC: 'Gly',
    GTA: 'Val',  GCA: 'Ala',  GAA: 'Glu',  GGA: 'Gly',
    GTG: 'Val',  GCG: 'Ala',  GAG: 'Glu',  GGG: 'Gly'
  };

  var SHORT_AMINO = {
    ALA : 'A',
    ARG : 'R',
    ASN : 'N',
    ASP : 'D',
    CYS : 'C',
    GLN : 'Q',
    GLU : 'E',
    GLY : 'G',
    HIS : 'H',
    ILE : 'I',
    LEU : 'L',
    LYS : 'K',
    MET : 'M',
    PHE : 'F',
    PRO : 'P',
    SER : 'S',
    THR : 'T',
    TRP : 'W',
    TYR : 'Y',
    VAL : 'V',
    STOP: '#'
  };

  // Reverse Complement
  function revCmp(seq) {
    return seq.split('').reverse().map(cmp).join('');
  }

  // get reading frame starting at `pos`
  function frame(seq, pos) {
    return seq.slice(pos, pos + seq.length - seq.length % 3);
  }

  // translate a seq of 3 codons
  function t3(seq3) {
    return TRANSLATION_MAP[seq3];
  }

  // translate a seq into a list of amino acids
  function t(seq, spec) {
    if (spec.rev) {
      seq = revCmp(seq);
    }
    var res = [];
    for (var ii = spec.pos; ii <= seq.length - 3; ii += 3) {
      res.push(t3(seq.slice(ii, ii + 3)));
    }
    return res;
  }

  function shortenOne(aminoAcid) {
    return SHORT_AMINO[aminoAcid.toUpperCase()];
  }

  function shorten(aminoAcidList) {
    return aminoAcidList.map(shortenOne);
  }


  var FRAME_SPECS = [
    {title: "5'3' Frame 1", pos: 0, rev: false},
    {title: "5'3' Frame 2", pos: 1, rev: false},
    {title: "5'3' Frame 3", pos: 2, rev: false},
    {title: "3'5' Frame 1", pos: 0, rev: true },
    {title: "3'5' Frame 2", pos: 1, rev: true },
    {title: "3'5' Frame 3", pos: 2, rev: true }
  ];

  function processFrames(seq, ii) {
    setTimeout(function() {
      ii = ii || 0;
      var aminos = t(seq, FRAME_SPECS[ii]);
      if (Store.get('options', {}).short) {
        aminos = shorten(aminos);
      }
      Store.update('frames', [{
        title: FRAME_SPECS[ii].title,
        aminos: aminos
      }]);

      if (ii + 1 < FRAME_SPECS.length) {
        processFrames(seq, ii + 1);
      }
    }, 50);
  }

  Store.listen(function() {
    var seq = Store.get('seq');
    if (!seq) {
      Store.set('frames', null);
      return;
    }

    Store.set('frames', []);
    processFrames(clean(seq));
  }, 'seq', 'options');

  exports.Frames = {
    isStart: function(amino) {
      return amino === START || amino === SHORT_AMINO[START.toUpperCase()];
    },
    isStop: function(amino) {
      return amino === STOP || amino === SHORT_AMINO[STOP];
    }
  };
})(this);
