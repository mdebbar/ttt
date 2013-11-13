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

  function allFrames(seq) {
    var rev = revCmp(clean(seq));
    return [
      {title: "5'3' Frame 1", seq: frame(seq, 0)},
      {title: "5'3' Frame 2", seq: frame(seq, 1)},
      {title: "5'3' Frame 3", seq: frame(seq, 2)},
      {title: "3'5' Frame 1", seq: frame(rev, 0)},
      {title: "3'5' Frame 2", seq: frame(rev, 1)},
      {title: "3'5' Frame 3", seq: frame(rev, 2)}
    ];
  }

  // translate a seq of 3 codons
  function t3(seq3) {
    return TRANSLATION_MAP[seq3];
  }

  // translate a seq into a list of amino acids
  function t(seq) {
    var res = [];
    for (var ii = 0; ii <= seq.length - 3; ii += 3) {
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

  function translateFrames(seq) {
    return allFrames(clean(seq)).map(function(fr) {
      tr = t(fr.seq);
      return {
        title: fr.title,
        aminos: tr,
        short: shorten(tr),
        seq: fr.seq
      };
    });
  }

  exports.Frames = {
    translateFrames: translateFrames,
    isStart: function(amino) {
      return amino === START || amino === SHORT_AMINO[START.toUpperCase()];
    },
    isStop: function(amino) {
      return amino === STOP || amino === SHORT_AMINO[STOP];
    }
  };
})(this);
