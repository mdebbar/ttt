import json, os, sys
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.SeqRecord import SeqRecord
from django.conf import settings

E_VALUE_THRESH = 0.04

DATABASES = {
  'blastn' : 'nt',
  'blastp' : 'nr',
  'blastx' : 'nr',
  'tblastn': 'nr',
  'tblastx': 'nt',
}


def sample():
  return fromFileHandle(open(
    os.path.join(settings.ROOT_DIR, 'blast/blastn.xml')
  ))

# algo: ['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx']
def blast(seq, algo):
  algo = algo or 'blastn' # default algo
  
  seqRec = SeqRecord(Seq(seq))
  result_handle = NCBIWWW.qblast(algo, DATABASES[algo], seqRec.seq)
  return fromFileHandle(result_handle)
  
def fromFileHandle(file_handle):
  # result = {'records': [
  #   {'alignments': [
  #     {
  #     'title': '...',
  #     'hsps': [
  #       {'query': '...', 'match': '...', 'subject': '...'},
  #       # ...
  #     ]}
  #   ]},
  #   # ...
  # ]}
  blast_records = NCBIXML.parse(file_handle)
  records = []
  for rec in blast_records:
    alignments = []
    for alignment in rec.alignments:
      hsps = []
      for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
          hspObj = {
            'query': hsp.query,
            'match': hsp.match,
            'subject': hsp.sbjct
          }
          hsps.append(hspObj)
      hsps and alignments.append({
        'title': alignment.hit_def,
        'length': alignment.length,
        'hsps': hsps
      })
    alignments and records.append({
      'alignments': alignments
    })
  return {'records': records}


if __name__ == '__main__':
  seq = sys.argv[-1]
  algorithm = sys.argv[-2]
  print json.dumps(blast(seq, algorithm))