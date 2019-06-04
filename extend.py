#!/usr/bin/env python

XML_OUT = 'blast_out.xml'
FASTA_OUT = 'blast_out.fa'

from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML

blastp_cline = NcbiblastpCommandline(query='./input-z2.fasta', db='nr', outfmt=5, out=XML_OUT, remote=True)

stdout, stderr = blastp_cline()

blast_results = NCBIXML.parse(open(XML_OUT))

from pprint import pprint
from Bio import SeqIO, SeqRecord, Seq

ppdir = lambda x: pprint(dir(x))
pt = lambda x: print(type(x))
# TODO potrzebne?
#  ungapped https://stackoverflow.com/questions/7430679/how-to-get-ungapped-sequences-from-blast-output
#  https://www.biostars.org/p/287383/

records = []
# subjects = set()
# queries = set()
# matches = set()
# i = 0
for hit in blast_results:
    # records.append(hit)
    aligns = hit.alignments
    ppdir = lambda x: pprint(dir(x))
    pt = lambda x: print(type(x))
    pt(hit.alignments[0])
    for al in aligns:
        for hsp in al.hsps:
            nowySeq = Seq.Seq(hsp.sbjct)
            nowySeqRecord = SeqRecord.SeqRecord(nowySeq, al.title)
            records.append(nowySeqRecord)
            # i += 1
            # subjects.add(hsp.sbjct)
            # queries.add(hsp.query)
            # matches.add(hsp.match)
            # print(hsp.sbjct)

# print(i, len(subjects), len(queries), len(matches), len(queries.intersection(matches)))

SeqIO.write(records, FASTA_OUT, "fasta")
