#!/usr/bin/env python

XML_OUT = 'blast_out.xml'
FASTA_OUT = 'blast_out.fa'

from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO, SeqRecord, Seq


from Bio import Entrez
Entrez.email = 'bieganski.m@wp.pl'

def do_things(evalue, identity):
    blastp_cline = NcbiblastpCommandline(query='./input-z2.fasta',
                                     db='nr',
                                     outfmt=5,
                                     out=XML_OUT,
                                     remote=True,)

    stdout, stderr = blastp_cline()

    blast_results = NCBIXML.parse(open(XML_OUT))


    records = []
    for hit in blast_results:
        aligns = hit.alignments
        for al in aligns:
            for hsp in al.hsps:
                if hsp.expect > evalue:
                    continue
                if hsp.identities / hsp.align_length < identity:
                    continue
                records.append((al.hit_id.split('|')[1], al.hit_id))


    handle = Entrez.efetch(db='protein', id=','.join([seq for seq, title in records]), retmode='xml')
    entrez_records = list(Entrez.parse(handle))
    handle.close()
    recs = [SeqRecord.SeqRecord(Seq.Seq(seq['GBSeq_sequence']), title) for seq, title in zip(entrez_records, [t for _, t in records])]
    SeqIO.write(recs, FASTA_OUT, "fasta")

import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--evalue', help='Max evalue', required=False, default='10e-10')
    parser.add_argument('-i', '--identity', help='Min identity', required=False, default='0.9')
    args = parser.parse_args()
    try:
        evalue = float(args.evalue)
        identity = float(args.identity)
    except ValueError:
        print("Złe argumenty wywołania! Mają być floaty")
        exit(1)

    do_things(evalue, identity)
    exit(1)