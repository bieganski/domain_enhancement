#!/usr/bin/env python
from Bio import SeqIO

sekwencje = list(SeqIO.parse('./input-z2.fasta', format='fasta'))
# sekwencje

from Bio.Blast.Applications import NcbiblastpCommandline


# inaczej :
from Bio.Blast import NCBIWWW
from Bio import SearchIO

blast_xml = 'alignments.xml'
fasta_string = open('./input-z2.fasta')
blast_out = 'alignments.fa'

# with NCBIWWW.qblast("blastp", "nr", fasta_string.read()) as result_handle:
#     with open(blast_xml, 'w') as xml_file:
#         xml_file.write(result_handle.read())
# parse xml and write to fasta
blast_qresult = SearchIO.read(blast_xml, 'blast-xml')
records = []
for hit in blast_qresult:
    records.append(hit[0].hit)
SeqIO.write(records, blast_out, "fasta")

# blastp_cline = NcbiblastpCommandline(query='./input-z2.fasta', db='nr', outfmt=5, out='my_output.xml', remote=True)
# print(blastp_cline)
# stdout, stderr = blastp_cline()
# print(stdout)
# print(stderr)

