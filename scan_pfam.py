#!/usr/bin/env python

from Bio import SeqIO

FASTA_OUT = 'blast_out.fa'

sekwencje = list(SeqIO.parse(FASTA_OUT, format='fasta'))

import urllib.request, urllib.parse, urllib.error, urllib.request, urllib.error, urllib.parse
from Bio import SeqRecord
from io import StringIO
import json

# send a GET request to the server

# print out the results
from pprint import pprint
pdir = lambda x : pprint(dir(x))
pt = lambda x : print(type(x))

class HMMER_wrapper:
    # install a custom handler to prevent following of redirects automatically.
    class SmartRedirectHandler(urllib.request.HTTPRedirectHandler):
        def http_error_302(self, req, fp, code, msg, headers):
            return headers

    OUT_FORMAT = 'tsv'
    # OUT_FORMAT = 'json'
    NUM_RECORDS = 50
    HMM_SERVER_ADDR = 'https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan'

    def __init__(self):
        opener = urllib.request.build_opener(self.SmartRedirectHandler())
        urllib.request.install_opener(opener)

    def create_query_params(self, seqrecord):
        f = StringIO()
        SeqIO.write(seqrecord, f, 'fasta')
        f.seek(0)
        parameters = {
            'hmmdb': 'pfam',
            'seq': '>{}\n{}'.format(seqrecord.name, seqrecord.seq)
        }
        parameters = {
            'hmmdb': 'pfam',
            'seq': '>2abl_A mol:protein length:163  ABL TYROSINE KINASE\nMGPSENDPNLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGEWCEAQTKNGQGWVPSNYITPVNSLEKHSWYHGPVSRNAAEYLLSSGINGSFLVRESESSPGQRSISLRYEGRVYHYRINTASDGKLYVSSESRFNTLAELVHHHSTVADGLITTLHYPAP'
        }
        pprint(parameters)
        enc_params = urllib.parse.urlencode(parameters)
        enc_params = enc_params.encode('ascii')
        # enc_params.encode('utf-8')

        # post the seqrch request to the server
        request = urllib.request.Request(self.HMM_SERVER_ADDR, enc_params)

        # get the url where the results can be fetched from
        results_url = urllib.request.urlopen(request).get('location')

        # modify the range, format and presence of alignments in your results here
        res_params = {
            'format': self.OUT_FORMAT,
            'range': '1,{}'.format(self.NUM_RECORDS)
        }

        # add the parameters to your request for the results
        enc_res_params = urllib.parse.urlencode(res_params)
        self.full_url = results_url.replace('results','download') + '?' + enc_res_params

    def __call__(self, seqrecord):
        assert isinstance(seqrecord, SeqRecord.SeqRecord)
        self.create_query_params(seqrecord)
        results_request = urllib.request.Request(self.full_url)
        data = urllib.request.urlopen(results_request)
        # pprint(str(data.read()).splitlines())
        from Bio.SearchIO import HmmerIO
        from Bio import SearchIO
        # hmmscan3-domtab

        # lol = SearchIO.parse(f, 'hmmscan3-tab')
        # # print(f.read())

        table = str(data.read().decode('utf-8')).splitlines()[1:]
        domains = set()
        for row in table:
            domain, _ = self.parse_row(table[0])
            domains.add(domain)
        return list(domains)
        # f = StringIO()
        # f.write(str(data.read().decode('utf-8')))
        # f.seek(0)
        # print(f.read())
        # lol = SearchIO.parse(f, 'hmmscan3-domtab')
        # for el in lol:
        #     print(el)
        # print(st)

    def parse_row(self, line):
        els = line.split('\t')
        TABLE_FAMILY_ID = 1
        TABLE_EVALUE_ID = 10
        try:
            evalue = float(els[TABLE_EVALUE_ID])
        except ValueError:
            print("Bład wewnętrzny - problem z e-wartością ()".format(els[TABLE_EVALUE_ID]))
            exit(1)
        return els[TABLE_FAMILY_ID], evalue


import argparse

def get_filename():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', help='Fasta filename', required=True)
    return parser.parse_args().file

if __name__ == '__main__':
    # file = get_filename()
    file = FASTA_OUT # TODO delete
    with open(file, 'r') as f:
        blast_results = list(SeqIO.parse(f, 'fasta'))
        HMM = HMMER_wrapper()
        seq2domains = {}
        domainsset = set()
        for seq in blast_results:
            domains = HMM(seq)
            domainsset.update(domains)
            seq2domains[seq.name] = domains
            break
        names = [seq.name for seq in blast_results]
        num_seq = len(names)
        num_domains = len(domainsset)
        res = []
        res = [[None] * num_seq for _ in range(num_domains)]
        res[0] = names
        for k, v in seq2domains.items():
            res[]
