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

    # OUT_FORMAT = 'tsv'
    OUT_FORMAT = 'json'
    NUM_RECORDS = 10
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
            'hmmdb': 'gene3d',
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

        pass
    def __call__(self, seqrecord):
        assert isinstance(seqrecord, SeqRecord.SeqRecord)
        self.create_query_params(seqrecord)
        results_request = urllib.request.Request(self.full_url)
        data = urllib.request.urlopen(results_request)
        # pprint(str(data.read()).splitlines())
        pprint(data.read())

if __name__ == '__main__':
    from Bio import Seq, SeqRecord
    dummySeq = SeqRecord.SeqRecord(Seq.Seq("dsdsd"), name="dsds")
    HMM = HMMER_wrapper()
    HMM(dummySeq)
    # with open(FASTA_OUT, 'r') as f:
    #     blast_results = SeqIO.parse(f, 'fasta')
    #     HMM = HMMER_wrapper()
    #     for el in blast_results:
    #         # print(el)
    #         pt(el)
    #         HMM(el)
    #         break
