#!/usr/bin/env python

from Bio import SeqIO

# TODO tak naprawde to powinny sekwencje z outputu extend.py
sekwencje = list(SeqIO.parse('./input-z2.fasta', format='fasta'))

import urllib.request, urllib.parse, urllib.error, urllib.request, urllib.error, urllib.parse

# install a custom handler to prevent following of redirects automatically.
class SmartRedirectHandler(urllib.request.HTTPRedirectHandler):
    def http_error_302(self, req, fp, code, msg, headers):
        return headers
opener = urllib.request.build_opener(SmartRedirectHandler())
urllib.request.install_opener(opener)

parameters = {
    'seqdb':'pdb',
    'seq':'>Seq\nKLRVLGYHNGEWCEAQTKNGQGWVPSNYITPVNSLENSIDKHSWYHGPVSRNAAEY'
}
enc_params = urllib.parse.urlencode(parameters)

#post the seqrch request to the server
request = urllib.request.Request('https://www.ebi.ac.uk/Tools/hmmer/search/phmmer', enc_params.encode('ascii'))


#get the url where the results can be fetched from
results_url = urllib.request.urlopen(request).get('location')

# modify the range, format and presence of alignments in your results here
res_params = {
    'output': 'json',
    'range': '1,10'
}

# add the parameters to your request for the results
enc_res_params = urllib.parse.urlencode(res_params)
modified_res_url = results_url + '?' + enc_res_params

# send a GET request to the server
results_request = urllib.request.Request(modified_res_url)
data = urllib.request.urlopen(results_request)

# print out the results
from pprint import pprint
pt = lambda x : print(type(x))
import json

json_obj = json.loads(data.read().decode('utf-8'))
pprint(json_obj)