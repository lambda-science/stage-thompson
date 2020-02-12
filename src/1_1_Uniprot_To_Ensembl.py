# Script: Convert Uniprot protein ID to Ensembl Transcript ID
# Input:
# 		arg1: File containing uniprot ID one per row
#
# Output:
# 		out1: (stdout) tabulated response from Uniprot

import requests
import pandas as pd
import sys

if __name__ == "__main__":
    URL = "https://www.uniprot.org/uploadlists/"

    ID_list = [line.rstrip('\n') for line in open(sys.argv[1])]
    params = {
        "from": "ACC+ID",
        "to": "ENSEMBL_TRS_ID",
        "format": "tab",
        "query": ' '.join(ID_list)
    }

    r = requests.post(URL, data=params)
    print(r.text)
