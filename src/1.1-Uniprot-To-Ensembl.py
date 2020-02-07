import requests
import pandas as pd
import sys
URL = "https://www.uniprot.org/uploadlists/"

ID_list = [line.rstrip('\n') for line in open(sys.argv[1])]
# defining a params dict for the parameters to be sent to the API
params = {
    "from": "ACC+ID",
    "to": "ENSEMBL_TRS_ID",
    "format": "tab",
    "query": ' '.join(ID_list)
}

# sending get request and saving the response as response object
r = requests.post(URL, data=params)
print(r.text)
