# Script: Retrieve CDS, Genomic and Exon map for HUMAN protein where this is a mismatch w/ primate
# Input:
# 		arg1: ID_File (str) path
# Output:
# 		out1:
# Description:
#
# %%
from src.Retrieve_genomic_CDS_and_Exon_3 import *

# %%
if __name__ == "__main__":
    ID_file = pd.read_csv(sys.argv[1], sep="\t")

    my_response = makeAsyncEnsemblSeqRequest(ID_file, "cds")
    print(my_response)
    writeAsyncEnsemblResponse(my_response, ID_file, sys.argv[2], "CDS")

    my_response2 = makeAsyncEnsemblSeqRequest(ID_file, "genomic")
    print(my_response2)
    writeAsyncEnsemblResponse(my_response2, ID_file, sys.argv[3], "GENOMIC")

    my_response3 = makeAsyncEnsemblExonmapRequest(ID_file)
    print(my_response3)
    writeAsyncEnsemblExonMapResposne(my_response3, sys.argv[4])
