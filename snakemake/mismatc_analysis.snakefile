#################################################
####### Analyse des alignement à deletion #######
#################################################
# input:
#      dossier raw/
#
# output:
#       ../data/mismatch-analysis2/correction-pairwise/translation_match.tab
#       ../data/mismatch-analysis2/uniprot-translation-correction/translation_match.tab
#
# processus :
#     1/ Récupération des erreurs
#     2/ Récupération des ID des erreurs
#     2.1 / Récupérer les séquences des ID à mismatch
#     3/ Correspondance Uniprot ID -> Ensembl Transcript
#     4/ Récupération séquences génomique, CDS, info exon
#     5/ Création de l'exon_map
#     ((6/ Création cDNA))
#     7/ Correction mismatch 100% identitié
#     8/ Colocalization mismatch exon/intron
#     9/ Correction mismatch par similarité
#
#################################################
#################################################

# Fichiers en sortie du workflow
rule target:
    input:
        "../data/mismatch-analysis2/Exon_map.tab",
        "../data/mismatch-analysis2/Intron_map.tab"
#############################################################################################
#  5/ Création de l'exon_map
rule create_exon_map:
    input:
        "../data/mismatch-analysis2/transcript_ensembl.tab", "../data/mismatch-analysis2/exon_json_dump.json", "../data/mismatch-analysis2/genomic_all.fasta"
    output:
        "../data/mismatch-analysis2/Exon_map.tab", "../data/mismatch-analysis2/Intron_map.tab"
    message:
        "Generation de l'exon et de l'intron map"
    shell:
        "python ../src/4-Generate_Exon_Map.py {input[0]} {input[1]} {input[2]} {output[0]} {output[1]}"

#  4/ Récupération séquences génomique, CDS, info exon
rule get_genomic_CDS_exon_info:
    input:
        "../data/mismatch-analysis2/transcript_ensembl_corrected2.tab"
    output:
        "../data/mismatch-analysis2/CDS_all.fasta", "../data/mismatch-analysis2/genomic_all.fasta", "../data/mismatch-analysis2/exon_json_dump.json"
    message:
        "Récupération seq génomique, seq CDS et json dump des exons"
    shell:
        "python ../src/3-Retrieve_genomic_CDS_and_Exon.pyled-1.py {input[0]} {output[0]} {output[1]} {output[2]}"
