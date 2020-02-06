#################################################
####### Analyse des alignement à deletion #######
#################################################
# input:
#      dossier raw/
#
# output:
#       ../data/deletion-analysis/correction-pairwise/translation_match.tab
#       ../data/deletion-analysis/uniprot-translation-correction/translation_match.tab
#
# processus :
#     1/ Récupération des erreurs
#     2/ Récupération des ID des erreurs
#     2.1 / Récupérer les séquences des ID à mismatch
#     3/ Correspondance Uniprot ID -> Ensembl Transcript
#     4/ Récupération séquences génomique, CDS, info exon
#     5/ Création de l'exon_map
#     6/ Création cDNA
#     7/ Correction mismatch 100% identitié
#     8/ Colocalization mismatch exon/intron
#     9/ Correction mismatch par similarité
#
#################################################
#################################################

# Fichiers en sortie du workflow
# rule target:
#    input:
#        "../data/deletion-analysis/correction-pairwise/translation_match.tab",
#        "../data/deletion-analysis/uniprot-translation-correction/translation_match.tab",

#############################################################################################

# 2.1/ Récupérer les séquences des ID à erreur
# rule get_SEQ_errors:
#    input:
#        "../data/deletion-analysis/mismatch.id"
#    output:
#        "../data/deletion-analysis/test_mkdir/mismatch.id.fasta"
#    message:
#        "Récupération des ID des erreurs"
#    shell:
#        "/biolo/blast/bin/blastdbcmd - entry_batch {{input}} -db / commun/bics/DB-Corentin/uniprot >> {{output}}"

# 2/ Récupération des ID des erreurs
rule get_ID_errors:
    input:
        "../data/deletion-analysis/uniprot_errors_type2.txt"
    output:
        "../data/deletion-analysis/mismatch.id"
    message:
        "Récupération des ID des erreurs"
    shell:
        "cat {input} | cut -d ' ' -f3 | uniq > {output}"

# 1/ Récupération des erreurs TYPE 2
rule get_errors:
    input:
        "../data/raw/uniprot_new_errors.txt"
    output:
        "../data/deletion-analysis/uniprot_errors_type2.txt"
    message:
        "Récupération des erreurs"
    shell:
        "cat {input} | grep \"SEQ_ERROR2\" > {output}"
