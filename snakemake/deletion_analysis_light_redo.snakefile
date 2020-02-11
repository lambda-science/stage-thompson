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
#     7/ Correction mismatch 100% identitié
#     9/ Correction mismatch par similarité

#################################################
#################################################

rule target:
    input:
        "../data/deletion-analysis/translation_match_100percent.tab",
        "../data/deletion-analysis/translation_match_similarity.tab",

#############################################################################################
# 9/ Correction erreur type2 par traduction et 80% similarité
rule error_correction_similarity_translation:
    input:
        "../data/deletion-analysis/uniprot_errors_type2.txt", "../data/deletion-analysis/genomic_all.fasta"
    params:
        out_folder="../data/deletion-analysis/similarity_correction/",
        mafft_path="/biolo/mafft/inst/bin/mafft",
        relatif="alignement"
    output:
        "../data/deletion-analysis/translation_match_similarity.tab"
    message:
        "Traduction des séquences et recherche du peptide human avec 80% similarité et <10 gap"
    shell:
        "python ../src/8-Mismatch_correction_similarity.py {input[0]} {input[1]} {params.out_folder} {params.mafft_path} {output} {params.relatif}"

# 7/ Correction de l'erreur type2 par traduction et recherche de 100% d'identité
rule error_correction_identity_translation:
    input:
        "../data/deletion-analysis/uniprot_errors_type2.txt", "../data/deletion-analysis/genomic_all.fasta"
    params:
        relatif="alignement"
    output:
        "../data/deletion-analysis/translation_match_100percent.tab"
    message:
        "Traduction des séquences génomique et recherche du peptide humain pour correction"
    shell:
        "python ../src/6-Mismatch_correction_translation.py {input[0]} {input[1]} {output} {params.relatif}"
