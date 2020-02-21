# Fichiers en sortie du workflow
rule target:
    input:
        "../data/kirsley-analysis2/kirsley_error_type3_seq.txt",
        "../data/kirsley-analysis2/Exon_map.tab",
        "../data/kirsley-analysis2/Intron_map.tab",
        "../data/kirsley-analysis2/CDS_all_filt.fasta",
        "../data/kirsley-analysis2/mismatch_exon_pos.tab",
        "../data/kirsley-analysis2/mismatch_exon_seq.tab",
        "../data/kirsley-analysis2/mismatch_intron_seq.tab",
        "../data/kirsley-analysis2/mismatch_CDS_seq.tab",
        "../data/kirsley-analysis2/mismatch_genomic_seq.tab",
        "../data/kirsley-analysis2/tblastn/all_couple.txt",
        "../data/kirsley-analysis2/tblastn/match.out"

#############################################################################################
# 8.3/ Processing des résultats tBlastn
rule process_tblastn:
    input:
        "../data/kirsley-analysis2/kirsley_error_type3.txt",
        "../data/kirsley-analysis2/tblastn/all_couple.txt",
        "../data/kirsley-analysis2/tblastn/tblastn.done"
    output:
        "../data/kirsley-analysis2/tblastn/match.out"
    message:
        "Lecture et processing resultats tblastn"
    shell:
    # PENSEZ A CHANGER LE CHEMIN OUTPUT
        "python ../src/kirsley/tblastn_process_9_1.py {output[0]} {input[0]} ../data/kirsley-analysis2/tblastn/"

# 8.2/ Execution du tBlastn
rule exec_tblastn:
    input:
        "../data/kirsley-analysis2/kirsley_error_type3.txt",
        "../data/kirsley-analysis2/genomic_all_filt.fasta",
        "../data/kirsley-analysis2/tblastn/all_couple.txt"
    output:
        "../data/kirsley-analysis2/tblastn/tblastn.done"
    message:
        "Execution des tblastn"
    shell:
    # PENSEZ A CHANGER LE CHEMIN OUTPUT
        "./../bin/tblastn.sh ../data/kirsley-analysis2/"

# 8.1/ Préparation des fichier pour le tBlastn
rule prep_tblastn:
    input:
        "../data/kirsley-analysis2/kirsley_error_type3.txt",
        "../data/kirsley-analysis2/genomic_all_filt.fasta",
    output:
        "../data/kirsley-analysis2/tblastn/all_couple.txt"
    message:
        "Préparation des sequence pour le tblastn"
    shell:
        "python ../src/kirsley/tblastn_seq_9.py {input[0]} {input[1]} ../data/raw/kirsley/ ../data/kirsley-analysis2/tblastn/"

# 7.1/ Colocaliser les mismatch sur la CDS - genomic seq
rule colocalise_on_CDS_genomic:
    input:
        "../data/kirsley-analysis2/mismatch_exon_seq.tab",
        "../data/kirsley-analysis2/mismatch_intron_seq.tab"
    output:
        "../data/kirsley-analysis2/mismatch_CDS_seq.tab",
        "../data/kirsley-analysis2/mismatch_genomic_seq.tab"
    message:
        "Localisation des mismatch au niveau de la CDS et de la sequence genomique"
    shell:
        "python ../src/kirsley/Colocalize_CDS_genomic_mismatch_7_1.py {input[0]} {input[1]} {output[0]} {output[1]}"

# 7/ Colocaliser les mismatch sur les exons/introns
rule colocalise_on_exon_introns:
    input:
        "../data/kirsley-analysis2/kirsley_error_type3.txt", 
        "../data/kirsley-analysis2/CDS_all_filt.fasta",
        "../data/kirsley-analysis2/Exon_map.tab",
        "../data/kirsley-analysis2/Intron_map.tab"
    output:
        "../data/kirsley-analysis2/mismatch_exon_pos.tab",
        "../data/kirsley-analysis2/mismatch_exon_seq.tab",
        "../data/kirsley-analysis2/mismatch_intron_seq.tab"
    message:
        "Localisation des mismatch au niveau des exons et introns"
    shell:
        "python ../src/kirsley/Colocalize_exon_intron_mismatch_7.py {input[0]} {input[1]} {input[2]} {output[0]} {input[3]} {output[1]} {output[2]}"

#  6/ Création de l'exon_map
rule create_exon_map:
    input:
        "../data/kirsley-analysis2/transcript_ensembl_corrected2.tab", "../data/kirsley-analysis2/exon_json_dump.json", "../data/kirsley-analysis2/genomic_all_filt.fasta"
    output:
        "../data/kirsley-analysis2/Exon_map.tab", "../data/kirsley-analysis2/Intron_map.tab"
    message:
        "Generation de l'exon et de l'intron map"
    shell:
        "python ../src/Generate_Exon_Map_4.py {input[0]} {input[1]} {input[2]} {output[0]} {output[1]}"

#  5.2/ Filtrer CDS et genomic
rule filt_cds_genom:
    input:
        "../data/kirsley-analysis2/transcript_ensembl_corrected2.tab", "../data/kirsley-analysis2/CDS_all.fasta", "../data/kirsley-analysis2/genomic_all.fasta"
    output:
        "../data/kirsley-analysis2/CDS_all_filt.fasta", "../data/kirsley-analysis2/genomic_all_filt.fasta"
    message:
        "Filtrer CDS et genomic pour n'avoir aucun dup par ID uniprot"
    shell:
        "cut -f2 {input[0]} | grep -f - -A 1 {input[1]} | grep '\-\-' -v > {output[0]} & " 
        "cut -f2 {input[0]} | grep -f - -A 1 {input[2]} | grep '\-\-' -v > {output[1]}"

# 5.1/ Retirer les dup des ID ensembl
rule filter_cds:
    input:
        "../data/kirsley-analysis2/transcript_ensembl_corrected.tab"
    output:
        "../data/kirsley-analysis2/transcript_ensembl_corrected2.tab"
    message:
        "Retirer les ID ensembl dupliques"
    shell:
        "cat {input} | sort -u -k1,1 | grep -P 'From\tTo' -v | sed '1 i\From\tTo' > {output}"

# 5/ Selection des bonnes CDS/Genomiques ID (long)
rule select_correct_ensemblID:
    input:
        "../data/kirsley-analysis2/transcript_ensembl.tab", "../data/kirsley-analysis2/CDS_all.fasta", "../data/raw/kirsley/all_sequence.fasta"
    output:
        "../data/kirsley-analysis2/transcript_ensembl_corrected.tab"
    message:
        "Selection des ID ensembl correspondant aux protéines"
    shell:
        "python ../src/kirsley/Select_correct_transcript_ensembl.py {input[0]} {input[1]} {input[2]} {output}"

#  4/ Récupération séquences génomique, CDS, info exon
rule get_genomic_CDS_exon_info:
    input:
        "../data/kirsley-analysis2/transcript_ensembl.tab"
    output:
        "../data/kirsley-analysis2/CDS_all.fasta", "../data/kirsley-analysis2/genomic_all.fasta", "../data/kirsley-analysis2/exon_json_dump.json"
    message:
        "Récupération seq génomique, seq CDS et json dump des exons"
    shell:
        "python ../src/Retrieve_genomic_CDS_and_Exon_3.py {input} {output[0]} {output[1]} {output[2]}"

# 3/ Uniprot -> Ensembl
rule get_SEQ_and_ID_errors:
    input:
        "../data/kirsley-analysis2/mismatch.id"
    output:
        "../data/kirsley-analysis2/transcript_ensembl.tab"
    message:
        "Récupération identifiant Trasncript Ensembl pour chaque ID uniprot"
    shell:
        "python ../src/Uniprot_To_Ensembl_1_1.py {input} > {output}"

# 2.1/ Extraction des ID
rule get_ID_errors:
    input:
        "../data/kirsley-analysis2/kirsley_error_type3.txt"
    output:
        "../data/kirsley-analysis2/mismatch.id"
    message:
        "Récupération des ID des erreurs"
    shell:
        "cat {input} | cut -d ' ' -f3 | grep -P -o '^[^_]*' | uniq > {output}"

# 2.1 / Selection des Erreur type 3
rule get_errors_and_seq:
    input:
        "../data/kirsley-analysis2/kirsley_error_type3.txt"
    output:
        "../data/kirsley-analysis2/kirsley_error_type3_seq.txt"
    message:
        "Récupération des sequence dans fichier d'erreur"
    shell:
        "python ../src/kirsley/get_seq_in_error_file.py {input} {output}"

# 2/ Selection des Erreur type 3
rule get_errors:
    input:
        "../data/raw/kirsley_error.txt"
    output:
        "../data/kirsley-analysis2/kirsley_error_type3.txt"
    message:
        "Récupération des erreurs"
    shell:
        "cat {input} | grep \"SEQ_ERROR3\" > {output}"
