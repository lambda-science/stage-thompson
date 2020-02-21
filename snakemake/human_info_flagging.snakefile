# Fichiers en sortie du workflow
rule target:
    input:
        "../data/mismatch-flagging/human_Exon_map.tab",
        "../data/mismatch-flagging/human_Intron_map.tab",
        "../data/mismatch-flagging/human_CDS_all_filt.fasta",
        "../data/mismatch-flagging/human_mismatch_exon_pos.tab",
        "../data/mismatch-flagging/human_mismatch_exon_seq.tab",
        "../data/mismatch-flagging/human_mismatch_intron_seq.tab",
        "../data/mismatch-flagging/human_mismatch_CDS_seq.tab",
        "../data/mismatch-flagging/human_mismatch_genomic_seq.tab",

#############################################################################################

# 7.1/ Colocaliser les mismatch sur la CDS - genomic seq
rule colocalise_on_CDS_genomic:
    input:
        "../data/mismatch-flagging/human_mismatch_exon_seq.tab",
        "../data/mismatch-flagging/human_mismatch_intron_seq.tab"
    output:
        "../data/mismatch-flagging/human_mismatch_CDS_seq.tab",
        "../data/mismatch-flagging/human_mismatch_genomic_seq.tab"
    message:
        "Localisation des mismatch au niveau de la CDS et de la sequence genomique"
    shell:
        "python ../src/Colocalize_CDS_genomic_mismatch_7_1.py {input[0]} {input[1]} {output[0]} {output[1]}"

# 7/ Colocaliser les mismatch sur les exons/introns
rule colocalise_on_exon_introns:
    input:
        "../data/mismatch-flagging/uniprot_errors_type3.txt", 
        "../data/mismatch-flagging/human_CDS_all_filt.fasta",
        "../data/mismatch-flagging/human_Exon_map.tab",
        "../data/mismatch-flagging/human_Intron_map.tab"
    output:
        "../data/mismatch-flagging/human_mismatch_exon_pos.tab",
        "../data/mismatch-flagging/human_mismatch_exon_seq.tab",
        "../data/mismatch-flagging/human_mismatch_intron_seq.tab"
    message:
        "Localisation des mismatch au niveau des exons et introns"
    shell:
        "python ../src/mismatch-flagging/Colocalize_exon_intron_mismatch_7.py {input[0]} {input[1]} {input[2]} {output[0]} {input[3]} {output[1]} {output[2]}"

#  6/ Création de l'exon_map
rule create_exon_map:
    input:
        "../data/mismatch-flagging/human_transcript_ensembl_corrected2.tab", 
        "../data/mismatch-flagging/human_exon_json_dump.json", 
        "../data/mismatch-flagging/human_genomic_all_filt.fasta"
    output:
        "../data/mismatch-flagging/human_Exon_map.tab", "../data/mismatch-flagging/human_Intron_map.tab"
    message:
        "Generation de l'exon et de l'intron map"
    shell:
        "python ../src/Generate_Exon_Map_4.py {input[0]} {input[1]} {input[2]} {output[0]} {output[1]}"

#  5.2/ Filtrer CDS et genomic
rule filt_cds_genom:
    input:
        "../data/mismatch-flagging/human_transcript_ensembl_corrected2.tab", 
        "../data/mismatch-flagging/human_CDS_all.fasta", 
        "../data/mismatch-flagging/human_genomic_all.fasta"
    output:
        "../data/mismatch-flagging/human_CDS_all_filt.fasta", 
        "../data/mismatch-flagging/human_genomic_all_filt.fasta"
    message:
        "Filtrer CDS et genomic pour n'avoir aucun dup par ID uniprot"
    shell:
        "cut -f2 {input[0]} | grep -f - -A 1 {input[1]} | grep '\-\-' -v > {output[0]} & " 
        "cut -f2 {input[0]} | grep -f - -A 1 {input[2]} | grep '\-\-' -v > {output[1]}"

# 5.1/ Retirer les dup des ID ensembl
rule filter_cds:
    input:
        "../data/mismatch-flagging/human_transcript_ensembl_corrected.tab"
    output:
        "../data/mismatch-flagging/human_transcript_ensembl_corrected2.tab"
    message:
        "Retirer les ID ensembl dupliques"
    shell:
        "cat {input} | sort -u -k1,1 | grep -P 'From\tTo' -v | sed '1 i\From\tTo' > {output}"

# 5/ Selection des bonnes CDS/Genomiques ID (long)
rule select_correct_ensemblID:
    input:
        "../data/mismatch-flagging/human_transcript_ensembl.tab",
        "../data/mismatch-flagging/human_CDS_all.fasta", 
        "../data/raw/uniprot-sequence/all_sequence.fasta"
    output:
        "../data/mismatch-flagging/human_transcript_ensembl_corrected.tab"
    message:
        "Selection des ID ensembl correspondant aux protéines"
    shell:
        "python ../src/Select_correct_transcript_ensembl.py {input[0]} {input[1]} {input[2]} {output}"

#  4/ Récupération séquences génomique, CDS, info exon
rule get_genomic_CDS_exon_info:
    input:
        "../data/mismatch-flagging/human_transcript_ensembl.tab"
    output:
        "../data/mismatch-flagging/human_CDS_all.fasta",
        "../data/mismatch-flagging/human_genomic_all.fasta",
        "../data/mismatch-flagging/human_exon_json_dump.json"
    message:
        "Récupération seq génomique, seq CDS et json dump des exons"
    shell:
        "python ../src/Retrieve_genomic_CDS_and_Exon_3.py {input} {output[0]} {output[1]} {output[2]}"

# 3/ Uniprot -> Ensembl
rule get_SEQ_and_ID_errors:
    input:
        "../data/mismatch-flagging/mismatch_human.id"
    output:
        "../data/mismatch-flagging/human_transcript_ensembl.tab"
    message:
        "Récupération identifiant Trasncript Ensembl pour chaque ID uniprot"
    shell:
        "python ../src/Uniprot_To_Ensembl_1_1.py {input} > {output}"

# 2.1/ Extraction des ID
rule get_ID_errors:
    input:
        "../data/mismatch-flagging/uniprot_errors_type3.txt"
    output:
        "../data/mismatch-flagging/mismatch_human.id"
    message:
        "Récupération des ID des erreurs"
    shell:
        "cat {input} | cut -d ' ' -f1 | sed 's/^\/.*\/.*\/.*\/\(.*\)\.id\.fasta\.mafft/\\1/' | uniq > {output}"

# 2/ Selection des Erreur type 3
rule get_errors:
    input:
        "../data/raw/uniprot_new_errors.txt"
    output:
        "../data/mismatch-flagging/uniprot_errors_type3.txt"
    message:
        "Récupération des erreurs"
    shell:
        "cat {input} | grep \"SEQ_ERROR3\" > {output}"
