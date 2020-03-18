FILES=~/stage-thompson/data/raw/uniprot-blast/*.fasta
database="callithrix_jacchus_uni.fasta  chlorocebus_sabaeus_uni.fasta macaca_fascicularis_uni.fasta nomascus_leucogenys_uni.fasta  pan_troglodytes_uni.fasta pongo_abelii_uni.fasta gorilla_gorilla_uni.fasta macaca_mulatta_uni.fasta papio_anubis_uni.fasta otolemur_garnettii_uni.fasta"
for f in $FILES
do
	cat $f >> $f.results
    for db in $database
	do
		/biolo/blast/bin/blastp -db /commun/bics/DB-Corentin/uniprot-per-primate/$db -evalue 0.005 -outfmt '6 sseqid sseq' -max_target_seqs 1 -query $f | awk 'BEGIN{FS="\t"; OFS="\n"}{gsub(/-/, "", $2); print ">"$1,$2}' >> $f.results
	done
done
