# Base de données Uniprot
FILES=/home/meyer/stage-thompson/data/raw/uniprot-blast/*.id
database="callithrix_jacchus_uni.fasta  chlorocebus_sabaeus_uni.fasta macaca_fascicularis_uni.fasta nomascus_leucogenys_uni.fasta  pan_troglodytes_uni.fasta pongo_abelii_uni.fasta gorilla_gorilla_uni.fasta macaca_mulatta_uni.fasta papio_anubis_uni.fasta otolemur_garnettii_uni.fasta"
for f in $FILES
do
	cat $f >> $f.results
    for db in $database
	do
		/biolo/blast/bin/blastp -db /commun/bics/DB-Corentin/uniprot-per-primate/$db -num_threads 48 -evalue 0.005 -outfmt '6 sseqid sseq' -max_target_seqs 1 -query $f | awk -v title=$db 'BEGIN{FS="\t"; OFS="\n"}{gsub(/-/, "", $2); print ">"$1 " " title,$2}' >> $f.results
	done
sed -i 's/>.*|\(.*\)|/>\1 /' $f.results
done

# Base de données RefSeq
FILES=/home/meyer/stage-thompson/data/raw/refseq-blast/*.id
database="callithrix_jacchus_ref  chlorocebus_sabaeus_ref macaca_fascicularis_ref nomascus_leucogenys_ref  pan_troglodytes_ref pongo_abelii_ref gorilla_gorilla_ref macaca_mulatta_ref papio_anubis_ref otolemur_garnettii_ref"
for f in $FILES
do
	cat $f >> $f.results
    for db in $database
	do
		/biolo/blast/bin/blastp -db /commun/bics/DB-Corentin/refseq-per-primate/$db -num_threads 48 -evalue 0.005 -outfmt '6 sseqid sseq' -max_target_seqs 1 -query $f | awk -v title=$db 'BEGIN{FS="\t"; OFS="\n"}{gsub(/-/, "", $2); print ">"$1 " " title,$2}' >> $f.results
	done
sed -i 's/>.*|\(.*\)|/>\1 /' $f.results
done
