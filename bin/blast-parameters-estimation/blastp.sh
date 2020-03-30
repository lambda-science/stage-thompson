# Base de données Uniprot
FILES=/home/meyer/stage-thompson/data/raw/blastp-parameter-estimation/*.id
database="callithrix_jacchus_uni.fasta  chlorocebus_sabaeus_uni.fasta macaca_fascicularis_uni.fasta nomascus_leucogenys_uni.fasta  pan_troglodytes_uni.fasta pongo_abelii_uni.fasta gorilla_gorilla_uni.fasta macaca_mulatta_uni.fasta papio_anubis_uni.fasta otolemur_garnettii_uni.fasta"
for f in $FILES
do
	cat $f | grep ">" >> $f.results
    for db in $database
	do
		/biolo/blast/bin/blastp -db /commun/bics/DB-Corentin/uniprot-per-primate/$db -num_threads 48 -evalue 0.00001 -outfmt '6 sseqid pident' -max_target_seqs 1 -query $f | awk -v title=macaca_mulatta_uni.fasta '{if ($2>=80) print ">"$1 " " title}' | sed 's/>.*|\(.*\)|/>\1 /' >> $f.results
	done
done