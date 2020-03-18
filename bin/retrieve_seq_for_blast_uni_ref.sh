FILES=~/stage-thompson/data/raw/uniprot-blast/*.id
for f in $FILES
do
	/biolo/blast/bin/blastdbcmd -entry_batch $f -db /commun/bics/DB-Corentin/uniprot >> $f.fasta
done

FILES=~/stage-thompson/data/raw/refseq-blast/*.id
for f in $FILES
do
	/biolo/blast/bin/blastdbcmd -entry_batch $f -db /commun/bics/DB-Corentin/uniprot >> $f.fasta
done
