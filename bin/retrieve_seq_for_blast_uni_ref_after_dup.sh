FILES=~/stage-thompson/data/raw/uniprot-blast/results/*.results
for f in $FILES
do
	/biolo/blast/bin/blastdbcmd -entry_batch $f -db /commun/bics/DB-Corentin/uniprot >> $f.fasta
done

FILES=~/stage-thompson/data/raw/refseq-blast/results/*.results
for f in $FILES
do
	/biolo/blast/bin/blastdbcmd -entry_batch $f -db /commun/bics/DB-Corentin/refseq-full-prim/all_refseq_prim.fasta >> $f.fasta
done