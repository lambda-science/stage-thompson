FILES=~/stage-thompson/temp/isoform/all_isoform.id
/biolo/blast/bin/blastdbcmd -entry_batch $FILES -db /commun/bics/DB-Corentin/uniprot >> $FILES.fasta
