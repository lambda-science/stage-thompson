FILES=~/stage-thompson/data/raw/uniprot-blast/results/*.fasta
for f in $FILES
do
	/biolo/mafft/inst/bin/mafft $f > $f.mafft
done

FILES=~/stage-thompson/data/raw/refseq-blast/results/*.fasta
for f in $FILES
do
	/biolo/mafft/inst/bin/mafft $f > $f.mafft
done
