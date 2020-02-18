FILES=~/stage-thompson/data/raw/kirsley/*.fasta
for f in $FILES
do
	/biolo/mafft/inst/bin/mafft $f > $f.mafft
done
