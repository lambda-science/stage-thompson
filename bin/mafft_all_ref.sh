FILES=~/stage-thompson/data/raw/refseq-blast/*.results
for f in $FILES
do
	/biolo/mafft/inst/bin/mafft $f > $f.mafft
done
