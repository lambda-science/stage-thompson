FILES=~/stage-thompson/data/raw/uniprot-blast/*.results
for f in $FILES
do
	/biolo/mafft/inst/bin/mafft $f > $f.mafft
done
