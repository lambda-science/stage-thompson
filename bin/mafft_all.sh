FILES=~/stage-thompson/raw/uniprot-sequence/*.fasta
for f in $FILES
do
	/biolo/mafft/inst/bin/mafft $f > $f.mafft
done
