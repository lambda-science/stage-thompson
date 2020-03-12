FILES=~/stage-thompson/data/raw/uniprot-sequence-v2/*.id.fasta
for f in $FILES
do
	/biolo/mafft/inst/bin/mafft $f > $f.mafft
done
