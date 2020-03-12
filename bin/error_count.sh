FILES=~/stage-thompson/data/raw/uniprot-sequence-v2/*.mafft
for f in $FILES
do
	/home/julie/ALN_UTILS/seqerrs $f >> ~/stage-thompson/data/raw/uniprot_errors_V2.txt
done
