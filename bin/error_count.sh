FILES=~/stage-thompson/raw/uniprot-sequence/*.mafft
for f in $FILES
do
	/home/julie/ALN_UTILS/seqerrs $f >> ~/stage-thompson/raw/uniprot-exon-map/uniprot_errors.txt
done
