FILES=~/stage-thompson/data/raw/refseq-blast/results/*.mafft
for f in $FILES
do
	/home/julie/ALN_UTILS/seqerrs $f >> ~/stage-thompson/data/raw/refseq_blast_errors.txt
done

FILES=~/stage-thompson/data/raw/uniprot-blast/results/*.mafft
for f in $FILES
do
	/home/julie/ALN_UTILS/seqerrs $f >> ~/stage-thompson/data/raw/uniprot_blast_errors.txt
done
