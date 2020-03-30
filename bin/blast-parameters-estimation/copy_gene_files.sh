while read l; do
  cp ../../data/raw/uniprot-blast/$l ../../data/raw/blastp-parameter-estimation/$l
done < ../../data/raw/500_genes_samples.txt