taxonID=(60711 9541 9544 9555 9595 9598 9601 61853 9483 30611)
name=(C.sabaeus M.fascicularis M.mulatta P.anubis G.gorilla H.sapiens P.troglodytes P.abellii N.leucogenys C.jacchus O.garnettii)
human=9606
for i in ${taxonID[@]}
do
    curl -X GET "https://lbgi.fr/orthoinspectorv3/api/Eukaryota/species/$human/orthologs/$i" -H  "accept: application/json" > ../data/raw/orthoinspector-json-v2/$i.json
done
