mkdir -p DEanalysis/
mkdir -p DEanalysis/counts/
for a in mapping/*/*.tsv; do name=$(dirname $a | cut -f2 -d '/'); sed -e "s/est_counts/$name/g" -e "s/tpm/$name\_tpm/g" $a > DEanalysis/counts/$name.tsv; done
