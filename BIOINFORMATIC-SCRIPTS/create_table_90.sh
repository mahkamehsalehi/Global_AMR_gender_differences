#!/bin/bash
grep '>' /scratch/project_2008149/DATABASES/resfinder_db/all.fsa | sed 's/>//g' | sort > all_resfinder_genenames
while read i
do
	name=($i)
	awk '{print $2}' $name"_ARGs_90.mapstat" | grep -v '#' | grep -v '[a-z]' > $name"_counts"
	awk '{print $1}' $name"_ARGs_90.mapstat" | grep -v '#' > $name"_genenames"
        paste $name"_genenames" $name"_counts" | sort -s -k1,1  > $name"_sorted"
	# Join the sorted count files to the list of all taxa
	join --nocheck-order -a1 -a2 -e "0"  $name"_sorted" all_resfinder_genenames > $name"_mat.txt"
	echo -e $name > $name"_counts"
        awk '{print $2}' $name"_mat.txt" >> $name"_counts"
        sed 's/^$/0/gi' $name"_counts" > $name"_temp"
        echo -e 'GENE' > uniq_taxa_allsamples_header
	rm -f $name"_counts"
	rm -f $name"_mat.txt"
	rm -f $name"_sorted"
	rm -f  $name"_genenames"
	cat all_resfinder_genenames >> uniq_taxa_allsamples_header
	done < names_all
	paste uniq_taxa_allsamples_header  *_temp > mapstat_90_mat.txt
	rm -f *_temp

