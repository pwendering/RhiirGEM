#!/bin/bash
########################################################################
# get-kcats.sh organismName(opt)
# Retrieves turnover rates for all enzymes for a given organism or the 
# entire set contained in the biochemical databases BRENDA and SABIO-RK
# Author: Philipp Wendering, Bioinformatics, Universitaet Potsdam
########################################################################


# ------------------------------- Input ------------------------------ #
echo ""
printf %"$COLUMNS"s |tr " " "-"

# absolute path of output file destination
outDir="iri1572/data/kcats"

# organism name
organismName=$1

# combined output file
if [[ -z $organismName ]]; then
	organismName="-"
	outFile=$outDir/kcats.tsv
else
	outFile=$outDir/kcats-$(echo $organismName | sed 's/ /-/g').tsv
fi

# ------------------------------ BRENDA ------------------------------ #
# retrieve kcats from BRENDA
echo -e "\nretrieving all available kcats from BRENDA..."

# result file for BRENDA request
if [[ $organismName == "-" ]]; then
	tmpOutFile=$outDir/kcats-BRENDA.tsv
else
	tmpOutFile=$outDir/kcats-BRENDA-$(echo $organismName | sed 's/ /-/g').tsv
fi

# use BRENDA API to get turnover numbers for the given organism or the entire database
getKcatsBRENDA.py $organismName $tmpOutFile

zcat $tmpOutFile.gz > $tmpOutFile

# exclude results for mutants and extract temperature information from commentary column
grep -v mutant $tmpOutFile \
	| awk -F"\t" -v OFS="\t" '{\
		if ($5 ~ /&deg;C/) {\
			match($5,/[0-9\.]+&deg;C/);\
			print $1,$2,$3,$4,substr($5,RSTART,RLENGTH)}\
		else if (NR==1)\
			print;\
		}' \
	| sed 's/&deg;C//g' \
	| sed 's/commentary/T/'\
	> tmp; mv tmp $tmpOutFile

# select rows, which only contain a temperature in the 5th column
awk -F"\t" '{if ((NR==1||$4 ~ /^[0-9\.]+$/)&& $5 ~ /^[0-9\.]+$/) print}' $tmpOutFile \
	> tmp; mv tmp $tmpOutFile

# initialize combined output file
mv $tmpOutFile $outFile

printf %"$COLUMNS"s |tr " " "-"

# ----------------------------- SABIO-RK ----------------------------- #
# retrieve kinetic parameters from SABIO-RK
echo -e "\nretrieving kinetic parameters from SABIO-RK..."

# result file for SABIO-RK request
if [[ $organismName == "-" ]]; then
	tmpOutFile=$outDir/kcats-SABIORK.tsv
else
	tmpOutFile=$outDir/kcats-SABIORK-$(echo $organismName | sed 's/ /-/g').tsv
fi

# use SABIO-RK API to get turnover numbers for the given organism or the entire database
getKineticParamsSABIORK.py $organismName $tmpOutFile

zcat $tmpOutFile.gz > $tmpOutFile

# exclude results for mutant enzymes, search for rows with only kcat as
# parameter and adjust format to: 
# EC \t substrate \t organism \t kcat \t temperature
# separate substrates with '|' instead of ';'
grep -v mutant $tmpOutFile \
	| awk -F"\t" -v OFS="\t" '{if ($4 ~ /^kcat$/ && $9 ~ /s\^\(-1\)/ && $10 ~ /^[0-9\.]+$/ && $6 ~ /^[0-9\.]+$/) print $1,$2,$3,$6,$10}' \
	| sed 's/;/|/g' > tmp; mv tmp $tmpOutFile

# ------------------------- combine results -------------------------- #
cat $tmpOutFile >> $outFile
# take the unique set of entries
cat $outFile | sed '1d'| sort | uniq > tmp; mv tmp $outFile
rm $tmpOutFile

# -------------------------- add lineages ---------------------------- #
echo "adding lineages to species names (this step will take a while)"
$finalKcatFile=$(echo $outFile | sed 's/\.tsv//')-full-lineages.tsv
retrieveLineages.py $outFile $finalKcatFile

# add species name at the end of the lineage (would be better with awk)
paste <(cut -f1,2,3 $finalKcatFile) <(cut -f3 $outFile) <(cut -f4,5 $finalKcatFile) --delimiters=";\t" > tmp ; mv tmp $finalKcatFile

printf %"$COLUMNS"s |tr " " "-"
