#!/bin/bash

### version history and description

## current version :
# 1.7 (17/02/2021,  13/01/2021 to 07/04/2021, 28/07/2021, 09/02/2022) : respectively:
	# - polyploid genotyping
	# - making options available to call as missing loci with more than two alleles (-M --max-noise-pct allelefinderr parameter [FDSTools default at 10; previous SSRseq default at 75]) and contaminated individuals ( -x --max-noisy allelefinder parameter [FDSTools default at 2; previous SSrseq default at 10])
	# - Also correcting partial match of allele sequence to count for allele occurance across individuals.
	# - correcting error in allelefinder_cpSSRfiltered file where the header line was missing (ParameterCOmparison scrip only : FinalGenotyping script didn't have the error)
	
## previous versions : 
# 1.6 (04/06/2021) : correcting of a bug when counting the number of mismatch leading to an underestimation of error.
# 1.5 (28/01/2020): adding option for chloroplastic SSR data analysis (filtering of the allele with the higest number of repeat after accounting for strong stuttering, allelic error rate computed based on a single allele called per individuals)
# 1.4 (11/10/19): adapting for use with FDStools 1.2 (-q option previously used by tssv to inticate a fastq input file is not in use anymore in the last FDSTools 1.2
# 1.3 (29/03/2019) : correcting a bug when counting allele occurance due to the sorting routine.
# version : 1.2 (06/02/2019) : correcting a bug in the cuont of the number of alleles for the LocusCharacteristic output file.
# version 1.1 (14/01/2019) : corroige le format de sortie dans le fichier LocusCharacteristics et AlleleInformations afin de ne pas inclure d'espaces additionnels entre les champs mais seulement une tabulation.
# version 1.0 (16/11/20018) : bug corrected and improvement of genotypic table format and loci characteristic report
# version beta (end of october 2018) : first compilation of the script for automatic genotyping and parameters comparison (detection of small bugs due to issue in input file names; identification of output file format to improve)
# version alpha (thoughout 2017 and 2018) : typical script file without automation; was used for data exploration and result generation in a manual way

ANALYSIS=${1:-nSSR}
TYPE=${2:-FullLength}
PARAMSET1=${3:-ParameterSet1}
S1a=${4:-50}
S1b=${5:-10}
M1=${6:-15}
N1=${7:-20}
PARAMSET2=${8:-ParameterSet2}
S2a=${9:-70}
S2b=${10:-10}
M2=${11:-10}
N2=${12:-20}
MNP=${13:-75}
X=${14:-10}
A=${15:-2}


if [ $ANALYSIS = "nSSR" ] && [ $A -gt 2 ]
	then
echo "Polyploid"

# tssv
echo "***************FDSTools tssv analysis***************"
rm -f -r tssvResults; mkdir tssvResults
rm -f -r tssvReports; mkdir tssvReports
for file in ./samples/*.fast[a-q]
do
echo ${file:10:-6}
fdstools tssv --num-threads 4 -R ./tssvReports/${file:10:-6}_tssvReport.txt -m 0.08 ${ANALYSIS}_${TYPE}_FDSTools_InputFile.txt $file ./tssvResults/${file:10:-6}_tssv.txt
done
 
# extracting the number of reads by locus by individuals
echo "***************Computing coverage by locus by individual***************"
rm -f LocusCoverageperIndividual_${ANALYSIS}_${TYPE}.txt
paste <(echo "Locus") <(ls tssvReports | sed "s/_tssvReport.txt//g" | tr "\n" "\t" | sed -e '$a\') > LocusCoverageperIndividual_${ANALYSIS}_${TYPE}.txt
while read locus
do
	echo $locus
	echo $locus > tmpLocusCoverage.txt
	for file in ./tssvReports/*
	do
	#echo $file	
	cat $file | grep -w $locus | cut -f3 >> tmpLocusCoverage.txt
	done
cat tmpLocusCoverage.txt | tr "\n" "\t" | sed -e '$a\' >> LocusCoverageperIndividual_${ANALYSIS}_${TYPE}.txt
done < ${ANALYSIS}_LocusList.txt
rm tmpLocusCoverage.txt

######################### analysis with PARAMSET1###################################

# running stuttermark 
echo "***************FDSTools stuttermark ${PARAMSET1} analysis***************"
rm -f -r stuttermark; mkdir stuttermark
for file in ./tssvResults/*.txt
do
echo ${file:14:-4}
fdstools stuttermark -s=-1:${S1a},+1:${S1b} -m 2 -l ${ANALYSIS}_${TYPE}_FDSTools_InputFile.txt $file ./stuttermark/${file:14:-4}_stuttermarked.txt
done

# running allelefinder 
echo "***************FDSTools allelefinder ${PARAMSET1} analysis***************"
rm -f -r allelefinder; mkdir allelefinder
rm -f -r alellefinderReports; mkdir alellefinderReports
for file in ./stuttermark/*.txt
do
echo ${file:14:-4}
fdstools allelefinder -m ${M1} -n ${N1} -a ${A} -x ${X} -M ${MNP} -c annotation -l ${ANALYSIS}_${TYPE}_FDSTools_InputFile.txt $file -o ./allelefinder/${file:14:-4}_allelefinder.txt -R ./alellefinderReports/${file:14:-4}_AFreport.txt
done


# concatenate all resulting allele across individuals (except the firt line which is a header)
echo "***************Creating a unique genotype file***************"
echo -e "sample\tmarker\treads\tallele" > AllAlleles.txt
awk 'FNR>1' ./allelefinder/*.txt >> AllAlleles.txt

#### convertir les séquence des allèles en code allele standards et formater les genotypes

# creating allele naming convention file
echo "***************Listing allele sequence, allele size, allele code and allele count***************"
rm -f AlleleNamingConvention_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt
while read locus
do
	echo $locus
	nall=$(awk -v l=${locus} '$2==l { print $4 }' AllAlleles.txt | sort | uniq | wc -l)
	echo $nall
	allmax=$(expr 100 + $nall)
	echo $allmax
	all=$(awk -v l=${locus} '$2==l { print $4 }' AllAlleles.txt | sort --version-sort | uniq)	
	paste <(yes $locus | head -n $nall) <(echo "${all}") <(eval echo {101..$allmax} | tr " " "\n") <(for alleleSeq in $(awk -v l=${locus} '$2==l { print $4 }' AllAlleles.txt | sort --version-sort | uniq); do grep -w $locus AllAlleles.txt | cut -f4 | grep -x $alleleSeq | wc -l; done) >> AlleleNamingConvention_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt
done < ${ANALYSIS}_LocusList.txt
wait
sed -i '/^\s\s*/d' AlleleNamingConvention_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt

# renaming the alleles
echo "***************Allele coding***************"
rm -f AlleleNamed.txt
while read locus
do
	echo $locus
	grep ^$locus AlleleNamingConvention_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt | cut -f2,3 | sed -e $'s!\t!$/!' | sed -e 's$^$s/^$' | sed -e 's?$?/g?' > tmp_seddScript.txt
	paste <(awk -v l=$locus '$2==l { print $1, $2 }' AllAlleles.txt) <(awk -v l=$locus '$2==l { print $4 }' AllAlleles.txt | sed -f tmp_seddScript.txt) | sed 's/\t/ /' >> AlleleNamed.txt
done < ${ANALYSIS}_LocusList.txt
rm tmp_seddScript.txt

# clearing the names of individuals
echo "***************Individual name cleaning***************"
sed -i -E 's%./stuttermark/%%g' AlleleNamed.txt
sed -i -E 's%_tssv_stuttermarked%%g' AlleleNamed.txt
sed -i -E 's%.assembled%%g' AlleleNamed.txt

# formating allele information file
echo "***************Retrieving allele information across locus***************"
rm -f AllAlleleSeq.txt
while read locus seq code n
do
counter=1
rm -f tmp_alleleSeq.txt
stg=$(echo $seq | sed "s/(/\t/g" | sed "s/)/\t/g" | sed "s/ /\t/g")
while [ $counter -le $(echo $stg | awk '{print NF}' | sort -nu | tail -n 1) ]
do
	seq=$(echo $stg | cut -d ' ' -f$counter)
	c2=$(expr $counter + 1)
	rep=$(echo $stg | cut -d ' ' -f$c2)
	yes $seq | head -n $rep | tr "\n" " " | sed "s/ //g" | sed -e '$a\' >> tmp_alleleSeq.txt
	((counter++))
	((counter++))
done
paste <(cat tmp_alleleSeq.txt | tr "\n" " " | sed "s/ //g" | sed -e '$a\') <(echo $(cat tmp_alleleSeq.txt | tr "\n" " " | sed "s/ //g" | sed -e '$a\' | wc -m)-1 | bc) >> AllAlleleSeq.txt
done < AlleleNamingConvention_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt

echo -e "Locus\tAlleleSequenceAnnotated\tAlleleSeqCode\tOccurancesAcrossIndivs\tAlleleSequence\tAlleleLength" > AlleleInformationFile_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt
paste <(cat AlleleNamingConvention_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt) <(cat AllAlleleSeq.txt) >> AlleleInformationFile_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt

rm AlleleNamingConvention_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt

# formating the genotypic table across all individuals / locus
echo "***************Genotypic table formating***************"
rm -f tmp.txt tmp_indiv.txt
paste <(echo "Indiv") <(cut -f1 ${ANALYSIS}_LocusList.txt | tr "\n" "\t" | sed "s/\t/$(for i in $(seq 1 ${A}); do printf "\t"; done)/g") > GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt
for file in ./samples/*.fast[a-q]
do
	ind=${file:10:-6}
	echo $ind > tmp_indiv.txt
	echo $ind
	while read locus
	do
	grep -w "^${ind}" AlleleNamed.txt | grep -w "$locus" > tmp.txt
	lines=$(wc -l tmp.txt | cut -f1 -d" ")
	if [ $lines -gt 0 ]
	then
	alleles=$(echo -e $(cat tmp.txt | cut -f3 -d" " | sort | tr "\n" "\t") $(for i in $(seq 1 $((A-lines))); do printf "\t."; done)) 
	elif [ $lines = 0 ]
	then	
	alleles=$(echo -e "NA" $(for i in $(seq 1 $((A-1))); do printf "\tNA"; done))
	else
	alleles=$(echo -e "Err" $(for i in $(seq 1 $((A-1))); do printf "\tErr"; done))
	fi
	echo $alleles >> tmp_indiv.txt
	done < ${ANALYSIS}_LocusList.txt
	cat tmp_indiv.txt | tr "\n" "\t" | sed -e '$a\' >> GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt
done
sed -i "s/ /\t/g" GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt
sed -i "s/_S[0-9]*\t/\t/g" GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt
rm tmp.txt tmp_indiv.txt


### Summarizing locus characteristics and comparing repeated genotypes to compute genotypic error rate 
echo "***************Comparing repeated genotypes to computing allelic error rate, and retrieving the number of alleles and missing data rate across the whole genotypic table ***************"
echo -e "Loci\tNallelesSequence\tNalleleSize\tMissingRate\tAllelicMismatches\tPolyploidGenotypesComparred\tGenotypicError" > LocusCharacteristics_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt
for loc in $(eval echo "{1..$(wc -l < ${ANALYSIS}_LocusList.txt)}")
do
	echo $loc
	locName=$(head GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt -n 1 | sed "s/$(for i in $(seq 1 ${A}); do printf "\t"; done)/\t/g" | cut -f$(($loc + 1)))
	echo $locName
	## Number of alleles based on sequence
	NallSeq=$(grep -w $locName AlleleInformationFile_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt | cut -f3 | uniq | wc -l)
	## Number of alleles based on size
	NallSize=$(grep -w $locName AlleleInformationFile_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt | cut -f6 | sort -g |  uniq |wc -l)
	## Missing data
	cl1=$(echo "($loc-1)*${A}+2" | bc )
	cl2=$(echo "$cl1+(${A}-1)" | bc )
	Nmiss=$(echo "scale=5; ($(sed 1d GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt | cut -f$cl1 | grep "NA" | wc -l) / $(sed 1d GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt | cut -f$cl1 | wc -l))*100" | bc -l 2> /dev/null)%
	## Error rates
	tmperr=0
	compar=0
	paste <(sed 1d GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt | cut -f1) <(sed 1d GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt | cut -f$cl1-$cl2) > tmpTable.txt
	cat tmpTable.txt | cut -f1 | sed "s/_[A-B]$//g" | uniq -d | (while read ind
	do
	#echo $ind	
	na=$(grep ${ind}_[A-B] tmpTable.txt | cut -f2-$(($A+1)) | grep -v "NA" | wc -l)
	if [[ -n $na ]] && [[ $na > 1 ]]
	then
	match1=$(grep ${ind}_A tmpTable.txt | cut -f2-$(($A+1)) | grep -v "NA")
	match2=$(grep ${ind}_B tmpTable.txt | cut -f2-$(($A+1)) | grep -v "NA")
	err=$(echo ${match1} ${match2} | tr ' ' '\n' | sort | uniq -u | grep -v "\." | wc -l)
	compar=$((compar+1))
	tmperr=$((tmperr+err))
	err=0
	fi
	done
	echo $tmperr genotypic error in $compar comparisons.
	echo -e "$locName\t$NallSeq\t$NallSize\t$Nmiss\t$tmperr\t$compar\t$(echo "scale=5; ($tmperr/($compar*$A))*100" | bc 2> /dev/null)% " >> LocusCharacteristics_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt )
done	
rm tmpTable.txt



######################## analysis with PARAMSET2 ###################################

# running stuttermark 
echo "***************FDSTools stuttermark ${PARAMSET2} analysis***************"
rm -f -r stuttermark; mkdir stuttermark
for file in ./tssvResults/*.txt
do
echo ${file:14:-4}
fdstools stuttermark -s=-1:${S2a},+1:${S2b} -m 2 -l ${ANALYSIS}_${TYPE}_FDSTools_InputFile.txt $file ./stuttermark/${file:14:-4}_stuttermarked.txt
done

# running allelefinder 
echo "***************FDSTools allelefinder ${PARAMSET2} analysis***************"
rm -f -r allelefinder; mkdir allelefinder
rm -f -r alellefinderReports; mkdir alellefinderReports
for file in ./stuttermark/*.txt
do
echo ${file:14:-4}
fdstools allelefinder -m ${M2} -n ${N2} -a ${A} -x ${X} -M ${MNP} -c annotation -l ${ANALYSIS}_${TYPE}_FDSTools_InputFile.txt $file -o ./allelefinder/${file:14:-4}_allelefinder.txt -R ./alellefinderReports/${file:14:-4}_AFreport.txt
done


# concatenate all resulting allele across individuals (except the firt line which is a header)
echo "***************Creating a unique genotype file***************"
echo -e "sample\tmarker\treads\tallele" > AllAlleles.txt
awk 'FNR>1' ./allelefinder/*.txt >> AllAlleles.txt

#### convertir les séquence des allèles en code allele standards et formater les genotypes

# creating allele naming convention file
echo "***************Listing allele sequence, allele size, allele code and allele count***************"
rm -f AlleleNamingConvention_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt
while read locus
do
	echo $locus
	nall=$(awk -v l=${locus} '$2==l { print $4 }' AllAlleles.txt | sort | uniq | wc -l)
	echo $nall
	allmax=$(expr 100 + $nall)
	echo $allmax
	all=$(awk -v l=${locus} '$2==l { print $4 }' AllAlleles.txt | sort --version-sort | uniq)	
	paste <(yes $locus | head -n $nall) <(echo "${all}") <(eval echo {101..$allmax} | tr " " "\n") <(for alleleSeq in $(awk -v l=${locus} '$2==l { print $4 }' AllAlleles.txt | sort --version-sort | uniq); do grep -w $locus AllAlleles.txt | cut -f4 | grep -x $alleleSeq | wc -l; done) >> AlleleNamingConvention_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt
done < ${ANALYSIS}_LocusList.txt
wait
sed -i '/^\s\s*/d' AlleleNamingConvention_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt

# renaming the alleles
echo "***************Allele coding***************"
rm -f AlleleNamed.txt
while read locus
do
	echo $locus
	grep ^$locus AlleleNamingConvention_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt | cut -f2,3 | sed -e $'s!\t!$/!' | sed -e 's$^$s/^$' | sed -e 's?$?/g?' > tmp_seddScript.txt
	paste <(awk -v l=$locus '$2==l { print $1, $2 }' AllAlleles.txt) <(awk -v l=$locus '$2==l { print $4 }' AllAlleles.txt | sed -f tmp_seddScript.txt) | sed 's/\t/ /' >> AlleleNamed.txt
done < ${ANALYSIS}_LocusList.txt
rm tmp_seddScript.txt

# clearing the names of individuals
echo "***************Individual name cleaning***************"
sed -i -E 's%./stuttermark/%%g' AlleleNamed.txt
sed -i -E 's%_tssv_stuttermarked%%g' AlleleNamed.txt
sed -i -E 's%.assembled%%g' AlleleNamed.txt

# formating allele information file
echo "***************Retrieving allele information across locus***************"
rm -f AllAlleleSeq.txt
while read locus seq code n
do
counter=1
rm -f tmp_alleleSeq.txt
stg=$(echo $seq | sed "s/(/\t/g" | sed "s/)/\t/g" | sed "s/ /\t/g")
while [ $counter -le $(echo $stg | awk '{print NF}' | sort -nu | tail -n 1) ]
do
	seq=$(echo $stg | cut -d ' ' -f$counter)
	c2=$(expr $counter + 1)
	rep=$(echo $stg | cut -d ' ' -f$c2)
	yes $seq | head -n $rep | tr "\n" " " | sed "s/ //g" | sed -e '$a\' >> tmp_alleleSeq.txt
	((counter++))
	((counter++))
done
paste <(cat tmp_alleleSeq.txt | tr "\n" " " | sed "s/ //g" | sed -e '$a\') <(echo $(cat tmp_alleleSeq.txt | tr "\n" " " | sed "s/ //g" | sed -e '$a\' | wc -m)-1 | bc) >> AllAlleleSeq.txt
done < AlleleNamingConvention_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt

echo -e "Locus\tAlleleSequenceAnnotated\tAlleleSeqCode\tOccurancesAcrossIndivs\tAlleleSequence\tAlleleLength" > AlleleInformationFile_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt
paste <(cat AlleleNamingConvention_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt) <(cat AllAlleleSeq.txt) >> AlleleInformationFile_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt

rm AlleleNamingConvention_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt

# formating the genotypic table across all individuals / locus
echo "***************Genotypic table formating***************"
rm -f tmp.txt tmp_indiv.txt
paste <(echo "Indiv") <(cut -f1 ${ANALYSIS}_LocusList.txt | tr "\n" "\t" | sed "s/\t/$(for i in $(seq 1 ${A}); do printf "\t"; done)/g") > GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt
for file in ./samples/*.fast[a-q]
do
	ind=${file:10:-6}
	echo $ind > tmp_indiv.txt
	echo $ind
	while read locus
	do
	grep -w "^${ind}" AlleleNamed.txt | grep -w "$locus" > tmp.txt
	lines=$(wc -l tmp.txt | cut -f1 -d" ")
	if [ $lines -gt 0 ]
	then
	alleles=$(echo -e $(cat tmp.txt | cut -f3 -d" " | sort | tr "\n" "\t") $(for i in $(seq 1 $((A-lines))); do printf "\t."; done)) 
	elif [ $lines = 0 ]
	then	
	alleles=$(echo -e "NA" $(for i in $(seq 1 $((A-1))); do printf "\tNA"; done))
	else
	alleles=$(echo -e "Err" $(for i in $(seq 1 $((A-1))); do printf "\tErr"; done))
	fi
	echo $alleles >> tmp_indiv.txt
	done < ${ANALYSIS}_LocusList.txt
	cat tmp_indiv.txt | tr "\n" "\t" | sed -e '$a\' >> GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt
done
sed -i "s/ /\t/g" GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt
sed -i "s/_S[0-9]*\t/\t/g" GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt
rm tmp.txt tmp_indiv.txt


### Summarizing locus characteristics and comparing repeated genotypes to compute genotypic error rate 
echo "***************Comparing repeated genotypes to computing allelic error rate, and retrieving the number of alleles and missing data rate across the whole genotypic table ***************"
echo -e "Loci\tNallelesSequence\tNalleleSize\tMissingRate\tAllelicMismatches\tPolyploidGenotypesComparred\tGenotypicError" > LocusCharacteristics_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt
for loc in $(eval echo "{1..$(wc -l < ${ANALYSIS}_LocusList.txt)}")
do
	echo $loc
	locName=$(head GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt -n 1 | sed "s/$(for i in $(seq 1 ${A}); do printf "\t"; done)/\t/g" | cut -f$(($loc + 1)))
	echo $locName
	## Number of alleles based on sequence
	NallSeq=$(grep -w $locName AlleleInformationFile_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt | cut -f3 | uniq | wc -l)
	## Number of alleles based on size
	NallSize=$(grep -w $locName AlleleInformationFile_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt | cut -f6 | sort -g |  uniq |wc -l)
	## Missing data
	cl1=$(echo "($loc-1)*${A}+2" | bc )
	cl2=$(echo "$cl1+(${A}-1)" | bc )
	Nmiss=$(echo "scale=5; ($(sed 1d GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt | cut -f$cl1 | grep "NA" | wc -l) / $(sed 1d GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt | cut -f$cl1 | wc -l))*100" | bc -l 2> /dev/null)%
	## Error rates
	tmperr=0
	compar=0
	paste <(sed 1d GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt | cut -f1) <(sed 1d GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt | cut -f$cl1-$cl2) > tmpTable.txt
	cat tmpTable.txt | cut -f1 | sed "s/_[A-B]$//g" | uniq -d | (while read ind
	do
	#echo $ind	
	na=$(grep ${ind}_[A-B] tmpTable.txt | cut -f2-$(($A+1)) | grep -v "NA" | wc -l)
	if [[ -n $na ]] && [[ $na > 1 ]]
	then
	match1=$(grep ${ind}_A tmpTable.txt | cut -f2-$(($A+1)) | grep -v "NA")
	match2=$(grep ${ind}_B tmpTable.txt | cut -f2-$(($A+1)) | grep -v "NA")
	err=$(echo ${match1} ${match2} | tr ' ' '\n' | sort | uniq -u | grep -v "\." | wc -l)
	compar=$((compar+1))
	tmperr=$((tmperr+err))
	err=0
	fi
	done
	echo $tmperr genotypic error in $compar comparisons.
	echo -e "$locName\t$NallSeq\t$NallSize\t$Nmiss\t$tmperr\t$compar\t$(echo "scale=5; ($tmperr/($compar*$A))*100" | bc 2> /dev/null)% " >> LocusCharacteristics_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt )
done	
rm tmpTable.txt




###############################################################################################
# ANALYSIS OF DIPLOID nSSR : 

elif [ $ANALYSIS = "nSSR" ] && [ $A -eq 2 ]
	then

# tssv
echo "***************FDSTools tssv analysis***************"
rm -f -r tssvResults; mkdir tssvResults
rm -f -r tssvReports; mkdir tssvReports
# for file in ./samples/*.fast[a-q]
for file in ./samples/*.fast*
do
echo ${file:10:-6}
fdstools tssv --num-threads 4 -R ./tssvReports/${file:10:-6}_tssvReport.txt -m 0.08 ${ANALYSIS}_${TYPE}_FDSTools_InputFile.txt $file ./tssvResults/${file:10:-6}_tssv.txt
done
 
# extracting the number of reads by locus by individuals
echo "***************Computing coverage by locus by individual***************"
rm -f LocusCoverageperIndividual_${ANALYSIS}_${TYPE}.txt
paste <(echo "Locus") <(ls tssvReports | sed "s/_tssvReport.txt//g" | tr "\n" "\t" | sed -e '$a\') > LocusCoverageperIndividual_${ANALYSIS}_${TYPE}.txt
while read locus
do
	echo $locus
	echo $locus > tmpLocusCoverage.txt
	for file in ./tssvReports/*
	do
	#echo $file	
	cat $file | grep -w $locus | cut -f3 >> tmpLocusCoverage.txt
	done
cat tmpLocusCoverage.txt | tr "\n" "\t" | sed -e '$a\' >> LocusCoverageperIndividual_${ANALYSIS}_${TYPE}.txt
done < ${ANALYSIS}_LocusList.txt
rm tmpLocusCoverage.txt

######################### analysis with PARAMSET1###################################

# running stuttermark 
echo "***************FDSTools stuttermark ${PARAMSET} ${S} analysis***************"
rm -f -r stuttermark; mkdir stuttermark
for file in ./tssvResults/*.txt
do
echo ${file:14:-4}
fdstools stuttermark -s=-1:${S1a},+1:${S1b} -m 2 -l ${ANALYSIS}_${TYPE}_FDSTools_InputFile.txt $file ./stuttermark/${file:14:-4}_stuttermarked.txt
done

# running allelefinder 
echo "***************FDSTools allelefinder ${PARAMSET} ${M} ${N} analysis***************"
rm -f -r allelefinder; mkdir allelefinder
rm -f -r alellefinderReports; mkdir alellefinderReports
for file in ./stuttermark/*.txt
do
echo ${file:14:-4}
fdstools allelefinder -m ${M1} -n ${N1} -x ${X} -M ${MNP} -c annotation -l ${ANALYSIS}_${TYPE}_FDSTools_InputFile.txt $file -o ./allelefinder/${file:14:-4}_allelefinder.txt -R ./alellefinderReports/${file:14:-4}_AFreport.txt
done

# concatenate all resulting allele across individuals (except the firt line which is a header)
echo "***************Creating a unique genotype file***************"
echo -e "sample\tmarker\treads\tallele" > AllAlleles.txt
awk 'FNR>1' ./allelefinder/*.txt >> AllAlleles.txt

#### convertir les séquence des allèles en code allele standards et formater les genotypes

# creating allele naming convention file
echo "***************Listing allele sequence, allele size, allele code and allele count***************"
rm -f AlleleNamingConvention_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt
while read locus
do
	echo $locus
	nall=$(awk -v l=${locus} '$2==l { print $4 }' AllAlleles.txt | sort | uniq | wc -l)
	echo $nall
	allmax=$(expr 100 + $nall)
	echo $allmax
	all=$(awk -v l=${locus} '$2==l { print $4 }' AllAlleles.txt | sort --version-sort | uniq)	
	paste <(yes $locus | head -n $nall) <(echo "${all}") <(eval echo {101..$allmax} | tr " " "\n") <(for alleleSeq in $(awk -v l=${locus} '$2==l { print $4 }' AllAlleles.txt | sort --version-sort | uniq); do grep -w $locus AllAlleles.txt | cut -f4 | grep -x $alleleSeq | wc -l; done) >> AlleleNamingConvention_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt
done < ${ANALYSIS}_LocusList.txt
wait
sed -i '/^\s\s*/d' AlleleNamingConvention_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt

# renaming the alleles
echo "***************Allele coding***************"
rm -f AlleleNamed.txt
while read locus
do
	echo $locus
	grep ^$locus AlleleNamingConvention_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt | cut -f2,3 | sed -e $'s!\t!$/!' | sed -e 's$^$s/^$' | sed -e 's?$?/g?' > tmp_seddScript.txt
	paste <(awk -v l=$locus '$2==l { print $1, $2 }' AllAlleles.txt) <(awk -v l=$locus '$2==l { print $4 }' AllAlleles.txt | sed -f tmp_seddScript.txt) | sed 's/\t/ /' >> AlleleNamed.txt
done < ${ANALYSIS}_LocusList.txt
rm tmp_seddScript.txt

# clearing the names of individuals
echo "***************Individual name cleaning***************"
sed -i -E 's%./stuttermark/%%g' AlleleNamed.txt
sed -i -E 's%_tssv_stuttermarked%%g' AlleleNamed.txt
sed -i -E 's%.assembled%%g' AlleleNamed.txt

# formating allele information file
echo "***************Retrieving allele information across locus***************"
rm -f AllAlleleSeq.txt
while read locus seq code n
do
counter=1
rm -f tmp_alleleSeq.txt
stg=$(echo $seq | sed "s/(/\t/g" | sed "s/)/\t/g" | sed "s/ /\t/g")
while [ $counter -le $(echo $stg | awk '{print NF}' | sort -nu | tail -n 1) ]
do
	seq=$(echo $stg | cut -d ' ' -f$counter)
	c2=$(expr $counter + 1)
	rep=$(echo $stg | cut -d ' ' -f$c2)
	yes $seq | head -n $rep | tr "\n" " " | sed "s/ //g" | sed -e '$a\' >> tmp_alleleSeq.txt
	((counter++))
	((counter++))
done
paste <(cat tmp_alleleSeq.txt | tr "\n" " " | sed "s/ //g" | sed -e '$a\') <(echo $(cat tmp_alleleSeq.txt | tr "\n" " " | sed "s/ //g" | sed -e '$a\' | wc -m)-1 | bc) >> AllAlleleSeq.txt
done < AlleleNamingConvention_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt

echo -e "Locus\tAlleleSequenceAnnotated\tAlleleSeqCode\tOccurancesAcrossIndivs\tAlleleSequence\tAlleleLength" > AlleleInformationFile_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt
paste <(cat AlleleNamingConvention_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt) <(cat AllAlleleSeq.txt) >> AlleleInformationFile_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt

rm AlleleNamingConvention_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt

# formating the genotypic table across all individuals / locus
echo "***************Genotypic table formating***************"
rm -f tmp.txt tmp_indiv.txt
paste <(echo "Indiv") <(cut -f1 ${ANALYSIS}_LocusList.txt | tr "\n" "\t" | sed "s/\t/\t\t/g") > GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt
for file in ./samples/*.fast[a-q]
do
	ind=${file:10:-6}
	echo $ind > tmp_indiv.txt
	echo $ind
	while read locus
	do
	grep -w "^${ind}" AlleleNamed.txt | grep -w "$locus" > tmp.txt
	lines=$(wc -l tmp.txt | cut -f1 -d" ")
	if [ $lines = 2 ]
	then
	alleles=$(cat tmp.txt | cut -f3 -d" " | sort | tr "\n" "\t")
	elif [ $lines = 1 ]
	then
	alleles=$(paste <(cat tmp.txt | cut -f3 -d" ") <(cat tmp.txt | cut -f3 -d" "))
	elif [ $lines = 0 ]
	then	
	alleles=$(echo -e "NA\tNA")
	else
	alleles=$(echo -e "Err\tErr")
	fi
	echo $alleles >> tmp_indiv.txt
	done < ${ANALYSIS}_LocusList.txt
	cat tmp_indiv.txt | tr "\n" "\t" | sed -e '$a\' >> GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt
done
sed -i "s/ /\t/g" GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt
sed -i "s/_S[0-9]*\t/\t/g" GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt
rm tmp.txt tmp_indiv.txt

### Summarizing locus characteristics and comparing repeated genotypes to compute genotypic error rate 
echo "***************Comparing repeated genotypes to computing allelic error rate, and retrieving the number of alleles and missing data rate across the whole genotypic table ***************"
echo -e "Loci\tNallelesSequence\tNalleleSize\tMissingRate\tAllelicMismatches\tDiploidGenotypesComparred\tAllelicError" > LocusCharacteristics_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt
for loc in $(eval echo "{1..$(wc -l < ${ANALYSIS}_LocusList.txt)}")
do
	echo $loc
	locName=$(head GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt -n 1 | sed "s/\t\t/\t/g" | cut -f$(($loc + 1)))
	echo $locName
	## Number of alleles based on sequence
	NallSeq=$(grep -w $locName AlleleInformationFile_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt | cut -f3 | uniq | wc -l)
	## Number of alleles based on size
	NallSize=$(grep -w $locName AlleleInformationFile_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt | cut -f6 | sort -g |  uniq |wc -l)
	## Missing data
	cl1=$(expr $loc \* 2)
	cl2=$(expr $cl1 + 1)
	Nmiss=$(echo "scale=5; ($(sed 1d GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt | cut -f$cl1 | grep "NA" | wc -l) / $(sed 1d GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt | cut -f$cl1 | wc -l))*100" | bc -l 2> /dev/null)%
	## Error rates
	tmperr=0
	compar=0
	paste <(sed 1d GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt | cut -f1) <(sed 1d GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt | cut -f$cl1,$cl2) > tmpTable.txt
	cat tmpTable.txt | cut -f1 | sed "s/_[A-B]$//g" | uniq -d | (while read ind
	do
	#echo $ind	
	na=$(grep ${ind}_[A-B] tmpTable.txt | cut -f2,3 | grep -v "NA" | wc -l)
	if [[ -n $na ]] && [[ $na > 1 ]]
	then
	i1al1=$(grep ${ind}_[A-B] tmpTable.txt | cut -f2 | head -1)
	i1al2=$(grep ${ind}_[A-B] tmpTable.txt | cut -f3 | head -1)
	i2al1=$(grep ${ind}_[A-B] tmpTable.txt | cut -f2 | tail -1)
	i2al2=$(grep ${ind}_[A-B] tmpTable.txt | cut -f3 | tail -1)
	if [[ $i1al1 == $i2al1 ]] && [[ $i1al2 == $i2al2 ]]
	then
	err=0
	compar=$((compar+1))
	#echo case 1.1
	elif [[ $i1al1 != $i2al1 || $i1al1 != $i2al2 ]] && [[ $i1al2 == $i2al1 || $i1al2 == $i2al2 ]]
	then
	err=1
	compar=$((compar+1))
	#echo case 1.2
	elif [[ $i1al1 == $i2al1 || $i1al1 == $i2al2 ]] && [[ $i1al2 != $i2al1 || $i1al2 != $i2al2 ]]
	then err=1
	compar=$((compar+1))
	#echo case 1.3
	elif [[ $i1al1 != $i2al1 && $i1al2 != $i2al2 ]] && [[ $i1al1 != $i2al2 && $i1al2 != $i2al1 ]]
	then
	err=2
	compar=$((compar+1))
	#echo case 1.4
	else echo "not possible case - error"
	fi
	tmperr=$((tmperr+err))
	err=0
	fi
	done
	echo $tmperr allelic error in $compar comparisons.
	echo -e "$locName\t$NallSeq\t$NallSize\t$Nmiss\t$tmperr\t$compar\t$(echo "scale=5; ($tmperr/(2*$compar))*100" | bc 2> /dev/null)% " >> LocusCharacteristics_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt )
done	
rm tmpTable.txt

######################### analysis with PARAMSET2 ###################################

# running stuttermark 
echo "***************stuttermark ${PARAMSET} ${S} analysis***************"
rm -f -r stuttermark; mkdir stuttermark
for file in ./tssvResults/*.txt
do
echo ${file:14:-4}
fdstools stuttermark -s=-1:${S2a},+1:${S2b} -m 2 -l ${ANALYSIS}_${TYPE}_FDSTools_InputFile.txt $file ./stuttermark/${file:14:-4}_stuttermarked.txt
done

# running allelefinder 
echo "***************allelefinder ${PARAMSET} ${M} ${N} analysis***************"
rm -f -r allelefinder; mkdir allelefinder
rm -f -r alellefinderReports; mkdir alellefinderReports
for file in ./stuttermark/*.txt
do
echo ${file:14:-4}
fdstools allelefinder -m ${M2} -n ${N2} -x ${X} -M ${MNP} -c annotation -l ${ANALYSIS}_${TYPE}_FDSTools_InputFile.txt $file -o ./allelefinder/${file:14:-4}_allelefinder.txt -R ./alellefinderReports/${file:14:-4}_AFreport.txt
done

# concatenate all resulting allele across individuals (except the firt line which is a header)
echo "***************creating a single genotype file***************"
echo -e "sample\tmarker\treads\tallele" > AllAlleles.txt
awk 'FNR>1' ./allelefinder/*.txt >> AllAlleles.txt

#### convertir les séquence des allèles en code allele standards et formater les genotypes

# creating allele naming convention file
echo "***************listing allele sequence, allele code and allele count***************"
rm -f AlleleNamingConvention_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt
while read locus
do
	echo $locus
	nall=$(awk -v l=${locus} '$2==l { print $4 }' AllAlleles.txt | sort | uniq | wc -l)
	echo $nall
	allmax=$(expr 100 + $nall)
	echo $allmax
	all=$(awk -v l=${locus} '$2==l { print $4 }' AllAlleles.txt | sort --version-sort | uniq)	
	paste <(yes $locus | head -n $nall) <(echo "${all}") <(eval echo {101..$allmax} | tr " " "\n") <(for alleleSeq in $(awk -v l=${locus} '$2==l { print $4 }' AllAlleles.txt | sort --version-sort | uniq); do grep -w $locus AllAlleles.txt | cut -f4 | grep -x $alleleSeq | wc -l; done) >> AlleleNamingConvention_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt
done < ${ANALYSIS}_LocusList.txt
wait
sed -i '/^\s\s*/d' AlleleNamingConvention_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt

# renaming the alleles
echo "***************allele coding***************"
rm -f AlleleNamed.txt
while read locus
do
	echo $locus
	grep ^$locus AlleleNamingConvention_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt | cut -f2,3 | sed -e $'s!\t!$/!' | sed -e 's$^$s/^$' | sed -e 's?$?/g?' > tmp_seddScript.txt
	paste <(awk -v l=$locus '$2==l { print $1, $2 }' AllAlleles.txt) <(awk -v l=$locus '$2==l { print $4 }' AllAlleles.txt | sed -f tmp_seddScript.txt) | sed 's/\t/ /' >> AlleleNamed.txt
done < ${ANALYSIS}_LocusList.txt
rm tmp_seddScript.txt

# clearing the names of individuals
echo "***************individual name cleaning***************"
sed -i -E 's%./stuttermark/%%g' AlleleNamed.txt
sed -i -E 's%_tssv_stuttermarked%%g' AlleleNamed.txt
sed -i -E 's%.assembled%%g' AlleleNamed.txt

# formating allele information file
echo "***************Retrieving allele information across locus***************"
rm -f AllAlleleSeq.txt
while read locus seq code n
do
counter=1
rm -f tmp_alleleSeq.txt
stg=$(echo $seq | sed "s/(/\t/g" | sed "s/)/\t/g" | sed "s/ /\t/g")
while [ $counter -le $(echo $stg | awk '{print NF}' | sort -nu | tail -n 1) ]
do
	seq=$(echo $stg | cut -d ' ' -f$counter)
	c2=$(expr $counter + 1)
	rep=$(echo $stg | cut -d ' ' -f$c2)
	yes $seq | head -n $rep | tr "\n" " " | sed "s/ //g" | sed -e '$a\' >> tmp_alleleSeq.txt
	((counter++))
	((counter++))
done
paste <(cat tmp_alleleSeq.txt | tr "\n" " " | sed "s/ //g" | sed -e '$a\') <(echo $(cat tmp_alleleSeq.txt | tr "\n" " " | sed "s/ //g" | sed -e '$a\' | wc -m)-1 | bc) >> AllAlleleSeq.txt
done < AlleleNamingConvention_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt

echo -e "Locus\tAlleleSequenceAnnotated\tAlleleSeqCode\tOccurancesAcrossIndivs\tAlleleSequence\tAlleleLength" > AlleleInformationFile_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt
paste <(cat AlleleNamingConvention_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt) <(cat AllAlleleSeq.txt) >> AlleleInformationFile_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt

rm AlleleNamingConvention_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt

# formating the genotypic table across all individuals / locus
echo "***************Genotypic table formating***************"
rm -f tmp.txt tmp_indiv.txt
paste <(echo "Indiv") <(cut -f1 ${ANALYSIS}_LocusList.txt | tr "\n" "\t" | sed "s/\t/\t\t/g") > GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt
for file in ./samples/*.fast[a-q]
do
	ind=${file:10:-6}
	echo $ind > tmp_indiv.txt
	echo $ind
	while read locus
	do
	grep -w "^${ind}" AlleleNamed.txt | grep -w "$locus" > tmp.txt
	lines=$(wc -l tmp.txt | cut -f1 -d" ")
	if [ $lines = 2 ]
	then
	alleles=$(cat tmp.txt | cut -f3 -d" " | sort | tr "\n" "\t")
	elif [ $lines = 1 ]
	then
	alleles=$(paste <(cat tmp.txt | cut -f3 -d" ") <(cat tmp.txt | cut -f3 -d" "))
	elif [ $lines = 0 ]
	then	
	alleles=$(echo -e "NA\tNA")
	else
	alleles=$(echo -e "Err\tErr")
	fi
	echo $alleles >> tmp_indiv.txt
	done < ${ANALYSIS}_LocusList.txt
	cat tmp_indiv.txt | tr "\n" "\t" | sed -e '$a\' >> GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt
done
sed -i "s/ /\t/g" GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt
sed -i "s/_S[0-9]*\t/\t/g" GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt
rm tmp.txt tmp_indiv.txt

### Summarizing locus characteristics and comparing repeated genotypes to compute genotypic error rate 
echo "***************Comparing repeated genotypes to computing allelic error rate, and retrieving the number of alleles and missing data rate across the whole genotypic table ***************"
echo -e "Loci\tNallelesSequence\tNalleleSize\tMissingRate\tAllelicMismatches\tDiploidGenotypesComparred\tAllelicError" > LocusCharacteristics_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt
for loc in $(eval echo "{1..$(wc -l < ${ANALYSIS}_LocusList.txt)}")
do
	echo $loc
	locName=$(head GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt -n 1 | sed "s/\t\t/\t/g" | cut -f$(($loc + 1)))
	echo $locName
	## Number of alleles based on sequence
	NallSeq=$(grep -w $locName AlleleInformationFile_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt | cut -f3 | uniq | wc -l)
	## Number of alleles based on size
	NallSize=$(grep -w $locName AlleleInformationFile_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt | cut -f6 | sort -g | uniq |wc -l)
	## Missing data
	cl1=$(expr $loc \* 2)
	cl2=$(expr $cl1 + 1)
	Nmiss=$(echo "scale=5; ($(sed 1d GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt | cut -f$cl1 | grep "NA" | wc -l) / $(sed 1d GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt | cut -f$cl1 | wc -l))*100" | bc -l 2> /dev/null)%
	## Error rates
	tmperr=0
	compar=0
	paste <(sed 1d GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt | cut -f1) <(sed 1d GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt | cut -f$cl1,$cl2) > tmpTable.txt
	cat tmpTable.txt | cut -f1 | sed "s/_[A-B]$//g" | uniq -d | (while read ind
	do
	#echo $ind	
	na=$(grep ${ind}_[A-B] tmpTable.txt | cut -f2,3 | grep -v "NA" | wc -l)
	if [[ -n $na ]] && [[ $na > 1 ]]
	then
	i1al1=$(grep ${ind}_[A-B] tmpTable.txt | cut -f2 | head -1)
	i1al2=$(grep ${ind}_[A-B] tmpTable.txt | cut -f3 | head -1)
	i2al1=$(grep ${ind}_[A-B] tmpTable.txt | cut -f2 | tail -1)
	i2al2=$(grep ${ind}_[A-B] tmpTable.txt | cut -f3 | tail -1)
	if [[ $i1al1 == $i2al1 ]] && [[ $i1al2 == $i2al2 ]]
	then
	err=0
	compar=$((compar+1))
	#echo case 1.1
	elif [[ $i1al1 != $i2al1 || $i1al1 != $i2al2 ]] && [[ $i1al2 == $i2al1 || $i1al2 == $i2al2 ]]
	then
	err=1
	compar=$((compar+1))
	#echo case 1.2
	elif [[ $i1al1 == $i2al1 || $i1al1 == $i2al2 ]] && [[ $i1al2 != $i2al1 || $i1al2 != $i2al2 ]]
	then err=1
	compar=$((compar+1))
	#echo case 1.3
	elif [[ $i1al1 != $i2al1 && $i1al2 != $i2al2 ]] && [[ $i1al1 != $i2al2 && $i1al2 != $i2al1 ]]
	then
	err=2
	compar=$((compar+1))
	#echo case 1.4
	else echo "not possible case - error"
	fi
	tmperr=$((tmperr+err))
	err=0
	fi
	done
	echo $tmperr allelic error in $compar comparisons.
	echo -e "$locName\t$NallSeq\t$NallSize\t$Nmiss\t$tmperr\t$compar\t$(echo "scale=5; ($tmperr/(2*$compar))*100" | bc 2> /dev/null)% " >> LocusCharacteristics_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt )
done	
rm tmpTable.txt

rm AllAlleleSeq.txt tmp_alleleSeq.txt AllAlleles.txt AlleleNamed.txt


###############################################################################################
# Specific section for cpSSR analysis

elif [ $ANALYSIS = "cpSSR" ] # cpSSR analysis

then

# tssv
echo "***************FDSTools tssv analysis***************"
rm -f -r tssvResults; mkdir tssvResults
rm -f -r tssvReports; mkdir tssvReports
for file in ./samples/*.fast[a-q]
do
echo ${file:10:-6}
fdstools tssv --num-threads 4 -R ./tssvReports/${file:10:-6}_tssvReport.txt -m 0.08 ${ANALYSIS}_${TYPE}_FDSTools_InputFile.txt $file ./tssvResults/${file:10:-6}_tssv.txt
done
 
# extracting the number of reads by locus by individuals
echo "***************Computing coverage by locus by individual***************"
rm -f LocusCoverageperIndividual_${ANALYSIS}_${TYPE}.txt
paste <(echo "Locus") <(ls tssvReports | sed "s/_tssvReport.txt//g" | tr "\n" "\t" | sed -e '$a\') > LocusCoverageperIndividual_${ANALYSIS}_${TYPE}.txt
while read locus
do
	echo $locus
	echo $locus > tmpLocusCoverage.txt
	for file in ./tssvReports/*
	do
	#echo $file	
	cat $file | grep -w $locus | cut -f3 >> tmpLocusCoverage.txt
	done
cat tmpLocusCoverage.txt | tr "\n" "\t" | sed -e '$a\' >> LocusCoverageperIndividual_${ANALYSIS}_${TYPE}.txt
done < ${ANALYSIS}_LocusList.txt
rm tmpLocusCoverage.txt

######################### analysis with PARAMSET1###################################

# running stuttermark 
echo "***************FDSTools stuttermark ${PARAMSET} ${S} analysis***************"
rm -f -r stuttermark; mkdir stuttermark
for file in ./tssvResults/*.txt
do
echo ${file:14:-4}
fdstools stuttermark -s=-1:${S1a},+1:${S1b} -m 2 -l ${ANALYSIS}_${TYPE}_FDSTools_InputFile.txt $file ./stuttermark/${file:14:-4}_stuttermarked.txt
done

# running allelefinder 
echo "***************FDSTools allelefinder ${PARAMSET} ${M} ${N} analysis***************"
rm -f -r allelefinder; mkdir allelefinder
rm -f -r alellefinderReports; mkdir alellefinderReports
for file in ./stuttermark/*.txt
do
echo ${file:14:-4}
fdstools allelefinder -m ${M1} -n ${N1} -x 10 -M 75 -c annotation -l ${ANALYSIS}_${TYPE}_FDSTools_InputFile.txt $file -o ./allelefinder/${file:14:-4}_allelefinder.txt -R ./alellefinderReports/${file:14:-4}_AFreport.txt
done

# specific cpSSR calling section : when more than one allele called : select the one with the highest number of repeats

rm -f -r allelefinder_cpSSRfiltered; mkdir allelefinder_cpSSRfiltered
for file in ./allelefinder/*.txt
do
echo ${file:15:-4}
rm tmp_indiv.txt
while read locus
do
	rm tmp.txt
	grep -w "$locus" $file > tmp.txt
	lines=$(wc -l tmp.txt | cut -f1 -d" ")
	if [ $lines = 2 ]
	then
	num1=$(cat tmp.txt | head -1 | awk -vRS=")" -vFS="(" '{print $2}' | head -n -1 | tr '\n' '+')
	res1=$(echo $((${num1::-1})))
	num2=$(cat tmp.txt | tail -1 | awk -vRS=")" -vFS="(" '{print $2}' | head -n -1 | tr '\n' '+')
	res2=$(echo $((${num2::-1})))
	echo $num1
	echo $res1
	echo $num2
	echo $res2
	paste <(cat tmp.txt | head -1) <(echo $res1) > tmp_all.txt
	paste <(cat tmp.txt | tail -1) <(echo $res2) >> tmp_all.txt
	cat tmp_all.txt | sort -n -k 5 -r
	allele=$(cat tmp_all.txt | sort -n -k 5 -r | head -1 | awk 'NF{NF-=1};1')
	else
	allele=$(cat tmp.txt)
	fi
	echo $allele
	echo $allele >> tmp_indiv.txt
done < ${ANALYSIS}_LocusList.txt
echo -e "sample\tmarker\ttotal\tallele" > ./allelefinder_cpSSRfiltered/${file:15}
cat tmp_indiv.txt | sed -e 's/ /\t/g' >> ./allelefinder_cpSSRfiltered/${file:15}
done


# concatenate all resulting allele across individuals (except the firt line which is a header)
echo "***************Creating a unique genotype file***************"
echo -e "sample\tmarker\treads\tallele" > AllAlleles.txt
awk 'FNR>1' ./allelefinder_cpSSRfiltered/*.txt >> AllAlleles.txt


#### convertir les séquence des allèles en code allele standards et formater les genotypes

# creating allele naming convention file
echo "***************Listing allele sequence, allele size, allele code and allele count***************"
rm -f AlleleNamingConvention_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt
while read locus
do
	echo $locus
	nall=$(awk -v l=${locus} '$2==l { print $4 }' AllAlleles.txt | sort | uniq | wc -l)
	echo $nall
	allmax=$(expr 100 + $nall)
	echo $allmax
	all=$(awk -v l=${locus} '$2==l { print $4 }' AllAlleles.txt | sort --version-sort | uniq)	
	paste <(yes $locus | head -n $nall) <(echo "${all}") <(eval echo {101..$allmax} | tr " " "\n") <(for alleleSeq in $(awk -v l=${locus} '$2==l { print $4 }' AllAlleles.txt | sort --version-sort | uniq); do grep -w $locus AllAlleles.txt | cut -f4 | grep -x $alleleSeq | wc -l; done) >> AlleleNamingConvention_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt
done < ${ANALYSIS}_LocusList.txt
wait
sed -i '/^\s\s*/d' AlleleNamingConvention_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt

# renaming the alleles
echo "***************Allele coding***************"
rm -f AlleleNamed.txt
while read locus
do
	echo $locus
	grep ^$locus AlleleNamingConvention_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt | cut -f2,3 | sed -e $'s!\t!$/!' | sed -e 's$^$s/^$' | sed -e 's?$?/g?' > tmp_seddScript.txt
	paste <(awk -v l=$locus '$2==l { print $1, $2 }' AllAlleles.txt) <(awk -v l=$locus '$2==l { print $4 }' AllAlleles.txt | sed -f tmp_seddScript.txt) | sed 's/\t/ /' >> AlleleNamed.txt
done < ${ANALYSIS}_LocusList.txt
rm tmp_seddScript.txt

# clearing the names of individuals
echo "***************Individual name cleaning***************"
sed -i -E 's%./stuttermark/%%g' AlleleNamed.txt
sed -i -E 's%_tssv_stuttermarked%%g' AlleleNamed.txt
sed -i -E 's%.assembled%%g' AlleleNamed.txt

# formating allele information file
echo "***************Retrieving allele information across locus***************"
rm -f AllAlleleSeq.txt
while read locus seq code n
do
counter=1
rm -f tmp_alleleSeq.txt
stg=$(echo $seq | sed "s/(/\t/g" | sed "s/)/\t/g" | sed "s/ /\t/g")
while [ $counter -le $(echo $stg | awk '{print NF}' | sort -nu | tail -n 1) ]
do
	seq=$(echo $stg | cut -d ' ' -f$counter)
	c2=$(expr $counter + 1)
	rep=$(echo $stg | cut -d ' ' -f$c2)
	yes $seq | head -n $rep | tr "\n" " " | sed "s/ //g" | sed -e '$a\' >> tmp_alleleSeq.txt
	((counter++))
	((counter++))
done
paste <(cat tmp_alleleSeq.txt | tr "\n" " " | sed "s/ //g" | sed -e '$a\') <(echo $(cat tmp_alleleSeq.txt | tr "\n" " " | sed "s/ //g" | sed -e '$a\' | wc -m)-1 | bc) >> AllAlleleSeq.txt
done < AlleleNamingConvention_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt

echo -e "Locus\tAlleleSequenceAnnotated\tAlleleSeqCode\tOccurancesAcrossIndivs\tAlleleSequence\tAlleleLength" > AlleleInformationFile_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt
paste <(cat AlleleNamingConvention_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt) <(cat AllAlleleSeq.txt) >> AlleleInformationFile_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt

rm AlleleNamingConvention_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt

# formating the genotypic table across all individuals / locus
echo "***************Genotypic table formating***************"
rm -f tmp.txt tmp_indiv.txt
paste <(echo "Indiv") <(cut -f1 ${ANALYSIS}_LocusList.txt | tr "\n" "\t" | sed "s/\t/\t\t/g") > GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt
for file in ./samples/*.fast[a-q]
do
	ind=${file:10:-6}
	echo $ind > tmp_indiv.txt
	echo $ind
	while read locus
	do
	grep -w "^${ind}" AlleleNamed.txt | grep -w "$locus" > tmp.txt
	lines=$(wc -l tmp.txt | cut -f1 -d" ")
	if [ $lines = 2 ]
	then
	alleles=$(cat tmp.txt | cut -f3 -d" " | sort | tr "\n" "\t")
	elif [ $lines = 1 ]
	then
	alleles=$(paste <(cat tmp.txt | cut -f3 -d" ") <(cat tmp.txt | cut -f3 -d" "))
	elif [ $lines = 0 ]
	then	
	alleles=$(echo -e "NA\tNA")
	else
	alleles=$(echo -e "Err\tErr")
	fi
	echo $alleles >> tmp_indiv.txt
	done < ${ANALYSIS}_LocusList.txt
	cat tmp_indiv.txt | tr "\n" "\t" | sed -e '$a\' >> GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt
done
sed -i "s/ /\t/g" GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt
sed -i "s/_S[0-9]*\t/\t/g" GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt
rm tmp.txt tmp_indiv.txt

### Summarizing locus characteristics and comparing repeated genotypes to compute genotypic error rate 
echo "***************Comparing repeated genotypes to computing allelic error rate, and retrieving the number of alleles and missing data rate across the whole genotypic table ***************"
echo -e "Loci\tNallelesSequence\tNalleleSize\tMissingRate\tAllelicMismatches\tHaploidGenotypesComparred\tAllelicError" > LocusCharacteristics_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt
for loc in $(eval echo "{1..$(wc -l < ${ANALYSIS}_LocusList.txt)}")
do
	echo $loc
	locName=$(head GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt -n 1 | sed "s/\t\t/\t/g" | cut -f$(($loc + 1)))
	echo $locName
	## Number of alleles based on sequence
	NallSeq=$(grep -w $locName AlleleInformationFile_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt | cut -f3 | uniq | wc -l)
	## Number of alleles based on size
	NallSize=$(grep -w $locName AlleleInformationFile_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt | cut -f6 | sort -g |  uniq |wc -l)
	## Missing data
	cl1=$(expr $loc \* 2)
	cl2=$(expr $cl1 + 1)
	Nmiss=$(echo "scale=5; ($(sed 1d GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt | cut -f$cl1 | grep "NA" | wc -l) / $(sed 1d GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt | cut -f$cl1 | wc -l))*100" | bc -l 2> /dev/null)%
	## Error rates
	tmperr=0
	compar=0
	paste <(sed 1d GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt | cut -f1) <(sed 1d GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt | cut -f$cl1,$cl2) > tmpTable.txt
	cat tmpTable.txt | cut -f1 | sed "s/_[A-B]$//g" | uniq -d | (while read ind
	do
	#echo $ind	
	na=$(grep ${ind}_[A-B] tmpTable.txt | cut -f2,3 | grep -v "NA" | wc -l)
	if [[ -n $na ]] && [[ $na > 1 ]]
	then
	i1al1=$(grep ${ind}_[A-B] tmpTable.txt | cut -f2 | head -1)
	i1al2=$(grep ${ind}_[A-B] tmpTable.txt | cut -f3 | head -1)
	i2al1=$(grep ${ind}_[A-B] tmpTable.txt | cut -f2 | tail -1)
	i2al2=$(grep ${ind}_[A-B] tmpTable.txt | cut -f3 | tail -1)
	if [[ $i1al1 == $i2al1 ]] && [[ $i1al2 == $i2al2 ]]
	then
	err=0
	compar=$((compar+1))
	#echo case 1.1
	elif [[ $i1al1 != $i2al1 || $i1al1 != $i2al2 ]] && [[ $i1al2 == $i2al1 || $i1al2 == $i2al2 ]]
	then
	err=1
	compar=$((compar+1))
	#echo case 1.2
	elif [[ $i1al1 == $i2al1 || $i1al1 == $i2al2 ]] && [[ $i1al2 != $i2al1 || $i1al2 != $i2al2 ]]
	then err=1
	compar=$((compar+1))
	#echo case 1.3
	elif [[ $i1al1 != $i2al1 && $i1al2 != $i2al2 ]] && [[ $i1al1 != $i2al2 && $i1al2 != $i2al1 ]]
	then
	err=2
	compar=$((compar+1))
	#echo case 1.4
	else echo "not possible case - error"
	fi
	tmperr=$((tmperr+err))
	err=0
	fi
	done
	echo $tmperr allelic error in $compar comparisons.
	echo -e "$locName\t$NallSeq\t$NallSize\t$Nmiss\t$tmperr\t$compar\t$(echo "scale=5; ($tmperr/($compar))*100" | bc 2> /dev/null)% " >> LocusCharacteristics_${ANALYSIS}_${TYPE}_${PARAMSET1}_sa${S1a}_sb${S1b}_m${M1}_n${N1}.txt )
done	
rm tmpTable.txt

######################### analysis with PARAMSET2 ###################################

# running stuttermark 
echo "***************stuttermark ${PARAMSET} ${S} analysis***************"
rm -f -r stuttermark; mkdir stuttermark
for file in ./tssvResults/*.txt
do
echo ${file:14:-4}
fdstools stuttermark -s=-1:${S2a},+1:${S2b} -m 2 -l ${ANALYSIS}_${TYPE}_FDSTools_InputFile.txt $file ./stuttermark/${file:14:-4}_stuttermarked.txt
done

# running allelefinder 
echo "***************allelefinder ${PARAMSET} ${M} ${N} analysis***************"
rm -f -r allelefinder; mkdir allelefinder
rm -f -r alellefinderReports; mkdir alellefinderReports
for file in ./stuttermark/*.txt
do
echo ${file:14:-4}
fdstools allelefinder -m ${M2} -n ${N2} -x 10 -M 75 -c annotation -l ${ANALYSIS}_${TYPE}_FDSTools_InputFile.txt $file -o ./allelefinder/${file:14:-4}_allelefinder.txt -R ./alellefinderReports/${file:14:-4}_AFreport.txt
done

# specific cpSSR calling section : when more than one allele called : select the one with the highest number of repeats

rm -f -r allelefinder_cpSSRfiltered; mkdir allelefinder_cpSSRfiltered
for file in ./allelefinder/*.txt
do
echo ${file:15:-4}
rm tmp_indiv.txt
while read locus
do
	rm tmp.txt
	grep -w "$locus" $file > tmp.txt
	lines=$(wc -l tmp.txt | cut -f1 -d" ")
	if [ $lines = 2 ]
	then
	num1=$(cat tmp.txt | head -1 | awk -vRS=")" -vFS="(" '{print $2}' | head -n -1 | tr '\n' '+')
	res1=$(echo $((${num1::-1})))
	num2=$(cat tmp.txt | tail -1 | awk -vRS=")" -vFS="(" '{print $2}' | head -n -1 | tr '\n' '+')
	res2=$(echo $((${num2::-1})))
	echo $num1
	echo $res1
	echo $num2
	echo $res2
	paste <(cat tmp.txt | head -1) <(echo $res1) > tmp_all.txt
	paste <(cat tmp.txt | tail -1) <(echo $res2) >> tmp_all.txt
	cat tmp_all.txt | sort -n -k 5 -r
	allele=$(cat tmp_all.txt | sort -n -k 5 -r | head -1 | awk 'NF{NF-=1};1')
	else
	allele=$(cat tmp.txt)
	fi
	echo $allele
	echo $allele >> tmp_indiv.txt
done < ${ANALYSIS}_LocusList.txt
echo -e "sample\tmarker\ttotal\tallele" > ./allelefinder_cpSSRfiltered/${file:15}
cat tmp_indiv.txt | sed -e 's/ /\t/g' >> ./allelefinder_cpSSRfiltered/${file:15}
done


# concatenate all resulting allele across individuals (except the firt line which is a header)
echo "***************Creating a unique genotype file***************"
echo -e "sample\tmarker\treads\tallele" > AllAlleles.txt
awk 'FNR>1' ./allelefinder_cpSSRfiltered/*.txt >> AllAlleles.txt

#### convertir les séquence des allèles en code allele standards et formater les genotypes

# creating allele naming convention file
echo "***************listing allele sequence, allele code and allele count***************"
rm -f AlleleNamingConvention_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt
while read locus
do
	echo $locus
	nall=$(awk -v l=${locus} '$2==l { print $4 }' AllAlleles.txt | sort | uniq | wc -l)
	echo $nall
	allmax=$(expr 100 + $nall)
	echo $allmax
	all=$(awk -v l=${locus} '$2==l { print $4 }' AllAlleles.txt | sort --version-sort | uniq)	
	paste <(yes $locus | head -n $nall) <(echo "${all}") <(eval echo {101..$allmax} | tr " " "\n") <(for alleleSeq in $(awk -v l=${locus} '$2==l { print $4 }' AllAlleles.txt | sort --version-sort | uniq); do grep -w $locus AllAlleles.txt | cut -f4 | grep -x $alleleSeq | wc -l; done) >> AlleleNamingConvention_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt
done < ${ANALYSIS}_LocusList.txt
wait
sed -i '/^\s\s*/d' AlleleNamingConvention_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt

# renaming the alleles
echo "***************allele coding***************"
rm -f AlleleNamed.txt
while read locus
do
	echo $locus
	grep ^$locus AlleleNamingConvention_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt | cut -f2,3 | sed -e $'s!\t!$/!' | sed -e 's$^$s/^$' | sed -e 's?$?/g?' > tmp_seddScript.txt
	paste <(awk -v l=$locus '$2==l { print $1, $2 }' AllAlleles.txt) <(awk -v l=$locus '$2==l { print $4 }' AllAlleles.txt | sed -f tmp_seddScript.txt) | sed 's/\t/ /' >> AlleleNamed.txt
done < ${ANALYSIS}_LocusList.txt
rm tmp_seddScript.txt

# clearing the names of individuals
echo "***************individual name cleaning***************"
sed -i -E 's%./stuttermark/%%g' AlleleNamed.txt
sed -i -E 's%_tssv_stuttermarked%%g' AlleleNamed.txt
sed -i -E 's%.assembled%%g' AlleleNamed.txt

# formating allele information file
echo "***************Retrieving allele information across locus***************"
rm -f AllAlleleSeq.txt
while read locus seq code n
do
counter=1
rm -f tmp_alleleSeq.txt
stg=$(echo $seq | sed "s/(/\t/g" | sed "s/)/\t/g" | sed "s/ /\t/g")
while [ $counter -le $(echo $stg | awk '{print NF}' | sort -nu | tail -n 1) ]
do
	seq=$(echo $stg | cut -d ' ' -f$counter)
	c2=$(expr $counter + 1)
	rep=$(echo $stg | cut -d ' ' -f$c2)
	yes $seq | head -n $rep | tr "\n" " " | sed "s/ //g" | sed -e '$a\' >> tmp_alleleSeq.txt
	((counter++))
	((counter++))
done
paste <(cat tmp_alleleSeq.txt | tr "\n" " " | sed "s/ //g" | sed -e '$a\') <(echo $(cat tmp_alleleSeq.txt | tr "\n" " " | sed "s/ //g" | sed -e '$a\' | wc -m)-1 | bc) >> AllAlleleSeq.txt
done < AlleleNamingConvention_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt

echo -e "Locus\tAlleleSequenceAnnotated\tAlleleSeqCode\tOccurancesAcrossIndivs\tAlleleSequence\tAlleleLength" > AlleleInformationFile_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt
paste <(cat AlleleNamingConvention_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt) <(cat AllAlleleSeq.txt) >> AlleleInformationFile_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt

rm AlleleNamingConvention_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt

# formating the genotypic table across all individuals / locus
echo "***************Genotypic table formating***************"
rm -f tmp.txt tmp_indiv.txt
paste <(echo "Indiv") <(cut -f1 ${ANALYSIS}_LocusList.txt | tr "\n" "\t" | sed "s/\t/\t\t/g") > GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt
for file in ./samples/*.fast[a-q]
do
	ind=${file:10:-6}
	echo $ind > tmp_indiv.txt
	echo $ind
	while read locus
	do
	grep -w "^${ind}" AlleleNamed.txt | grep -w "$locus" > tmp.txt
	lines=$(wc -l tmp.txt | cut -f1 -d" ")
	if [ $lines = 2 ]
	then
	alleles=$(cat tmp.txt | cut -f3 -d" " | sort | tr "\n" "\t")
	elif [ $lines = 1 ]
	then
	alleles=$(paste <(cat tmp.txt | cut -f3 -d" ") <(cat tmp.txt | cut -f3 -d" "))
	elif [ $lines = 0 ]
	then	
	alleles=$(echo -e "NA\tNA")
	else
	alleles=$(echo -e "Err\tErr")
	fi
	echo $alleles >> tmp_indiv.txt
	done < ${ANALYSIS}_LocusList.txt
	cat tmp_indiv.txt | tr "\n" "\t" | sed -e '$a\' >> GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt
done
sed -i "s/ /\t/g" GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt
sed -i "s/_S[0-9]*\t/\t/g" GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt
rm tmp.txt tmp_indiv.txt

### Summarizing locus characteristics and comparing repeated genotypes to compute genotypic error rate 
echo "***************Comparing repeated genotypes to computing allelic error rate, and retrieving the number of alleles and missing data rate across the whole genotypic table ***************"
echo -e "Loci\tNallelesSequence\tNalleleSize\tMissingRate\tAllelicMismatches\tHaploidGenotypesComparred\tAllelicError" > LocusCharacteristics_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt
for loc in $(eval echo "{1..$(wc -l < ${ANALYSIS}_LocusList.txt)}")
do
	echo $loc
	locName=$(head GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt -n 1 | sed "s/\t\t/\t/g" | cut -f$(($loc + 1)))
	echo $locName
	## Number of alleles based on sequence
	NallSeq=$(grep -w $locName AlleleInformationFile_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt | cut -f3 | uniq | wc -l)
	## Number of alleles based on size
	NallSize=$(grep -w $locName AlleleInformationFile_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt | cut -f6 | sort -g | uniq |wc -l)
	## Missing data
	cl1=$(expr $loc \* 2)
	cl2=$(expr $cl1 + 1)
	Nmiss=$(echo "scale=5; ($(sed 1d GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt | cut -f$cl1 | grep "NA" | wc -l) / $(sed 1d GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt | cut -f$cl1 | wc -l))*100" | bc -l 2> /dev/null)%
	## Error rates
	tmperr=0
	compar=0
	paste <(sed 1d GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt | cut -f1) <(sed 1d GenotypicTable_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt | cut -f$cl1,$cl2) > tmpTable.txt
	cat tmpTable.txt | cut -f1 | sed "s/_[A-B]$//g" | uniq -d | (while read ind
	do
	#echo $ind	
	na=$(grep ${ind}_[A-B] tmpTable.txt | cut -f2,3 | grep -v "NA" | wc -l)
	if [[ -n $na ]] && [[ $na > 1 ]]
	then
	i1al1=$(grep ${ind}_[A-B] tmpTable.txt | cut -f2 | head -1)
	i1al2=$(grep ${ind}_[A-B] tmpTable.txt | cut -f3 | head -1)
	i2al1=$(grep ${ind}_[A-B] tmpTable.txt | cut -f2 | tail -1)
	i2al2=$(grep ${ind}_[A-B] tmpTable.txt | cut -f3 | tail -1)
	if [[ $i1al1 == $i2al1 ]] && [[ $i1al2 == $i2al2 ]]
	then
	err=0
	compar=$((compar+1))
	#echo case 1.1
	elif [[ $i1al1 != $i2al1 || $i1al1 != $i2al2 ]] && [[ $i1al2 == $i2al1 || $i1al2 == $i2al2 ]]
	then
	err=1
	compar=$((compar+1))
	#echo case 1.2
	elif [[ $i1al1 == $i2al1 || $i1al1 == $i2al2 ]] && [[ $i1al2 != $i2al1 || $i1al2 != $i2al2 ]]
	then err=1
	compar=$((compar+1))
	#echo case 1.3
	elif [[ $i1al1 != $i2al1 && $i1al2 != $i2al2 ]] && [[ $i1al1 != $i2al2 && $i1al2 != $i2al1 ]]
	then
	err=2
	compar=$((compar+1))
	#echo case 1.4
	else echo "not possible case - error"
	fi
	tmperr=$((tmperr+err))
	err=0
	fi
	done
	echo $tmperr allelic error in $compar comparisons.
	echo -e "$locName\t$NallSeq\t$NallSize\t$Nmiss\t$tmperr\t$compar\t$(echo "scale=5; ($tmperr/($compar))*100" | bc 2> /dev/null)% " >> LocusCharacteristics_${ANALYSIS}_${TYPE}_${PARAMSET2}_sa${S2a}_sb${S2b}_m${M2}_n${N2}.txt )
done	
rm tmpTable.txt

rm AllAlleleSeq.txt tmp_alleleSeq.txt AllAlleles.txt AlleleNamed.txt

else

echo "Analysis type "nSSR or cpSSR" not properly specified. Please try again."

fi
