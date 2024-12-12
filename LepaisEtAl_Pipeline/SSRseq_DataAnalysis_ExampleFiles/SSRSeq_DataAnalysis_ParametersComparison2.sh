
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




for file in ./stuttermark/*.txt
do
echo ${file:14:-4}
# echo "***fdstools allelefinder -m ${M1} -n ${N1} -a ${A} -x ${X} -M ${MNP} -c annotation***"
echo "fdstools allelefinder -m ${M1} -n ${N1} -a ${A} -x ${X} -M ${MNP} -c annotation -l ${ANALYSIS}_${TYPE}_FDSTools_InputFile.txt $file -o ./allelefinder/${file:14:-4}_allelefinder.txt -R ./alellefinderReports/${file:14:-4}_AFreport.txt"
done