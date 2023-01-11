#!/usr/bin/sh

help() {
 
    echo "please do as follow:"
    # echo "注意 -r -u -f -c 参数请勿同时使用 "
    #echo "    -h          : help information"
    echo "    Userage : identifier.sh transcript.fasta model thread"
    echo "    transcript.fasta : the full-length transcript in fasta format"
    echo "    model : choosing from [arabidopsis, human], arabidopsis or rice for plant, human for animal"
    echo "    thread : the number of thread to construct the colored de Bruijn graph(cDBG)"
    echo "      "
    echo "    All the output files located in the path of input fasta file which prefix is the name of the fasta file"
    exit 1
}
 
if [[ $# == 0 || "$1" == "-h"  || "$1" == "-help"  || "$1" == "--help" ]]; then
    help
    exit 1
fi
if [[ $# > 3 ]]; then
    echo "We need three parameters, 1st for transcript.fasta, 2th for predicting model, 3rd for thread"
    exit 1
fi
filename=$1
if [[ $# == 2 ]]; then
    if [[ ${filename:0-2} != "fa" && ${filename:0-5} != "fasta" ]]; then
        echo "We need two parameters, 1st for transcript.fasta"
        exit 1
    fi
fi
if [[ $# == 2 && $2 != "arabidopsis" && $2 != "human" ]]; then
    echo "We need two parameters, 2th for predicting model, choosing from [arabidopsis, rice, human], arabidopsis or rice for plant, human for animal"
    exit 1
fi

echo "Step1: make blast database"
date
makeblastdb -in $1 -dbtype nucl

echo "Step2: sequence alignment using blastn"
date
blastn -query $1 -db $1 -strand plus -evalue 1E-10 -outfmt 6 -ungapped -num_threads 20 -out $1\_blastout.txt




echo "Step3: predict AS transcript pair"
rm $1*split* 
python3 unique.py $1\_blastout.txt $3
date

for i in `ls $1*split*`
do
    #echo $i
    python3 cdbg.py $i $1 > $i\_four_AS.seq &
done

touch $1\_done

echo ""
echo "Start building the cDBG and update the progress every minute"
t=0
while true
do
    s=`cat $1*done|wc -l`
    u=`cat $1*.unique*|wc -l`
    echo "It takes "$t" minutes. Completed:" $s/$u = `awk 'BEGIN{printf "%.1f%%\n",('$s'/'$u')*100}'`
    if [ $u == $s ]; then
        echo "done"
        break
    fi
    t=$(($t+1))
    sleep 60
done

cat $1*split*AF >$1\_AF.txt
cat $1*split*AL >$1\_AL.txt
cat $1*split*MX >$1\_MX.txt
cat $1*split*snp >$1\_snp
cat $1*split*bubble >$1\_bubble
cat $1*split*_four_AS.seq >$1\_four_AS.seq
rm $1*split*

echo "Step4: get features"
date

python getfeature.py $1\_four_AS.seq $1\_four_AS.txt > $1\_four_AS.feature


echo "Step5: classify AS transcript pair"
date

python predict.py $1\_four_AS.feature $2 $1\_four_AS.txt $1\_four_AS_type.txt

echo "done, thank u"
date
