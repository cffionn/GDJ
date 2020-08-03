#!/bin/bash

DATE=`date +%Y%m%d`

mkdir -p logs/$DATE

inFile=input/ntupleToMBHist.config
outFile=${inFile%.config}
jetR=(2 4)

for i in ${jetR[@]}
do    
    cp $inFile "$outFile"_R$i.config
    sed -i -e "s@INRVAL@$i@g" "$outFile"_R$i.config

    ./bin/gdjNTupleToMBHist.exe "$outFile"_R$i.config >& logs/$DATE/ntupleToMBHist_R"$i"_"$DATE".log &
done

wait
echo "runNTupleToMBHist.sh Complete!"
