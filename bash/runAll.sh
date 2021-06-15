#!/bin/bash

#JETR,INRVAL

DATE=`date +%Y%m%d`

mkdir -p logs/$DATE/

rs=(4 2)
files=(ntupleToHist_PbPbMC)

dPhi=("pi/2" "2pi/3" "3pi/4" "5pi/6")
multiJtDPhi=("pi/2" "2pi/3" "3pi/4" "5pi/6" "7pi/8")

#mixJtDR=(0.4 0.8 1.0)
mixJtDR=(0.8)

#ptMins=(25 30 36)
#nPtBins=(16 15 14)

#genMinPt=(15 25)
genMinPt=(15)

nomNPtBins=16
customBinsBase=25,30,36,44,53,64,77,93,112,135,163,196,236,285,344,415,500
#customBinsBase=40,44,53,64,77,93,112,135,163,196,236,285,344,415,500


#Overriding the nominal for quick tests
files=(ntupleToHist_PPData ntupleToHist_PPMC ntupleToHist_PbPbData ntupleToHist_PbPbMC)
rs=(4)
#ptMins=(30 36)
ptMins=(30 36)
dPhi=("pi/2")
#multiJtDPhi=("pi/2" "7pi/8")
multiJtDPhi=("pi/2" "7pi/8")

for i in ${rs[@]}
do
    echo "R=0.$i"

    for j in ${dPhi[@]}
    do
	dPhiName=$(echo $j | sed -e "s@/@Over@g")
	echo " dPhi=$j, $dPhiName"
	for k in ${multiJtDPhi[@]}
	do
	    multiJtDPhiName=$(echo $k | sed -e "s@/@Over@g")
	    echo "  multiJtDPhi=$k, $multiJtDPhiName"

	    for n in ${mixJtDR[@]}
	    do
		mixJtDRName=$(echo $n | sed -e "s@\.@p@g")
	    
		echo "  mixJtDR=$n, $mixJtDRName"

		for p in ${genMinPt[@]}
		do
		    
		    for l in ${ptMins[@]}
		    do
			nPtBinsTemp=$nomNPtBins
			ptBins=$customBinsBase
			
			startVal=${ptBins#*,}
			startVal=${ptBins%,$startVal}


			while [[ $startVal -ne $l ]]
			do
			    ptBins=${ptBins#*,}
			    nPtBinsTemp=$((nPtBinsTemp - 1))

#			    startVal=${ptBins%,*}
			    startVal=${ptBins#*,}
			    startVal=${ptBins%,$startVal}

			done

			echo "   ptMin=$l, nPtBins=$nPtBinsTemp, ptBins=$ptBins"
			pos=$((pos + 1))
			
			for m in ${files[@]}
			do
			    newFile=input/ntupleToHist/"$m"_R"$i"_"$dPhiName"_"$multiJtDPhiName"_DR"$mixJtDRName"_GenMin"$p"_"$l".config
			    logFile=logs/$DATE/"$m"_R"$i"_"$dPhiName"_"$multiJtDPhiName"_DR"$mixJtDRName"_GenMin"$p"_"$l".log
			    
			    cp input/ntupleToHist/$m.config $newFile
			    sed -i "s@INNJTPTBINS@$nPtBinsTemp@g" $newFile
			    sed -i "s@INRVAL@$i@g" $newFile
			    sed -i "s@INGENMINPT@$p@g" $newFile
			    sed -i "s@INJTPTBINSLOW@$l@g" $newFile
			    sed -i "s@INJTPTBINSCUSTOM@$ptBins@g" $newFile
			    sed -i "s@INGAMMAJTDPHINAME@$dPhiName@g" $newFile
			    sed -i "s@INGAMMAMULTIJTDPHINAME@$multiJtDPhiName@g" $newFile
			    sed -i "s@INGAMMAJTDPHI@$j@g" $newFile
			    sed -i "s@INGAMMAMULTIJTDPHI@$k@g" $newFile
			    sed -i "s@INMIXJETEXCLUSIONDRNAME@$mixJtDRName@g" $newFile
			    sed -i "s@INMIXJETEXCLUSIONDR@$n@g" $newFile
			
			    ./bin/gdjNTupleToHist.exe $newFile &> $logFile &
#			    exit 1
			done		   
		    done
		    wait
		done
	    done
	done
    done
done

exit 1

#files=(ntupleToHist_PbPbData ntupleToHist_PPData)
#files=(ntupleToHist_PPMC ntupleToHist_PPData)
#files=(ntupleToHist_PPData)

for i in ${rs[@]}
do
    for j in ${files[@]}
    do
	cp input/ntupleToHist/$j.config  input/ntupleToHist/"$j"_R$i.config 
	sed -i -e "s@INRVAL@$i@g" input/ntupleToHist/"$j"_R$i.config

	echo "./bin/gdjNTupleToHist.exe input/ntupleToHist/"$j"_R$i.config &> logs/$DATE/"$j"_R"$i"_$DATE.log &"
    done    
    wait
done

echo "RUNALL.SH COMPLETE!"
