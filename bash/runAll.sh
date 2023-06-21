#!/bin/bash

#JETR,INRVAL

DATE=`date +%Y%m%d`

mkdir -p logs/$DATE/

rs=(2 4)
files=(ntupleToHist_PPMC)

dPhi=("pi/2" "2pi/3" "3pi/4" "5pi/6")
multiJtDPhi=("pi/2" "2pi/3" "3pi/4" "5pi/6" "7pi/8")

#mixJtDR=(0.4 0.8 1.0)

#ptMins=(25 30 36)
#nPtBins=(16 15 14)

#genMinPt=(15 25)
genMinPt=(20)

nomNPtBins=26
#customBinsBase=25,30,36,44,53,64,77,93,112,135,163,196,236,285,344,415,500
#customBinsBase=40,44,53,64,77,93,112,135,163,196,236,285,344,415,500
#YJ Bins
customBinsBase=15,20,25,30,35,40,45,50,56,63,70,79,89,100,112,125,141,158,177,199,223,251,281,316,354,398,501
nomNSubPtBins=2
customSubBinsBase=15,30,501

#Overriding the nominal for quick tests
#files=(ntupleToHist_PPData)
#files=(ntupleToHist_PbPbMC ntupleToHist_PbPbData ntupleToHist_PPMC ntupleToHist_PPData)
files=(ntupleToHist_PPMC)
rs=(2)
#ptMins=(30 36)
#ptMinsR2=(30)
#ptMinsR4=(36)
ptMinsR2=(15)
ptMinsR2Reco=(30)
ptMinsR2RecoSyst=(35)

ptMinsR4=(35)
mixJtDRR2=(0.4)
mixJtDRR4=(0.8)

dPhi=("pi/2")
#multiJtDPhi=("pi/2" "7pi/8")
multiJtDPhi=("7pi/8")
#multiJtDPhi=("pi/2")

for i in ${rs[@]}
do
    echo "R=0.$i"

    ptMins=()
    ptMinsReco=()
    ptMinsRecoSyst=()

    if [[ $i -eq 2 ]]
    then
	pos=0
	for tempI in  ${ptMinsR2[@]}
	do
 	    ptMins+=($tempI)
	    ptMinsReco+=(${ptMinsR2Reco[$pos]})
	    ptMinsRecoSyst+=(${ptMinsR2RecoSyst[$pos]})
	    pos=$((pos+1))
	done
    elif [[ $i -eq 4 ]]
    then
	for tempI in  ${ptMinsR4[@]}
	do
	    ptMins+=($tempI)
	done	
    fi

    mixJtDR=()
    if [[ $i -eq 2 ]]
    then
	for tempI in  ${mixJtDRR2[@]}
	do
	    mixJtDR+=($tempI)
	done
    elif [[ $i -eq 4 ]]
    then
	for tempI in  ${mixJtDRR4[@]}
	do
	    mixJtDR+=($tempI)
	done	
    fi


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

		drVal=$n
		if [[ $i == 4 ]] 
		then
		    drVal=0.8
		fi

		mixJtDRName=$(echo $drVal | sed -e "s@\.@p@g")
	    
		echo "  mixJtDR=$drVal, $mixJtDRName"

		for p in ${genMinPt[@]}
		do

		    pos=0
		    for l in ${ptMins[@]}
		    do
			nPtBinsTemp=$nomNPtBins
			nSubPtBinsTemp=$nomNSubPtBins

			ptBins=$customBinsBase
			subPtBins=$customSubBinsBase
			
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

			if [[ $startVal -ne $l ]]
			then
			    echo "Requested ptMin='$l' not found in custom ptbins array $ptBins. return"
			    exit 1
			fi

			startVal=${subPtBins#*,}
			startVal=${subPtBins%,$startVal}
 
			while [[ $startVal -ne $l ]]
			do
			    subPtBins=${subPtBins#*,}
			    nSubPtBinsTemp=$((nSubPtBinsTemp - 1))

#			    startVal=${ptBins%,*}
			    startVal=${subPtBins#*,}
			    startVal=${subPtBins%,$startVal}
			done

			if [[ $startVal -ne $l ]]
			then
			    echo "Requested ptMin='$l' not found in custom subPtbins array $subPtBins. return"
			    exit 1
			fi

			ptMax=$ptBins		     
			while [[ $ptMax == *","* ]]
			do
			    ptMax=${ptMax#*,}
			done

			ptMaxReco=${ptBins%,$ptMax}
			while [[ $ptMaxReco == *","* ]]
			do
			    ptMaxReco=${ptMaxReco#*,}
			done

			ptMinReco=${ptMinsReco[$pos]}
			ptMinRecoSyst=${ptMinsRecoSyst[$pos]}


			echo "   ptMin=$l, ptMax=$ptMax, nPtBins=$nPtBinsTemp, ptBins=$ptBins"
			pos=$((pos + 1))
#			exit 1
			
			for m in ${files[@]}
			do
			    newFile=input/ntupleToHist/"$m"_R"$i"_"$dPhiName"_"$multiJtDPhiName"_DR"$mixJtDRName"_GenMin"$p"_"$l".config
			    logFile=logs/$DATE/"$m"_R"$i"_"$dPhiName"_"$multiJtDPhiName"_DR"$mixJtDRName"_GenMin"$p"_"$l".log
			    
			    cp input/ntupleToHist/$m.config $newFile
			    sed -i "s@INNJTPTBINS@$nPtBinsTemp@g" $newFile
			    sed -i "s@INNSUBJTPTBINS@$nSubPtBinsTemp@g" $newFile
			    sed -i "s@INRVAL@$i@g" $newFile
			    sed -i "s@INGENMINPT@$p@g" $newFile
			    sed -i "s@INJTPTBINSLOWRECOSYST@$ptMinRecoSyst@g" $newFile
			    sed -i "s@INJTPTBINSLOWRECO@$ptMinReco@g" $newFile
			    sed -i "s@INJTPTBINSHIGHRECO@$ptMaxReco@g" $newFile
			    sed -i "s@INJTPTBINSLOW@$l@g" $newFile
			    sed -i "s@INJTPTBINSHIGH@$ptMax@g" $newFile
			    sed -i "s@INJTPTBINSCUSTOM@$ptBins@g" $newFile
			    sed -i "s@INSUBJTPTBINSCUSTOM@$subPtBins@g" $newFile
			    sed -i "s@INGAMMAJTDPHINAME@$dPhiName@g" $newFile
			    sed -i "s@INGAMMAMULTIJTDPHINAME@$multiJtDPhiName@g" $newFile
			    sed -i "s@INGAMMAJTDPHI@$j@g" $newFile
			    sed -i "s@INGAMMAMULTIJTDPHI@$k@g" $newFile
			    sed -i "s@INMIXJETEXCLUSIONDRNAME@$mixJtDRName@g" $newFile
			    sed -i "s@INMIXJETEXCLUSIONDR@$n@g" $newFile
			
			    echo "./bin/gdjNTupleToHist.exe $newFile &> $logFile &"
			    exit 1
			done
			wait
		    done
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
