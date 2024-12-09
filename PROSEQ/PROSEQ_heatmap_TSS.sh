#!/bin/bash

#check args
if [ $# -ne 3 ]; then
    echo "Error: Invalid number of arguments!"
    usage
fi

# Description: Pipeline for PRO-seq plotheatmap using deeptools
# bwGroupinfo.txt (n* n)
# control1 treat1 * * * *
# control2 treat2 * * * *

# bedinfo.txt (n * 2)
# bedG1_fwd.bed bedG1_rev.bed
# bedG2_fwd.bed bedG2_rev.bed

# bwType: fulllength / singlebase

## Path
workDir=$1
trackDir=$workDir/03_spikein_scaled_bw_track
heatmapDir=$workDir/Profile_Heatmap
sampleGroupinfo=${heatmapDir}/bwGroupinfo.txt
bedGroupDir=$2
bedGroupinfo=${bedGroupDir}/bedinfo.txt
bwType=$3

#plot arguments 
upNum=1000
downNum=1000
binsize=20
refPointLable="TSS"
## Compute matrix
#addtreat="sort_byL"
#rmpart1="-dTAG-DLD1_231122"
#rmpart2="hg19.ucsc.refseq_proteincoding_fwd_sense_"
cat $sampleGroupinfo | while read group;
do
    arr=($group)
    control1=${arr[0]}
    fwdbwFiles=()
    revbwFiles=()
    for element in "${arr[@]}"; 
    do
        fwdbwFiles+=("${trackDir}/${element}_${bwType}_fwd.bw")
        revbwFiles+=("${trackDir}/${element}_${bwType}_rev.bw")
    done
    fwdbwArgs="${fwdbwFiles[@]}"
    revbwArgs="${revbwFiles[@]}"


	cat $bedGroupinfo | while read beds;
	do
	    echo "==================== Samples running: ========================"
	    echo "computeMatrix $group in $beds" 
	    echo "=============================================================="
	    arr=($beds)
	    fwd=${arr[0]}
	    bedFwd=${bedGroupDir}/${arr[0]}
	    bedRev=${bedGroupDir}/${arr[1]}
	    outPrefix=${control1}-${fwd}
	#    outPrefix=${control1/$rmpart1/}-${fwd#"$rmpart2"}
	    if [ ! -f $heatmapDir/${outPrefix}_TSS_fwd-rev_sense.gz ]
	    then
	    # forward_sense
		computeMatrix  reference-point --referencePoint TSS -p 24  -b ${upNum} -a ${downNum} \
		    -R $bedFwd \
		    -S ${fwdbwArgs} \
		    --binSize ${binsize} --missingDataAsZero --skipZeros \
		    -o $heatmapDir/${outPrefix}_TSS_fwd_sense.gz

	    # reverse_sense
		computeMatrix  reference-point --referencePoint TSS  -p 24  -b ${upNum} -a ${downNum} \
		    -R $bedRev \
		    -S ${revbwArgs} \
		    --binSize ${binsize} --missingDataAsZero --skipZeros \
		    -o $heatmapDir/${outPrefix}_TSS_rev_sense.gz

	    # merge
		computeMatrixOperations rbind \
		    -m $heatmapDir/${outPrefix}_TSS_fwd_sense.gz $heatmapDir/${outPrefix}_TSS_rev_sense.gz \
		    -o $heatmapDir/${outPrefix}_TSS_fwd-rev_sense.gz
	    fi
	done
done


## Plot
cat $sampleGroupinfo | while read group;
do
 
    arr1=($group)
    control1=${arr1[0]}
   # outPrefix=${control/$rmpart1/}-${addtreat}
	cat $bedGroupinfo | while read beds;
	do
	    echo "==================== Samples running: ========================"
	    echo "plot $group in  $beds" 
	    echo "=============================================================="
	    arr=($beds)
	    fwd=${arr[0]}
	    outPrefix=${control1}-${fwd}
	#    outPrefix=${control1/$rmpart1/}-${bed#"$rmpart2"}
	    matrix=$heatmapDir/${outPrefix}_TSS_fwd-rev_sense.gz

	    # Heatmap
	    plotHeatmap -m $matrix \
		-out $heatmapDir/${outPrefix}_TSS_Heatmap.pdf \
		--missingDataColor "white" \
		--samplesLabel ${arr1[@]} \
		--sortUsingSamples 1 \
		--refPointLabel ${refPointLable} \
		--heatmapHeight 12 \
		--plotFileFormat pdf  --dpi 720
	    # Profile
	    plotProfile -m $matrix \
		--perGroup \
		--samplesLabel ${arr1[@]} \
		--refPointLabel ${refPointLable} \
		-out $heatmapDir/${outPrefix}_TSS_Profile.pdf \
		--plotHeight 16  --plotWidth 20 

	done

done
