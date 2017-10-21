#!/bin/bash

## ENTER VARIABLES HERE
genomeWindows=../../genomeWindows/hg19_100kbBins_partition_excludingGaps_100kbBinsOnly_noChrYnoChrM.bed
FILES=""

## Make subdirectories to put new files in
mkdir bedGraph
mkdir counts
mkdir features
mkdir features/bedGraph
mkdir features/counts

## Make bedgraphs of peaks
numLinesWindowFile=`cat $genomeWindows | wc -l`

for file in $FILES
do
name=`stringEdit $file .bed .bedGraph | stringEdit - ../`
echo Working with $file >> bedGraphloop.log
echo Naming it $name >> bedGraphloop.log
coverageBed -counts -a $file -b $genomeWindows | sortBed -i - > bedGraph/$name
numLinesFile=`cat $file | wc -l`
numLinesCounts=`cat bedGraph/$name | wc -l`
echo There are $numLinesFile lines in the original file >> bedGraphloop.log
echo ...and $numLinesCounts lines in the counts file >> bedGraphloop.log
echo ...and $numLinesWindowFile in the genomic window file >> bedGraphloop.log
echo >> bedGraphloop.log
done

### counts
for file in bedGraph/*
do
name=`stringEdit $file .bedGraph .txt | stringEdit - bedGraph/`
echo Working with $file >> bedGraphloop.log
echo Naming it $name >> bedGraphloop.log
awk '{print $4}' $file  > counts/$name
numLinesFile=`cat $file | wc -l`
numLinesCounts=`cat counts/$name | wc -l`
echo There are $numLinesFile lines in the original file >> bedGraphloop.log
echo ...and $numLinesCounts lines in the counts file >> bedGraphloop.log
echo ...and $numLinesWindowFile in the genomic window file >> bedGraphloop.log
echo >> bedGraphloop.log
done

for file in counts/*
do
echo $file >> peakCountsFile.txt
done

#######################

### Other genomic Features
G4=$NSNB/2013/03_04-March_April/peakFilesAndIGVtraces_hg19/G4s/genomeWide_GquadruplexMotifs_ABCC_NBDdb_hg19.bed
Genes=$NSNB/2013/07_July/peakFiles/genomeWindows/genes/refSeqGenes.bed 
## These are uniq genes
ORC=$NSNB/2013/03_04-March_April/peakFilesAndIGVtraces_hg19/OtherPeakSets/hg19/DellinoORC1/Dellino_GSM922790_ChIPseq_Orc1_GradientHela_enrichRegions_liftoverHg18ToHg19.bed
CpGislands=$NSNB/2013/03_04-March_April/peakFilesAndIGVtraces_hg19/CpGislands/hg19/CpGIslands_hg19.bed
AladjemMCF7=$NSNB/2013/03_04-March_April/peakFilesAndIGVtraces_hg19/OtherPeakSets/hg19/Aladjem/Aladjem_hg19_originalData_MCF7_Peaks_FDR15_col1-4only.bed
AladjemK562=$NSNB/2013/03_04-March_April/peakFilesAndIGVtraces_hg19/OtherPeakSets/hg19/Aladjem/Aladjem_hg19_originalData_K562_Peaks_FDR15_col1-4only.bed
MesnerBTseq=$NSNB/2013/03_04-March_April/peakFilesAndIGVtraces_hg19/OtherPeakSets/LarryMesnerBTseq_preRelease/liftOver/GM_3reps_RD_bubbles_latest_liftOverHg18ToHg19.bed

echo copying
## move feature files to features subdir
for file in $G4 $Genes $ORC $CpGislands $AladjemMCF7 $AladjemK562 $MesnerBTseq; do
cp $file features/
done

echo bdg for features
## bedgraphs of feature files
for file in features/*.bed
do
name=`stringEdit $file .bed .bedGraph | stringEdit - features/`
echo $name
echo Working with $file >> bedGraphloop.log
echo Naming output $name >> bedGraphloop.log
coverageBed -counts -a $file -b $genomeWindows | sortBed -i - > features/bedGraph/$name
numLinesFile=`cat $file | wc -l`
numLinesCounts=`cat features/bedGraph/$name | wc -l`
echo There are $numLinesFile lines in the original file >> bedGraphloop.log
echo ...and $numLinesCounts lines in the counts file >> bedGraphloop.log
echo ...and $numLinesWindowFile in the genomic window file >> bedGraphloop.log
echo >> bedGraphloop.log
done

echo renaming feature.bdgs
## rename feature files
mv features/bedGraph/genomeWide_GquadruplexMotifs_ABCC_NBDdb_hg19.bedGraph features/bedGraph/G4.bedGraph
########## genes are named fine as is
mv features/bedGraph/Dellino_GSM922790_ChIPseq_Orc1_GradientHela_enrichRegions_liftoverHg18ToHg19.bedGraph features/bedGraph/ORC.bedGraph
mv features/bedGraph/CpGIslands_hg19.bedGraph features/bedGraph/CpGIslands.bedGraph
mv features/bedGraph/Aladjem_hg19_originalData_MCF7_Peaks_FDR15_col1-4only.bedGraph features/bedGraph/Aladjem_MCF7.bedGraph
mv features/bedGraph/Aladjem_hg19_originalData_K562_Peaks_FDR15_col1-4only.bedGraph features/bedGraph/Aladjem_K562.bedGraph
mv features/bedGraph/GM_3reps_RD_bubbles_latest_liftOverHg18ToHg19.bedGraph features/bedGraph/MesnerBTseq.bedGraph


echo counts for features
### counts
for file in features/bedGraph/*.bedGraph
do
name=`stringEdit $file features/bedGraph/ | stringEdit - .bedGraph .txt`
echo $name
echo Working with $file >> bedGraphloop.log
echo Naming it $name >> bedGraphloop.log
awk '{print $4}' $file  > features/counts/$name
numLinesFile=`cat $file | wc -l`
numLinesCounts=`cat features/counts/$name | wc -l`
echo There are $numLinesFile lines in the original file >> bedGraphloop.log
echo ...and $numLinesCounts lines in the counts file >> bedGraphloop.log
echo ...and $numLinesWindowFile in the genomic window file >> bedGraphloop.log
echo >> bedGraphloop.log
done

echo count file for features
for file in features/counts/* 
do 
echo $file >> featureCountsFile.txt
done

