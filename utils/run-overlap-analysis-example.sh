#!/bin/bash

## Writing this script 07/06/2017 to help Miiko

## abspath.py part of sciaratools
## intersectBedLoop.sh now part of LexoNSseq2015/utils
## 

## CHANGE ALL VARIABLES AS NEEDED
WD=analysis
FILES="lg0peaks/*"
NARROWPEAK=/Users/johnurban/searchPaths/github/LexoNSseq2015/utils/narrowPeakAnalyzer.R
OVLPSTATS=/Users/johnurban/searchPaths/github/LexoNSseq2015/utils/ovl-analysis.R
GENOME=hg19.genome
GAPS=gaps_hg19_sorted.genome


## IF DIR ALREADY EXISTS, DO NOT OVERWRITE -- EXIT
if [ -d $WD ]; then echo "Analysis directory already exists... Exiting..." ; exit; fi

## ELSE CARRY ON
mkdir ${WD}
abspath.py ${FILES} > ${WD}/paths.txt
for f in ${FILES}; do basename $f .narrowPeak; done > ${WD}/names.txt

## ENSURE ChrY and chrM are taken out of this analysis
grep -v -E 'chrY|chrM' ${GENOME} | sort -k1,1 > ${WD}/hg19.genome
GENOME=`abspath.py ${WD}/hg19.genome`
grep -v -E 'chrY|chrM|_' ${GAPS} | sortBed -i - > ${WD}/hg19.gaps
GAPS=`abspath.py ${WD}/hg19.gaps`

## get approximate number of connected components by doing complementBed on the gapsFile and counting the number of regions that are reported (which are regions broken up by gaps)
NCONCOMP=`complementBed -i $GAPS -g $GENOME | wc -l`

cd ${WD}
intersectBedLoop.sh paths.txt counts > counts.txt

OVLPTABLE=`abspath.py counts.txt`
NAMES=`abspath.py names.txt`
PATHS=`abspath.py paths.txt`
MINOVLP=1

Rscript $OVLPSTATS $NARROWPEAK $GENOME $GAPS $NCONCOMP $OVLPTABLE $NAMES $PATHS $MINOVLP

