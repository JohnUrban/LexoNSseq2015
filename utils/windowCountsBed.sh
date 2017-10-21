#!/bin/bash

## John Urban, July 6, 2017; v1

## for BEDtools >=  v2.26.0
function help {
    echo "Usage: scriptname -w:s:g:e:f:o:c:bvh
        -w with argument = windowsize (e.g. 100000)
        -s with argument = stepsize (e.g. 100000)
        -g with argument = .genome file for BEDtools (e.g. hg19.genome)
        -e with argument = excluded regions BED file (e.g. excludedRegions.bed) often containing gaps and/or entire chromomes to exclude (chr 0 len)
        -f with argument = all files to work on. if more than one file, put all inside quotes. if doing wild card such as *, also put inside quotes.
        -o with argument = outdirectory (will be made if non-existent); default is pwd
        -c with argument = /path/to/dir that has cor-analysis.R (e.g. LexoNSseq2015/utils); this results in also running correlation analysis (only makes sense with multiple files)
        -b save intermediate bedGraphs
        -v verbose
        -h help - returns this message; also returns this when no arguments given
"
}


if [ $# -eq 0 ]; then help; exit; fi

##Defaults:
BDG=false
VERBOSE=false
HELP=false
OUTDIR="."
COR=false

while getopts "w:s:g:e:f:o:c:bvh" arg; do
    case $arg in
        w) WINDOW=$OPTARG;;
        s) STEP=$OPTARG;;
        g) G=$OPTARG;;
        e) EXCLUDE=$OPTARG;;
        f) FILES=$OPTARG;;
        o) OUTDIR=$OPTARG;;
        b) BDG=true;;
        c) CORDIR=$OPTARG; COR=true;;
        v) VERBOSE=true;;
        h) HELP=true;;
        *) help;;
    esac
done

## HELP
if ${HELP}; then help; exit; fi

## OUTDIR?
if [ ! -z ${OUTDIR} ]; then if [ ! -d ${OUTDIR} ]; then mkdir ${OUTDIR}; fi; fi

## ENTER VARIABLES HERE
#FILES="G4sfromquadParser_hg19.bed"
#EXCLUDE=excludedRegions.bed
#G=hg19.genome
#WINDOW=100000
#STEP=100000

## Could do in one command as:
function single_nobdg {
  bedtools makewindows -g <(sort -k1,1 $G) -w 100000 -s 100000 | intersectBed -sorted -v -a - -b <(sortBed -i $EXCLUDE ) | awk -v "KEEP=$WINDOW" '$3-$2==KEEP' | coverageBed -counts -sorted -a - -b <( sortBed -i $BEDFILE ) | sortBed -i - | awk '{print $4}' > $COUNTS.txt
}

## Could also store intermediate bedGraph with tee
function single_withbdg {
  bedtools makewindows -g <(sort -k1,1 $G) -w 100000 -s 100000 | intersectBed -sorted -v -a - -b <(sortBed -i $EXCLUDE ) | awk -v "KEEP=$WINDOW" '$3-$2==KEEP' | coverageBed -counts -sorted -a - -b <( sortBed -i $BEDFILE ) | sortBed -i - | tee ${COUNTS}.bedGraph | awk '{print $4}' > ${COUNTS}.txt
}

## if have windows bed file used on coverageBed and a bunch of different features counted in those windows
## then better to make the windows file to avoid redundant computation
## really don't need bedGrapg b/c can just store the counts and paste windows and counts together later if needed: paste windows.bed counts.txt
function makewindows {
  bedtools makewindows -g <(sort -k1,1 $G) -w 100000 -s 100000 | intersectBed -sorted -v -a - -b <(sortBed -i $EXCLUDE ) | awk -v "KEEP=$WINDOW" '$3-$2==KEEP' | sortBed -i - 
}

function get_counts_nobdg {
 WINDOWS=${1}
 coverageBed -counts -sorted -a ${WINDOWS} -b <( sortBed -i $BEDFILE ) | sortBed -i - | awk '{print $4}' > ${COUNTS}.txt
}

function get_counts_withbdg {
 WINDOWS=${1}
 coverageBed -counts -sorted -a ${WINDOWS} -b <( sortBed -i $BEDFILE ) | sortBed -i - | tee ${COUNTS}.bedGraph | awk '{print $4}' > ${COUNTS}.txt
}

function verbose {
    msg=${1}
    if ${VERBOSE}; then echo ${msg}; fi
}

## testing
##for file in $FILES; do echo $file; done

## HOW MANY BED FILES ANALYZING?
NUMFILES=`echo $FILES | awk '{print NF}'`

if [ $NUMFILES -eq 1 ]; then
  BEDFILE=$FILES
  BASE=`basename ${BEDFILE}`
  COUNTS=${OUTDIR}/${BASE}.counts
  verbose "single, BDG=${BDG}, ${BASE}"
  if $BDG; then
    single_withbdg
  else
    single_nobdg
  fi
else
  WINDOWSFILE=windows.w${WINDOW}.s${STEP}.bed
  makewindows > ${WINDOWSFILE}
  for BEDFILE in $FILES; do
    BASE=`basename ${BEDFILE}`
    COUNTS=${OUTDIR}/${BASE}.counts
    verbose "multi, BDG=${BDG}, ${BASE}"
    if $BDG; then
      get_counts_withbdg ${WINDOWSFILE}
    else
      get_counts_nobdg ${WINDOWSFILE}
    fi
  done
fi

if $COR; then
  verbose "Running correlation analyses..."
  paste ${OUTDIR}/*counts.txt | cat <(ls ${OUTDIR}/*counts.txt | paste -sd"\t" - ) - > ${OUTDIR}/counts.dataframe.txt
  for method in pearson spearman; do
    Rscript ${CORDIR}/cor-analysis.R ${OUTDIR}/counts.dataframe.txt ${method} ${OUTDIR}/cor.${method}
  done
fi

