#!/bin/bash
#### Made 9/3/13
#### This was made to somewhat replace "allByAllIntersectBed"
#### This only looks at 100%-vs-100% overlaps
#### Another function will look at ordered sets (e.g. top40-vs-all)
#### This assumes BED files of any denomination (only requires first 3 cols to be BED style)
###### i.e. 0-based start


### Define functions
function help {
if [ $# -eq 0 ] || [ $1 == "--help" ]; then
  echo
  echo "Usage:  scriptname  pathFile   counts|percents"
  echo "pathFile is file that has all paths to BED files being analyzed"
  echo
  exit
fi
}

function overlapCounts {
while read fileA; do
 ROW=""
 while read fileB; do
  COUNT=`intersectBed -u -a $fileA -b $fileB | wc -l`
  ROW=`echo -e $ROW "\t" $COUNT`
 done < $1
 echo $ROW
done < $1
}

function overlapPercents {
while read fileA; do
 ROW=""
 numPeaks=`cat $fileA | wc -l`
 while read fileB; do
  COUNT=`intersectBed -u -a $fileA -b $fileB | wc -l`
  PERCENT=`echo 100*${COUNT}/${numPeaks} | bc -l`
  ROW=`echo -e $ROW "\t" $PERCENT`
 done < $1
 echo $ROW
done < $1
}

### Execute commandline 
help $@
if [ "$2" == "counts" ]; then overlapCounts $1; exit; fi
if [ "$2" == "percents" ]; then overlapPercents $1; exit; fi
