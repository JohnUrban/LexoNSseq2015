#!/usr/bin/env python
from collections import defaultdict
import sys

def baseCount(DNAstring):
    """Takes in DNAstring, returns dictionary of base:count"""
    baseDict = defaultdict(int)
    DNAstring = DNAstring.upper()
    for base in DNAstring:
        baseDict[base] = baseDict[base] + 1
    return baseDict
    
        
def extractReadFromSamLine(line):
    """Take in and evaluates line from a SAM file.
    Assumes line is not a header line.
    Returns sequence of read from column 10"""
    return line.split()[9]

def makeGCdistDict(samFileConnection, readLen=50, numN=0):
    """Opens up SAM file for reading and datafile for writing
    Returns  dictionary of how many reads had GC count of x
    GC counts can be divided by readlength
    Read counts can be divided by the sum of counted reads
        (counted reads are reads satisfying minAllowedATGCcount rule established by numN)"""
    sam = samFileConnection
    minAllowedATGCcount = readLen-numN
    GCdict = defaultdict(int)
    for i in range(0, readLen+1): ## ensure stable output (i.e. for any read pile, it returns count for all possible GCcounts even if GCcount was 0)
        GCdict[i] = 0
    for line in sam:
        read = extractReadFromSamLine(line)
        if read != None:
            baseDict = baseCount(read)
            if baseDict['A'] + baseDict['C'] + baseDict['G'] + baseDict['T'] >= minAllowedATGCcount:
                GCcount = baseDict['G'] + baseDict['C']
                GCdict[GCcount] += 1
    sam.close()
    return GCdict


def writeOutGCDictHist(GCdict):
    """Takes in the GCdict and writes it out as file with 2 columns: numGCs and numReadsWithThisNumGCs"""
    for key in sorted(GCdict):
        sys.stdout.write(str(key) + "\t" + str(GCdict[key]) + "\n")

defaultmsg='''
Usage: GCcontentInSamReads.py inputSAMFile readLength numN

Try 'GCcontentInSamReads.py -h' for more information.
'''
helpmsg='''
Author: John M. Urban (Brown University)

Usage: GCcontentInSamReads.py inputSAMFile readLength numN


Parameters:
    inputSAMFile - file of READS in SAM format; assumes no header lines.
    readLength - the length of reads in file, assumes all same length (e.g. 50 bp)
    numN - number of non(ACGT) letters (e.g. N, R, Y) allowed to occur in read -- any more and the read is ignored.
        It is recommended to make numN = 0 (i.e. ignore reads with any non-ATGC letters, which is usually a very small proportion in Illumina data)

Quick Start Example: 
    $ GCcontentInSamReads.py file.sam 50 0

Examples in different scenarios:
    If starting with a SAM file that does NOT have header lines, just do:
        $ GCcontentInSamReads.py file.sam readlen num
        
    If starting with a SAM file that DOES have header lines, just do:
        $ grep -v ^@ file.sam | GCcontentInSamReads.py - readlen numN
        
    If starting with a BAM file, just do:
        $ samtools view file.bam | GCcontentInSamReads.py - readlen numN

Output has 2 columns: 
    1. G+C count from 0 to readLength
    2. numReads with given G+C count
 Multiply each number in the first column by 100, then divide by readlength to get percent GC.
 Multiply each number in the second column by 100, then divide by the sum of second column to get the percent of total counts for each

Try at commandline:
    $ readlen=50
    $ numN=0
    $ samtools view example-GC.bam | ./GCcontentInSamReads.py - $readlen $numN > numGCvsNumRead.txt
    $ total=`awk '{SUM+=$2} END {print SUM}' numGCvsNumRead.txt`
    $ awk -v readlen=$readLen -v total=$total 'OFS=\"\\t\" {print 100*$1/readlen, 100*$2/total}' numGCvsNumRead.txt > pctGCvsPctCountedReads.txt

Visualize in R:
    - Either resulting file above (numGCvsNumRead.txt or pctGCvsPctCountedReads.txt)
        can be brought into R for visualization
    - To start, in R try (while in directory with given file):
        > data <- read.table("numGCvsNumRead.txt")
        > plot(data$V1, data$V2)
        '''
helpsignals = ["-h", "--help", "-help"]

##If run from the linux commandline then execute the following
#### See 6.1.1 here for more info on how this works: http://docs.python.org/2/tutorial/modules.html

if __name__ == "__main__":

    if len(sys.argv) == 1:
        print defaultmsg

    elif sys.argv[1] in helpsignals:
        print helpmsg
        quit()

    else:
        ##Parameters
        samFile = sys.argv[1]
        if samFile == "-" or samFile == "stdin":
            samFileConnection = sys.stdin
        else:
            samFileConnection = open(samFile, 'r')
        readLen = int(sys.argv[2])
        numN = int(sys.argv[3])
        
        ##execute
        GCdict = makeGCdistDict(samFileConnection, readLen=readLen, numN=numN)
        writeOutGCDictHist(GCdict)


    quit()
