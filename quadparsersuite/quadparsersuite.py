#!/usr/bin/env python
## this has been modified from original (on June 9, 2014)
## Original script Copy/pasted from:
## http://bioinformatics-misc.googlecode.com/svn-history/r16/trunk/quadparser.py
## on 4/9/13 (11:42PM)
## Also see: http://code.google.com/p/bioinformatics-misc/source/browse/trunk/quadparser.py?spec=svn16&r=16
## Script by dario.beraldi
## Last updated Dec 25, 2011 by DB.

import re
import sys
import string
import argparse

parser = argparse.ArgumentParser(description="""

DESCRIPTION
    
    This script descends from Dario Beraldi's quadparser.py located here:
        http://bioinformatics-misc.googlecode.com/svn-history/r16/trunk/quadparser.py
    Major functional differences are that quadparsersuite.py:
    - is able to handle FASTQ files (in addition to FASTA files)
    - was designed with giant FASTQ files in mind
    - allows more control over output
    - defaults to BED output (unless --counts is used)
    - allows more control over the G4-motifs
        - minimum number of Gs per G-tract (--minG, default 3)
        - maximum number of nucleotides in between G-tracts (--maxN, default, 7)

    Notes about things under the hood and output:
        - sequences are read in memory one at a time. Thus, if the human genome is given, then 1
            chromosome at a time. This works fine on a regular Mac Book Pro with 4-8 GB RAM,
            though may have issues if there is not a lot of free memory. Species with larger
            chromosomes may require more resources. Alternatively, modifications to the script
            can allow it to parse through smaller chunks of each sequence at a time.
        - A difference in this script from the original quadparser.py is it writes out entries
            as it encounters them rather than sorting them first. Since it processes the forward regex
            followed by the complement regex (when applicable) it will write out all positive strand entries
            of a given sequence in order followed by all negative strand entries in order. Thus, positive
            and negative strand motifs are not sorted with respect to each other.
            - One way to sort is if name,start,end are used as first 3 columns, sortBed from BEDtools can sort the
                BED file if needed. For example:
                $ quadparsersuite.py -i hg19.fa -o name,start,end,strand | sortBed -i - > hg19G4.bed
            - The unix command "sort" is also fine if not BED-like. For example, to sort on name folowed by start positions:
                $ quadparsersuite.py -i hg19.fa -o name,start,strand | sort -k,k -k2,2n > hg19G4starts.txt

    Modifications to quadparser.py to derive quadparsersuite.py including this description and
    command-line protocols relevant to the 2015 paper "Characterizing and controlling intrinsic
    biases of Lambda exonuclease in nascent strand sequencing reveals phasing between nucleosomes
    and G-quadruplex motifs around a subset of human replication origins", were made by John M. Urban
    (Brown University).

    Search for matches to any regular expression (regex) in a FASTA or FASTQ file and return a
    BED/BED-like file with the coordinates of the match and the matched sequence itself
    (optional --reportseq). If the regex includes letters other than ATGCU, it is advised to use "--noreverse"
    as complement rules no longer apply. Note that this can take a FASTA file with any sequence line
    width rule (not just 2 line entry format of >seqname,seq pairs). However, it will only properly
    analyze FASTQ files with the standard 4 line entry format (seqname,seq,+seqname,qscores).
    
    With defaults, quadparsersuite.py searches for putative intra-strand quadruplexes on forward
    and reverse strand using the quadruplex rule described at
    http://en.wikipedia.org/wiki/G-quadruplex#Quadruplex_prediction_techniques
    and in:
        Huppert JL, Balasubramanian S. 2005.
        Prevalence of quadruplexes in the human genome.
        Nucleic Acids Res. 33: 2908-2916.
    Search only the given sequence ("forward strand") by using --noreverse.
    
    The defualt regex is '([gG]{3,}\w{1,7}){3,}[gG]{3,}' and its complement '([cC]{3,}\w{1,7}){3,}[cC]{3,}'
    produce the same output as in http://www.quadruplex.org/?view=quadbaseDownload
    and/or https://nonb-abcc.ncifcrf.gov/apps/site/default
    
    Output BED file has default columns:
    1. Name of fasta sequence (e.g. chromosome)
    2. Start of the match
    3. End of the match
    4. Strand 
    5. Optional Matched sequence (--reportseq/-s)

    These can be changed with --outformat/-o which allows you to report any subset of
        "name,start,end,strand,seq" in any order.
    For example "-o start" will result in only the start of each match being reported.

    If --counts is used, default columns are:
    1. name
    2. pos strand count
    3. neg strand count
    4. total count
        

    G4 Start Site Counts in Illumina reads:
    One way to obtain this at the command-line using quadparsersuite.py is:
        $ quadparsersuite.py -i $reads_file.fq --noreverse -o start | sort -n | uniq -c > g4StartSiteCounts.txt


    This returns counts, which can then be normalized by the number of reads searched and
    multiplied by 1 million to get the G4 Start Site Counts per Million Reads (CPMR).
    Note that "--noreverse" means to only count G4s on the strand given (i.e. the Illumina read),
    which corresponds to 5' ends of fragments.


    Create a BED file of G4s on both strands for your favorite genome - view in Browser,
    use in BEDtools, etc:
        $ quadparsersuite.py -i $reads_file.fq > g4s.bed


    Analyses on g4s.bed (created above) with BEDtools
    Intersect counts with peak sets
        How many peaks overlap at least one G4?
            $ intersectBed -u -a peaks.bed -b g4s.bed | wc -l
            
        How many G4s overlap at least one peak?
            $ intersectBed -u -a g4s.bed -b peaks.bed | wc -l

            
    What are G4 counts in 100 kb bins in the genome:
        Make bin file:
            $ bedtools makewindows -g target.genome -w 100000 > 100kbWindows.bed
            
        Or make bin file retaining only 100 kb bins that do not overlap gaps and are not on chrY or chrM.
            $ bedtools makewindows -g target.genome -w 100000 | intersectBed -v -a - -b gaps.genome.bed | awk '{if ($3-$2 == 100000) print}' | grep -v ^chrY | grep -v ^chrM > 100kbBinsFiltered.bed

        Get G4 counts in each bin:
            $ coverageBed -counts -a g4s.bed -b 100kbBins.bed | sortBed -i - > g4Counts.bedGraph
            
        Get peak counts in same bins:
            $ coverageBed -counts -a peaks.bed -b 100kbBins.bed |sortBed -i - > peakCounts.bedGraph
            
        Test correlation between peak counts and G4 counts in the same bins by bringing the
        bedGraph or just the counts in column 4 into R and doing:
            $ cor(g4Counts, peakCounts, method="pearson")
            $ cor(g4Counts, peakCounts, method="spearman")
            

    G4 CPMR:
    This metric offers a single number instead of a number of each possible start site position
    in a read. One can obtain it simply by summing the counts over each start site position.
    For example:
        $ quadparsersuite.py -i $reads_file.fq --noreverse -o start | wc -l

    Another way to get both metrics would be:
        $ quadparsersuite.py -i $reads_file.fq --noreverse > g4s.bed
        $ wc -l g4s.bed     #G4 CPMR
        $ cut -f 2 g4s.bed | sort -n | uniq -c   #G4 start site CPMR

    
    Supp notes:
    hg19 gaps can be found at UCSC Table Browser:
    http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=411961815_wAeC9o7jDZHSWZKpzI8tIrVsQMS4&clade=mammal&org=Human&db=hg19&hgta_group=map&hgta_track=gap&hgta_table=0&hgta_regionType=genome&position=chr21%3A33031597-33041570&hgta_outputType=primaryTable&hgta_outFileName=
    

    NOTE Nov 2017:
    New i-motif paper (featuring Julian Hupper) gives a regex for i-motifs:
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5605235/
    motif = C5(N1-19C5)3
    written as regex = [Cc]{5,}(\w{1,19}[Cc]{5,}){3,}

    """, formatter_class= argparse.RawTextHelpFormatter)


parser.add_argument('--regex', '-r',
                   type= str,
                   help='''Regex to be searched in the fasta input.
Matches to this regex will have + strand. This string passed to python
re.compile(). The default regex is '([gG]{3,}\w{1,7}){3,}[gG]{3,}' which searches
for G-quadruplexes.
                                   
                   ''',
                   default=False)

parser.add_argument('--regexrev', '-R',
                   type= str,
                   help='''The second regex to be searched in fasta input.
Matches to this regex will have - strand.
By default (None), --regexrev will be --regex complemented by replacing
'actguACTGU' with 'tgacaTGACA'.  
                                   
                   ''',
                   default= None)

parser.add_argument('--minG', '-g',
                   type= int,
                   help='''minG is the minimum number of Gs in a G tract.
A G4 is typically defined as: ([gG]{3,}\w{1,7}){3,}[gG]{3,}
As such, the default minG value is 3.
                   ''',
                   default=3)

parser.add_argument('--maxN', '-n',
                   type= int,
                   help='''maxN is the maximum number of number of Ns in loops between G tracts.
A G4 is typically defined as: ([gG]{3,}\w{1,7}){3,}[gG]{3,}
As such, the default maxN value is 7.
Recently people have also often used maxN=15 -- i.e. ([gG]{3,}\w{1,15}){3,}[gG]{3,}
Note, though, as G4 motifs get longer, they become less likely to actually form G4s,
especially in the context of dsDNA.
                   ''',
                   default=7)

parser.add_argument('--inputFile', '-i',
                   type= str,
                   help='''Input file in fasta format containing one or more
sequences. Use '-' to get the name of the file from stdin
                                   
                   ''',
                    required=True)
fileType = parser.add_mutually_exclusive_group()

fileType.add_argument('--fasta', '-a', default=False,
                   action= 'store_true',
                   help='''Input file in fasta (fa) format containing one or more sequences.''')

fileType.add_argument('--fastq', '-q',
                   action= 'store_true', default=False,
                   help='''Input file in fastq (fq) format containing one or more sequences.''')

parser.add_argument('--noreverse',
                   action= 'store_true',
                   help='''Do not search the reverse (-) strand. I.e. do not use
the complemented regex (or --regexrev/-R). Use this flag to search protein
sequences.
                                  
                   ''')
parser.add_argument('--reportseq', '-s',
                    action= 'store_true', default=False,
                    help='''Report sequence of reg exp in output.                                   
                   ''')
parser.add_argument('--outformat', '-o',
                    type=str, default='name,start,end,strand',
                    help='''Provide comma-separated list of desired output infomation.
Options are name (sequence name), start (start of match), end (end of match),
strand (strand of match +/-), seq (sequence of match).
Default = 'name,start,end,strand'. --reportSeq/-s option changes default to: 'name,start,end,strand,seq'
Any other combination can be provided.
When using --counts, defaults to name,pos,neg
                   ''')

parser.add_argument('--counts', '-c',
                    action= 'store_true', default=False,
                    help='''Report count for number of matches in each sequence in file. 
                   ''')

args = parser.parse_args()

" --------------------------[ Check and pare arguments ]---------------------- "


""" Identify filetype as fasta or fastq """
if not args.fasta and not args.fastq:
    if args.inputFile == "-":
        sys.exit('\nWhen using stdin, need to specify --fasta (-a) or --fastq (-q)\n')
    else:
        test = open(args.inputFile, 'r')
        filetypetest = test.next()[0]
        if filetypetest == '>':
            args.fasta = True
        elif filetypetest == '@':
            args.fastq = True
        else:
            test.close()
            sys.exit('\nCould not figure out file type. Try specifying --fasta (-a) or --fastq (-q)\n')
        test.close()
## create regex
if not args.regex:
    ##default='([gG]{3,}\w{1,7}){3,}[gG]{3,}'
    args.regex = '([gG]{' + str(args.minG) + ',}\w{1,' + str(args.maxN) + '}){3,}[gG]{' + str(args.minG) + ',}'


if args.reportseq:
    args.outformat = 'name,start,end,strand,seq'
if args.counts:
    args.outformat = 'name,pos,neg'

""" Reverse forward match """
intab=  'actguACTGU'
outtab= 'tgacaTGACA'
if args.regexrev is None:
    transtab = string.maketrans(intab, outtab)
    regexrev = args.regex.translate(transtab)
else:
    regexrev= args.regex

""" Opening file connection or dealing with stdin """
if args.inputFile == '-':
    ref_seq_fh = sys.stdin
    if len(args.inputFile) > 1:
        sys.exit('\nquadpareser.py: Only one input file at a time can be processed:\n--fasta/-f: %s\n' %(args.input))
else:
    ref_seq_fh = open(args.inputFile)

""" give values to variables """
psq_re_f= re.compile(args.regex)
psq_re_r= re.compile(regexrev)

    
" ------------------------------[  Functions ]--------------------------------- "

def parseSequence(ref_seq, ref_seq_name, fwd_re, rev_re, noreverse=False):
    for m in re.finditer(fwd_re, ref_seq):
        formatDict = {'name':ref_seq_name, 'start':m.start(), 'end':m.end(), 'strand':'+', 'seq':m.group(0)}
        print '\t'.join(str(x) for x in [formatDict[e] for e in args.outformat.split(',')])
    if noreverse is False:
        for m in re.finditer(rev_re, ref_seq):
            formatDict = {'name':ref_seq_name, 'start':m.start(), 'end':m.end(), 'strand':'-', 'seq':m.group(0)}
            print '\t'.join(str(x) for x in [formatDict[e] for e in args.outformat.split(',')])

def countSequence(ref_seq, ref_seq_name, fwd_re, rev_re, noreverse=False):
    poscount = 0
    negcount = 0
    for m in re.finditer(fwd_re, ref_seq):
        poscount += 1
    if noreverse is False:
        for m in re.finditer(rev_re, ref_seq):
            negcount += 1
    formatDict = {'name':ref_seq_name, 'pos':poscount, 'neg':negcount}
    print '\t'.join(str(x) for x in [formatDict[e] for e in args.outformat.split(',')])
    
def parseFasta(ref_seq_fh, psq_re_f, psq_re_r, function, noreverse=False):
    ref_seq = []
    line = (ref_seq_fh.readline()).strip()
    ref_seq_name = re.sub('^>', '', line)
    line = (ref_seq_fh.readline()).strip()
    while True:
        while line.startswith('>') is False:
            ref_seq.append(line)
            line= (ref_seq_fh.readline()).strip()
            if line == '':
                break
        ref_seq= ''.join(ref_seq)

        function(ref_seq, ref_seq_name, fwd_re=psq_re_f, rev_re=psq_re_r, noreverse=args.noreverse)
        
        ref_seq_name = re.sub('^>', '', line)
        ref_seq= []
        line= (ref_seq_fh.readline()).strip()
        if line == '':
            break


def parseFastq(ref_seq_fh, psq_re_f, psq_re_r, function, noreverse=False):
    line = (ref_seq_fh.readline()).strip()
    ref_seq_name = re.sub('^@', '', line)
    ref_seq = (ref_seq_fh.readline()).strip()
    while True:
        function(ref_seq, ref_seq_name, fwd_re=psq_re_f, rev_re=psq_re_r, noreverse=args.noreverse)
        ref_seq_fh.readline() ## skip + line
        ref_seq_fh.readline() ## skip base call quality line
        line = (ref_seq_fh.readline()).strip()   
        ref_seq_name = re.sub('^@', '', line)
        ref_seq = (ref_seq_fh.readline()).strip()
        if line == '':
            break

'''----------------------[ Execute ]-------------------------------'''

if args.counts:
    function = countSequence
else:
    function = parseSequence
    
if args.fasta:
    parseFasta(ref_seq_fh, psq_re_f, psq_re_r, function, noreverse=args.noreverse)
elif args.fastq:
    parseFastq(ref_seq_fh, psq_re_f, psq_re_r, function, noreverse=args.noreverse)

sys.exit()
