#!/usr/local/bin/python
import sys
from collections import defaultdict
import random

##if len(sys.argv) == 2:
##    method = sys.argv[1]
##else:
##    method = "sum"
if len(sys.argv) == 1:
    print "Method required as arg1: sum, tv, both, tvToSum"
method = sys.argv[1]

if method == "sum":
    sumSignal = defaultdict(int)
    numScoresForPosition = defaultdict(int)
    for line in sys.stdin:
        line = line[:-1].split() #[relPos, score]
        pos = int(line[0])
        score = float(line[1])
        sumSignal[pos] += score
        numScoresForPosition[pos] += 1
    for pos in sorted(sumSignal.keys()):
        entry = str(pos) + "\t" + str(sumSignal[pos]) + "\t" + str(sumSignal[pos]/numScoresForPosition[pos]) + "\t" + str(numScoresForPosition[pos]) + "\n" 
        sys.stdout.write(entry)
    ## Note: some scenarios may require pos to be a float -- right now it is an int
    
elif method == "tv":
    pass
    ## this will allow you to deal with median, min, max, sd, mean, etc in R
    ## it also lets you deal intelligently with missing data (not all missing data should be considered 0)
    targets = defaultdict(list)
##    numScoresForPosition = defaultdict(int)
    for line in sys.stdin:
        line = line[:-1].split() ## [peakNum, relPos, score]
        entryNumber = int(line[0])
        score = line[2] ##score kept as string purposely
        targets[entryNumber].append(score)
##        numScoresForPosition[pos] += 1
    for target in sorted(targets.keys()):
        entry = ','.join(targets[target]) + "\n"
        sys.stdout.write(entry)
elif method == "both":
    targets = defaultdict(list)
    sumSignal = defaultdict(int)
    numScoresForPosition = defaultdict(int)
    for line in sys.stdin:
        line = line[:-1].split() ## [peakNum, relPos, score]
        entryNumber = int(line[0])
        score = line[2] ##score kept as string purposely
        targets[entryNumber].append(score)
        pos = int(line[1])
        if score != ".":
            score = float(score) ##now score is a float
            sumSignal[pos] += score
            numScoresForPosition[pos] += 1       
    for pos in sorted(sumSignal.keys()):
        entry = str(pos) + "\t" + str(sumSignal[pos]) + "\t" + str(sumSignal[pos]/numScoresForPosition[pos]) + "\t" + str(numScoresForPosition[pos]) + "\n" 
        sys.stdout.write(entry)
    for target in sorted(targets.keys()):
        entry = ','.join(targets[target]) + "\n"
        sys.stdout.write(entry)
elif method == "tvToSum":
    ## assumes input lines are csv tv input
    ## assumes all same length
    ## sys.argv[2] should be window size to each side of point (summit)
    ## e.g. sys.argv[3] = 1000 -- this is subtracted from each index in list
    window = int(sys.argv[2])
    sumSignal = defaultdict(int)
    numScoresForPosition = defaultdict(int)
    for line in sys.stdin:
        line = line.rstrip().split(",")
        for i in range(len(line)):
            score = line[i]
            if score != ".":
                relPos = i-window
                score = float(score)
                sumSignal[relPos] += score
                numScoresForPosition[relPos] += 1
    for pos in sorted(sumSignal.keys()):
        entry = str(pos) + "\t" + str(sumSignal[pos]) + "\t" + str(sumSignal[pos]/numScoresForPosition[pos]) + "\t" + str(numScoresForPosition[pos]) + "\n" 
        sys.stdout.write(entry)

    ## tvToShuffSum will permute each tv screen set of numbers before moving on
    ## this keeps the nsome signal distribution the same but randomly shuffled
elif method == "scrambledTvToSum":
    ## virtually the same as tvToSum except it shuffles scores
    ## tvToShuffSum will permute each tv screen set of numbers before moving on
    ## this keeps the nsome signal distribution the same but randomly shuffled
    window = int(sys.argv[2])
    sumSignal = defaultdict(int)
    numScoresForPosition = defaultdict(int)
    for line in sys.stdin:
        line = line.rstrip().split(",")
        ##SHUFFLE
        random.shuffle(line)
        for i in range(len(line)):
            score = line[i]
            if score != ".":
                relPos = i-window
                score = float(score)
                sumSignal[relPos] += score
                numScoresForPosition[relPos] += 1
    for pos in sorted(sumSignal.keys()):
        entry = str(pos) + "\t" + str(sumSignal[pos]) + "\t" + str(sumSignal[pos]/numScoresForPosition[pos]) + "\t" + str(numScoresForPosition[pos]) + "\n" 
        sys.stdout.write(entry)
        
