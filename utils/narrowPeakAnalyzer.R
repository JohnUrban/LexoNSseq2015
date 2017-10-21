## Created by John Urban
## This suite contains:
##  __FUNCTION___________________________DATE_CREATED___
## - readNarrowPeak                     (early 8/2013)
## - readGenomeFile                     (early 8/2013)
## - readGapFile                        (early 8/2013)
## - readNarrowPeakPlusNucBed           (early 8/2013)
## - readBdg                            (June/2014)
## - writeBed                           (8/26/2013)
## - readBed3                           (early 8/2013)
## - readBed3PlusNucBed                 (early 8/2013)
## - readBed5                           (early 8/2013)
## - readBed5PlusNucBed                 (early 8/2013)
## - readFilePaths                      (9/3/13)
## - writeFilePaths                     (9/5/13)
## - histFE                             (early 8/2013)
## - histpVal                           (early 8/2013)
## - histqVal                           (early 8/2013)
## - histPeakLengths                    (early 8/2013)
## - mappableChromosomeSizes            (9/3/13)
## - mappableGenomeSize                 (9/3/13)
## - numPeaksByChrTable                 (early 8/2013, updated 9/3/13)
## - plotNumPeaksByChr                  (early 8/2013)
## - interPeakDistances                 (early 8/2013)
############## Closest feature analysis suite ###################
### No specific functions made for this analysis yet as the closestBed outputs are not necessarily standardized
### Do closestBed with summits, centers, and/or other single bp points (i.e starts)
### One just needs to read closestBed output in with read.table()
### Then to get distances subtract the start of the point (fileA) from the start of the feature (fileB)
###     e.g. ORC centers and NS summits might come in a 6 column table with: ORC[chr start end] NS[chr start end]
###           These columns are just names V1-V6
###           Do: distancesFromORC <- table$V5-table$V2
###           Anything to right (3') of ORC will be a positive distance. 
###           Anything to left (5') of ORC will be a negative distance. ORC is at 0.
############## Signal Around Points Suite #######################
## - signalAroundPoints                 (8/23-26/13) ##Needs work
## - extendPoints                       (8/26/2013)
## - removeRegionsThatFallOffChrStart   (8/26/2013)
## - removeRegionsThatFallOffChrEnd     (8/26/2013)
## - removeRegionsThatOverlapGaps       (8/26/2013)
## - nameBedRegions                     (8/26/2013)
## - getBinsInsideBedRegions            (8/26/2013)
## - bigWigAverageSignalOverBed         (8/26/2013) 
## - makeRowPile                        (8/26/2013) ##Needs work
############## Overlap Analysis Suite #######################
## - readNamesFile                      (9/3/13)
## - makeNamesDataFrame                 (9/3/13)
## - readOverlapTable                   (9/3/13)
## - binomialOverLapPval                (9/3/13)
## - writeAllTables                     (9/4/13)
## - pairwiseBinomialInfo               (9/3/13)
## - writePairwiseBinomialInfo          (9/4/13)
## - plotExpVsObsOverlap                (9/3/13)
## - plotObsExpOverlapRatio             (9/3/13)
## - plotSetSizeRatios                  (9/5/13)
## - plotSetSizeRatiosWithApproxMaxPercentMax (9/5/13)
## - plotApproxMaxWithPercentMax        (9/5/13)
## - plotlog10BinomialpVal              (9/5/13)
## - getAvgOverDiagonalMatrix           (9/3/13)
## - matrixPlot                         (9/3/13)
## - matrixCluster                      (9/3/13)
## - 
## - 
############## Density Correlation Suite #######################
## - 
## - 
## -
## Original version created in early August 2013
## Updated 8/23/13 -- signalAroundPoints started
## Updated 8/26/13 -- signalAroundPoints finished and 'module-ized' -- each step broken into its own function (see above function dates)
## Updated 9/13/13 -- added mappableChromosomeSizes, mappableGenomeSize, readNamesFile, makeNamesDataFrame, readOverlapTable, readFilePaths... Modified numPeaksByChrTable to use mappableChromosomeSizes (rather than the code written out).
##                 -- started binomialOverLapPval

#narrowPeak Format
#1 = chr
#2 = start
#3 = end
#4 = name
#5 = score (-10*log10(qval or pval)) in MACS2
#6 = strand
#7 = signalValue (FE)
#8 = -log10(pval) (-1 if none)
#9 = -log10(qvalue)
#10 = relative summit/point source

##Dependencies
library(lattice)

### Functions
readNarrowPeak <- function(narrowPeakFile){
  #Read Data in
  data <- read.table(narrowPeakFile)
  
  #Ensure class of each column is as it should be
  classes <- c("factor","integer","integer","factor","integer","factor","numeric","numeric","numeric","integer")
  for(i in 1:dim(data)[2]){class(data[,i]) <- classes[i]}  

  #Name columns according to narrowPeak column names
  colNames <- c("chr", "start", "end", "name", "score", "strand", "FE", "pValScore", "qValScore", "relSummit")
  colnames(data) <- colNames
  
  # Output data
  return(data)
}

readGenomeFile <- function(genomeFile){
  data <- read.table(genomeFile)
  classes <- c("factor","integer")
  for(i in 1:dim(data)[2]){class(data[,i]) <- classes[i]}
  colNames <- c("chr", "Length")
  colnames(data) <- colNames
  return(data)
}

readGapFile <- function(gapFile){
  data <- read.table(gapFile)
  ## ensure only dealing with BED3
  data <- data[,1:3]
  classes <- c("factor","integer", "integer")
  for(i in 1:dim(data)[2]){class(data[,i]) <- classes[i]}
  colNames <- c("chr", "start", "end")
  colnames(data) <- colNames
  return(data)
}

readNarrowPeakPlusNucBed <- function(narrowPeakPlusNucBedInfoFile){
  #Read Data in
  data <- read.table(narrowPeakPlusNucBedInfoFile)
  
  #Ensure class of each column is as it should be
  classes <- c("factor","integer","integer","factor","integer","factor","numeric","numeric","numeric","integer","numeric","numeric","integer","integer","integer","integer","integer","integer","integer")
  for(i in 1:dim(data)[2]){class(data[,i]) <- classes[i]}  
  
  #Name columns according to narrowPeak column names
  colNames <- c("chr", "start", "end", "name", "score", "strand", "FE", "pValScore", "qValScore", "relSummit", "percentAT", "percentGC", "A", "C", "G", "T", "N", "Other", "seqLen")
  colnames(data) <- colNames
  
  # Output data
  return(data)
}

writeBed <- function(bedObject,fileName){
  write.table(bedObject, file=fileName, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
}

readBed3 <- function(BED3File){
  #Read Data in
  data <- read.table(BED3File)
  
  #Name columns according to narrowPeak column names
  colNames <- c("chr", "start", "end")
  colnames(data) <- colNames
  
  # Output data
  return(data)
}

readBdg <- function(BdgFile){
  data <- read.table(BdgFile)
  colnames(data) <- c("chr", "start", "end", "score")
  return(data)
}
readBed3PlusNucBed <- function(BED3plusNucBedInfoFile){
  #Read Data in
  data <- read.table(BED3plusNucBedInfoFile)
  
  #Name columns according to narrowPeak column names
  colNames <- c("chr", "start", "end", "percentAT", "percentGC", "A", "C", "G", "T", "N", "Other", "seqLen")
  colnames(data) <- colNames
  
  # Output data
  return(data)
}

readBed5 <- function(BED5){
  #Read Data in
  data <- read.table(BED5)
  
  #Name columns according to narrowPeak column names
  colNames <- c("chr", "start", "end", "name", "score")
  colnames(data) <- colNames
  
  # Output data
  return(data)
}

readBed5PlusNucBed <- function(BED5plusNucBedInfoFile){
  #Read Data in
  data <- read.table(BED5plusNucBedInfoFile)
  
  #Name columns according to narrowPeak column names
  colNames <- c("chr", "start", "end", "name", "score", "percentAT", "percentGC", "A", "C", "G", "T", "N", "Other", "seqLen")
  colnames(data) <- colNames
  
  # Output data
  return(data)
}

readFilePaths <- function(filePathsTxtFile){
  ## reads in a txt file that has paths to files of interest (one per line)
  ## returns paths as a character vector
  files <- read.table(filePathsTxtFile)
  files <- as.matrix(files)
  files <- as.character(files)
  return(files)
}

writeFilePaths <- function(allPaths, outputFilename){
  ##allPaths is a vector of path strings --> c("pathTofIle1", "pathToFile2")
  ## writes a txt file with 1 file path per line
  allPaths <- data.frame(paths=allPaths)
  write.table(x=allPaths, file=outputFilename, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
}

histFE <- function(narrowPeakVar, ...){
  hist(narrowPeakVar$FE, ...)
}

histpVal <- function(narrowPeakVar, ...){
  hist(narrowPeakVar$pValScore, ...)
}

histqVal <- function(narrowPeakVar, ...){
  hist(narrowPeakVar$qValScore, ...)
}

histPeakLengths <- function(narrowPeakVar, ...){
  hist(narrowPeakVar$end-narrowPeakVar$start, ...)
}

mappableChromosomeSizes <- function(genomeFileVar, gapFileVar){
  ## To get the genomeFileVar and gapFileVar do:
  ## genomeFileVar <- readGenomeFile(genomeFile)
  ## gapFileVar <- readGapFile(gapFile)
  ## To get mappable genome size, sum the lengths of output
  ##---------------------------------------------
  ## initialize variable that will store sum of gap lengths for each chromosome
  gapLen <- c()
  ## Get sum of gap lengths along each chr
  for(chr in genomeFileVar$chr){
    ## 
    gapLen <- c(gapLen, sum(as.numeric(gapFileVar[gapFileVar$chr == chr, ][,3]) - as.numeric(gapFileVar[gapFileVar$chr == chr, ][,2])))
  }
  ## Subtract gap lengths for each chr
  mapLen <- genomeFileVar$Length - gapLen
  ## report mappable chromosome sizes
  return(mapLen)
}

mappableGenomeSize <- function(genomeFileVar, gapFileVar){
  return(sum(mappableChromosomeSizes(genomeFileVar, gapFileVar)))
}

numPeaksByChrTable <- function(narrowPeakVar, genomeFileVar, gapFileVar, byLength=TRUE, byMapLen=FALSE){
  ## Returns table with cols: chr len maplen numPeaks
  if(byLength && byMapLen){return("byLength and byMapLen cannot both be TRUE at same time.")}
  ##Calculate mappable lengths by subtracting gap lengths from total chr length
     #gapLen <- c()
     #for(chr in genomeFileVar$chr){gapLen <- c(gapLen, sum(as.numeric(gapFileVar[gapFileVar$chr == chr, ][,3]) - as.numeric(gapFileVar[gapFileVar$chr == chr, ][,2])))}
     #mapLen <- genomeFileVar$Length - gapLen
  ## 9/3/13 replaced above code with mappableChromosomeSizes function -- keep until test -- need to make sure I didnt break this fxn (not tested yet)
  mapLen <- mappableChromosomeSizes(genomeFileVar, gapFileVar)
  
  ##Obtain num peaks for each chr
  numPeaksByChr <- table(narrowPeakVar$chr)
  numPeaks <- c()
  for(chr in genomeFileVar$chr){numPeaks <- c(numPeaks, numPeaksByChr[chr])}
  
  #Glue mapLen and numPeak info to genome file info
  numPeaksByChrTable <- cbind(genomeFileVar, mapLen, numPeaks)
  colNames <- c("chr", "Length", "mappableLength", "numPeaks")
  colnames(numPeaksByChrTable) <- colNames

  if(byLength){
    lengthOrder <- order(numPeaksByChrTable$Length, decreasing=TRUE)
    return(numPeaksByChrTable[lengthOrder,])
  }
  else if(byMapLen){
    lengthOrder <- order(numPeaksByChrTable$mappableLength, decreasing=TRUE)
    return(numPeaksByChrTable[lengthOrder,])
  }
  else {return(numPeaksByChrTable)}
  
}

plotNumPeaksByChr <- function(numPeaksByChrTable, lenScale=1, peakScale=1, sideways=FALSE, mapLen=TRUE, ...){
  ##where lenScale is number to scale lengths by (e.g. 10+e6) and peakScale is that to scale numPeaks (e.g. 10e+4)
  ## mapLen needs to be true or False
    if(mapLen){lenVar <- numPeaksByChrTable$mappableLength/lenScale}
    else{lenVar <- numPeaksByChrTable$Length/lenScale}
    numPeaksByChrTable <- numPeaksByChrTable[complete.cases(numPeaksByChrTable), ]
    heightmax <- max(max(lenVar),max(numPeaksByChrTable$numPeaks/peakScale))
    if(sideways){
     barplot(rbind(lenVar, numPeaksByChrTable$numPeaks/peakScale),names.arg=numPeaksByChrTable$chr, beside=TRUE, las=1, xlim=c(0,heightmax), horiz=TRUE, ...)
    }
    else{
      barplot(rbind(lenVar, numPeaksByChrTable$numPeaks/peakScale),names.arg=numPeaksByChrTable$chr, las=2, beside=TRUE, ylim=c(0,heightmax), ...)
    }

 ##TODO -- make it possible to plot in order of length
  ## TODO -- plot chr as %longest (chr1 = 100%) and numFeatures as %max (order the bars by % longest or karyotype though)
}

interPeakDistances <- function(narrowPeakVar, genomeFileVar, gapFileVar, narrowPeakFormat=TRUE){
  ## any BED format var -- use narrowPeakFormat=TRUE if it truly is a narrowPeak (default)
  ## TO FINISH -- number distances seems right -- has some NAs though
  distances <- c()
  for(chr in genomeFileVar$chr){
    if(narrowPeakFormat){peakSummits <- narrowPeakVar$start[narrowPeakVar$chr == chr] + narrowPeakVar$relSummit[narrowPeakVar$chr == chr]}
    else {peakSummits <- round((narrowPeakVar$start+narrowPeakVar$end-1)/2)}
    gaps <- gapFileVar[gapFileVar$chr == chr, ]
    numGaps <- dim(gaps)[1]
    ##base -- i=1
    i <- 1
    gap <- gaps$gapStart[i]
    qualifiedSummits <- peakSummits[peakSummits < gap]
    numSummits <- length(qualifiedSummits)
    if(is.integer(numSummits) & numSummits > 1){ ##what is this? 8/11/2015
      distances <- c(distances, qualifiedSummits[2:(numSummits)]-qualifiedSummits[1:(numSummits-1)])
    }
    prevGap <- gaps$gapEnd[i]
    
    ##iter i=2:n
    for(i in 2:numGaps){
      gap <- gaps$gapStart[i]
      qualifiedSummits <- peakSummits[peakSummits > prevGap & peakSummits < gap]
      numSummits <- length(qualifiedSummits)
      if(is.integer(numSummits) & numSummits > 1){
        distances <- c(distances, qualifiedSummits[2:(numSummits)]-qualifiedSummits[1:(numSummits-1)])
      }
      prevGap <- gaps$gapEnd[i]
    }
   ## anything after last gap...?
    qualifiedSummits <- peakSummits[peakSummits > prevGap]
    numSummits <- length(qualifiedSummits)
    if(is.integer(numSummits) & numSummits > 1){
      distances <- c(distances, qualifiedSummits[2:(numSummits)]-qualifiedSummits[1:(numSummits-1)])
    }
  }
  return(distances)
}

##TODO --- chr st end distnace
# distanceToNextPeak <- function(narrowPeakVar, genomeFileVar, gapFileVar){
#   peakSummits <- narrowPeakVar$start[narrowPeakVar$chr == chr] + narrowPeakVar$relSummit[narrowPeakVar$chr == chr]
#   gaps <- gapFileVar[gapFileVar$chr == chr, ]
#   numGaps <- dim(gaps)[1]
#   ##base -- i=1
#   i <- 1
#   gap <- gaps$gapStart[i]
#   qualifiedSummits <- peakSummits[peakSummits < gap]
#   numSummits <- length(qualifiedSummits)
#   if(is.integer(numSummits) & numSummits > 1){
#     distances <- c(distances, qualifiedSummits[2:(numSummits)]-qualifiedSummits[1:(numSummits-1)])
#   }
#   prevGap <- gaps$gapEnd[i]
#   
#   ##iter i=2:n
#   for(i in 2:numGaps){
#     gap <- gaps$gapStart[i]
#     qualifiedSummits <- peakSummits[peakSummits > prevGap & peakSummits < gap]
#     numSummits <- length(qualifiedSummits)
#     if(is.integer(numSummits) & numSummits > 1){
#       distances <- c(distances, qualifiedSummits[2:(numSummits)]-qualifiedSummits[1:(numSummits-1)])
#     }
#     prevGap <- gaps$gapEnd[i]
#   }
#   ## anything after last gap...?
#   qualifiedSummits <- peakSummits[peakSummits > prevGap]
#   numSummits <- length(qualifiedSummits)
#   if(is.integer(numSummits) & numSummits > 1){
#     distances <- c(distances, qualifiedSummits[2:(numSummits)]-qualifiedSummits[1:(numSummits-1)])
#   }
# }
# return(distances)
# }







############Signal Around Point (e.g. summit or center) Analysis #####################
##TODO: Break up signalAroundPoints into separate functions (at least 1 fxn per step)
## This way if you already have the binned window files, for example, you can do the analysis from that step on...
##

signalAroundPoints <- function(pointsFile, bigWigFile, genomeFile, outputPrefix, gapFile, flankSize=2e+3, binSize=100, elimExtendedPointsThatOverlapGaps=TRUE, writeOut=TRUE){
  ## Assumptions:
  ## pointsFile is a BED3 file of summit locations, center locations, or any other collection of 1 bp points
  ## bigWigFile is the bigWig file with the signal track in it
  ## genome file (e.g. hg19.genome has 2 col: chr length)
  ##
  ## Paramters:
  ## flankSize is how much to add and subtract from each point to define boundaries: default=2kb
  ## binSize must be a factor of flanksize (can evenly divide it). It is used to partition the flanks into binSize sized bins.
  ## --> the average signal over these bins is extracted from bigWig file
  ##
  ## Output: returns a 'rowPile': a stack of rows (bins around points) where each column (bin) represents same distance away from point 
  ## outputPrefix -- should be a string (character)
  ## writeOut -- whether or not to write rowPile to a file or return to R
  ##          -- default is TRUE -- write it out to file (b/c bigWigAverageOverBed requires too much memory for my computer -- so I will do this on Oscar and transfer output back)
  ##
  ## Dependencies:
  ## This script depends on the following being in system environment:
  ##    - BEDtools
  ##    - kent utilities bigWigAverageOverBed
  ##    - readBed3 from this suite
  ##    - readBed5 from this suite
  ##
  ## Algorithm:
  ## -. Names each point (line) "point1" to "pointnumPoints" in column4
  ## -. Extends point by flankSize bp in both directions
  ## -. Uses bedtools makewindows to divide pointWindows up into binSize bins (keeping pointNames on each bin)
  ## -. Uses bigWigAverageOverBed to get average signal over each entry (cleans up unwanted files here)
  ## -. Takes in the signalScores BED file and creates a vector of 0s of numBinsPerPointWindow length
  ## -. For each pointName, it:
  ##              extracts entries 
  ##              ensures the correct number are present (if not, it moves on)
  ##              ensures entries are in correct order
  ##              rbinds the vector to pile of rowVectors
  ## -. returns pile of rowVectors (matrix)
  ## 
  ## -- One can then do a lot with the matrix:
  ##        - Get averages, SDs, SEMs over each position
  ##        - create lineplot of average and SD lines.
  ##        - do a line plot with confidence intervals, etc
  ########################################################################################
  ## Step 1: Read in points file AND extend flanks
  extendedPoints <- extendPoints(pointsFile, flankSize)
  
  ##Step 2: Filtering. Removal of points that enxtended beyond chromosome limits  (else BEDtools quits) -- optional rmeoval of those that overlap gaps (default is TRUE)
  msg <- paste("Eliminating points that fall off chromosome due to adding/subtracting ", flankSize, " from point (i.e. start coord. < 0 OR end coord > chrLength):", sep="")
  print(msg)
  ## ensure extendedPoint start>=0  (else BEDtools quits)
  extendedPoints <- removeRegionsThatFallOffChrStart(extendedPoints)
  ## ensure extendedPoint end <= chrLen
  extendedPoints <- removeRegionsThatFallOffChrEnd(extendedPoints, genomeFile)
  ##Remove points that overlap gaps (NNN signal can add noise to analysis -- if no signal (or wrong signal) over the NNNs)
  if(elimExtendedPointsThatOverlapGaps){
    extendedPoints <- removeRegionsThatOverlapGaps(extendedPoints, gapFile)
  }
  
  ### Step 3: giving names to each point regions
  print("Adding names to each point region that will be used in analysis....")
  extendedPoints <- nameBedRegions(BEDregions=extendedPoints, name="point")
  
  ## Step 4: Write out the pointWindows file
  pointWindowOutFile <- paste("pointWindows_", outputPrefix, ".bed", sep="")
  print(paste("Writing out ", pointWindowOutFile, "...", sep=""))
  writeBed(bedObject=extendedPoints, fileName=pointWindowOutFile)

  ##Step 5: Divide up Regions into bins, preserving pointName information on each associated bin
  binnedOutputPrefix <- paste0("pointBinnedWindows_", outputPrefix)
  pointBinnedWindowInFile <- getBinsInsideBedRegions(pointWindowOutFile, outputPrefix=binnedOutputPrefix, binSize)
      ## uses default pathToBEDtools  

  ##Step 4: Generate average signal over each bin in pointBinnedWindowInFile -- get path to its output
  bigWigAverageOverBedOutput <- bigWigAverageSignalOverBed(bigWigFile, pointBinnedWindowInFile, outputPrefix) 
        ## uses default pathToBigWigAverageOverBed

  ##Step 5: Read in the output of bigWigAverageOverBed
  ## out.tab file from bigWigAverageOverBed has following 6 columns: name, size of bin (bp), num bp in bin covered, sum, mean of entire bin, mean of only covered bases
  ## Will be interested in mean (col 5) over all bases (not just covered ones) -- alt. Sum (col4) over all bases 
  ## bins were named pointi_j 
  ## put j in its own column (binNumber)
  ## -- shortened all pointnames to pointi (by removing _j)
  ## Step 6: go through each point name group, ensure that:
  ##   - each has the correct number of bins (sum of binNums == sum(1:numBinsPerPoint) -- or length() == numBinsPerPoint)
  ##   - that the bins are in the correct order (order by binNum)
  ## Add each vector of signal values to row pile
    ## continue if correct number of bins present
      ## extract only data associated with current point
      ## ensure bins are in order from 1:numPointsPerBin
      ## Add pointOfInterest$meanBin vector to row pile
  rowPile <- makeRowPile(bigWigAverageOverBedOutput, flankSize, binSize)
  
  ## step 7: return the rowPile or write it out
  if(writeOut){
    ## Note a bed file, but writeBed has all the correct paramters
    rowPileOut <- paste0("rowPile_", outputPrefix, ".txt")
    writeBed(bedObject=rowPile, fileName=rowPileOut)
    return()
  }
  else {return(rowPile)}
  
  ## step8: To do with returned data
  ## Can take SUM or MEAN of columns
  ## plot Mean line with SD lines above and below... or confidence intervals
}



extendPoints <- function(pointsFile, flankSize){
  ##Read in points file, extend flanks
  ## Returns extendedPoints
  ### Read in
  print("Processing points file....extending points to flankSize in both directions")
  points <- readBed3(pointsFile)
  ### extend flanks
  points$end <- points$start + flankSize ##need to change end before change start (since it uses start)
  ## start+flank is intentional here (rather than end + flank)
  points$start <- points$start - flankSize
  return(points) ## returns extendedPoints
}

removeRegionsThatFallOffChrStart <- function(BEDregions){
  ## For signalAroundPoint analysis, BEDregions = extendedPoints
  BEDregionsWithStartsLessThanZero <- BEDregions[BEDregions$start < 0, ]
  numBEDregionsWithStartsLessThanZero <- length(BEDregionsWithStartsLessThanZero$start)
  msg <- paste("The following ", numBEDregionsWithStartsLessThanZero, " regions were eliminated because they extend beyond chromosome start (region start coord. less than 0):", sep="")
  print(msg)
  print(BEDregionsWithStartsLessThanZero)
  BEDregions <- BEDregions[BEDregions$start >= 0, ]
  return(BEDregions)
}

removeRegionsThatFallOffChrEnd <- function(BEDregions, genomeFile){
  ## For signalAroundPoint analysis, BEDregions = extendedPoints
  msg <- paste("The following regions were eliminated because they extended beyond chromosome end (region end coord. greater than chrLength):", sep="")
  print(msg)
  genomeInfo <- readGenomeFile(genomeFile)
  for(chr in names(table(BEDregions$chr))){
    ## Get chr length
    chrLength <- genomeInfo[genomeInfo$chr == chr, ]$Length
    ## extract entries for chr of interest
    ## take note of which and how many fall of chr
    BEDregionsThatFallOffChr <- BEDregions[BEDregions$chr == chr & BEDregions$end > chrLength, ]
    numBEDregionsThatFallOffChr <- length(BEDregionsThatFallOffChr$end)
    ## if any fall off, report it and modify set of points to exclude them
    if(numBEDregionsThatFallOffChr > 0){
      msg <- paste(numBEDregionsThatFallOffChr, " region(s) on ", chr, " -- chrLength = ", chrLength, sep="")
      print(msg)
      print(BEDregionsThatFallOffChr)
      BEDregionsToExtract <- !(BEDregions$chr == chr & BEDregions$end > chrLength)
      BEDregions <- BEDregions[BEDregionsToExtract, ]
    } ## END IF numPointsThatFallOffChr > 0
  } ## END For chr in....
  return(BEDregions)
}

removeRegionsThatOverlapGaps <- function(BEDregions, gapFile){
  ## For signalAroundPoint analysis, BEDregions = extendedPoints
  msg <- paste("The following regions were eliminated because they overlapped gaps:", sep="")
  print(msg)
  gapInfo <- readGapFile(gapFile)
  numRegionsBeforeFiltering <- dim(BEDregions)[1]
  for(chr in names(table(BEDregions$chr))){
    ## intialize variables
    numOverlapsOnThisChr <- 0
    theOverlapsOnThisChrAre <- c()
    ## Check each gap on that chr
    ## extract entries for chr of interest
    gapsOnChr <- gapInfo[gapInfo$chr == chr, ]
    numGaps <- dim(gapsOnChr)[1]
    for(i in 1:numGaps){
      gap_i <- gapsOnChr[i,]
      ## Need to see if start or end falls within gap (partial overlap of gap), or if entire gap falls within BED region (full overlap over gap)
      ## The simpler way to think about it is if a gap start or end falls in a bed region, the bed region overlaps the gap (whether partially or fully)
      overlapGapLogic <- BEDregions$chr == chr & ((BEDregions$start < gap_i$start & BEDregions$end > gap_i$start) | (BEDregions$start < gap_i$end & BEDregions$end > gap_i$end))
      BEDregionsThatOverlapGap <- BEDregions[overlapGapLogic, ]
      numBEDregionsThatOverlapGap <- length(BEDregionsThatOverlapGap$end)
      ## if any fall off, report it and modify set of points to exclude them
      if(numBEDregionsThatOverlapGap > 0){
        ##update count for this chr
        numOverlapsOnThisChr <- numOverlapsOnThisChr + numBEDregionsThatOverlapGap
        ##update regions list for this chr
        theOverlapsOnThisChrAre <- rbind(theOverlapsOnThisChrAre, BEDregionsThatOverlapGap)
        ## update the regions to keep by the negatied logic of which ones are violaters 
        BEDregionsToExtract <- !(overlapGapLogic)
        BEDregions <- BEDregions[BEDregionsToExtract, ]
      } ## END IF numPoints...
    } ## END for gap_i in 1:numGaps...
    ## Report how many violating regions on this chr
    msg <- paste(numOverlapsOnThisChr, " region(s) on ", chr, sep="")
    print(msg)
    ## Report the violating regions for this chr
    print(theOverlapsOnThisChrAre)
  } ## END For chr in....
  numRegionsAfterFiltering <- dim(BEDregions)[1]
  print("Finished filtering out regions that overlapped gaps....")
  print(paste0("Started with ",   numRegionsBeforeFiltering, " regions."))
  print(paste0("Filtered out ", numRegionsBeforeFiltering-numRegionsAfterFiltering, " regions that overlapped gaps."))
  print(paste0("Ended with ", numRegionsAfterFiltering, " regions."))
  ## return the BEDregions with all regions that overlapped gaps removed
  return(BEDregions)
}



nameBedRegions <- function(BEDregions, name){
  ## Adds name and number to each region and adds this as a new column (last postion) to the input
  ## So if a BED3, the names will be in column 4 (if BED5, then column 6, etc)
  ## The name can be any string. All entries get same name with row_number appended with no spaces or underscores
  ## Thus if you want an underscore (name# vs. name_#), simply add it in the name string.
  #
  ## For signalAroundPoint analysis, BEDregions = extendedPoints; name="point"
  #
  numBEDregions <- dim(BEDregions)[1]
  bedNames <- c()
  for(i in 1:numBEDregions){bedNames <- c(bedNames, paste(name, i, sep=""))}
  BEDregions <- cbind(BEDregions, bedNames)
  return(BEDregions)
}

getBinsInsideBedRegions <- function(pathToBedFile, outputPrefix, binSize, pathToBEDtools="~/searchPaths/bedtools/current/bedtools"){
  ## Where:
  ##    pathToBedFile tells R where to find bed file (in which, the regions will be partitioned into bins)
  ## Where:
  ##    pathToBEDtools tells R where to find bedtools
  ##      Default on my computer: "~/searchPaths/bedtools/current/bedtools"
  ##      Default on Oscar: "~/software/bedtools/current_bin/bedtools"
  ##
  ## Uses bedtools makewindows (dependency) to generate the file in the working directory
  ## Returns the path to output file
  #
  ## In the signalAroundPoints analysis, the input BEDfile (pathToBedFile) to use is: pointWindowOutFile
  ##    And the output is the: pointBinnedWindowInFile  ## so do pointBinnedWindowInFile <- getBinsInsideBedRegions(pathToBedFile=pointWindowOutFile,...)
  ## TODO
  ## For now -- it should be stated that this is made especially for signalAroundPointAnalysis
  ##   Therefore it is not yet flexible... but I should actually turn this into a completely flexible makewindows R wrapper
  ## Make sure to update signalAroundPoints() if change parameters and/or name
  pathToBinnedBedFile <- paste(outputPrefix, ".bed", sep="")
  print(paste("Partitioning BED regions into ", binSize, " bp bins. Writing to ", pathToBinnedBedFile, "....", sep=""))
  bedtoolWindowCmd <- paste(pathToBEDtools," makewindows -b ", pathToBedFile, " -w ", binSize, " -i srcwinnum > ", pathToBinnedBedFile)
  system(bedtoolWindowCmd)
  return(pathToBinnedBedFile)
}

bigWigAverageSignalOverBed <- function(bigWigFile, pointBinnedWindowInFile, outputPrefix, pathToBigWigAverageOverBed="~/searchPaths/kentResourceTree/bigWigAverageOverBed"){
  ## Where:
  ##     bigWigFile is the signal track of interest
  ##     pointBinnedWindowInFile is a BED file made by breaking up the BED regions in a "pointWindowOutFile 
  ##          (generated by extendPoints)" into bins (using breakPointWindowsIntoBins)
  ##              ...where pointWindowOutFile is the BED result of extending points to flankSize in each direction
  ##
  ## Where pathToBigWigAverageOverBed
  ##    Default on my computer: "~/searchPaths/kentResourceTree/bigWigAverageOverBed"
  ##    Default on Oscar: "~/software/kentTools/bigWigAverageOverBed"
  ## This function will use bigWigAverageOverBed (dependency) to generate the out.tab file
  ##    That file will be in the working directory
  ##    The function returns the path to it
  ##
  ## Get average signal over each bin in pointBinnedWindowInFile
  bigWigAverageOverBedOutput <- paste("pointBinnedWindows_", outputPrefix, ".tab", sep="")
  print(paste("Performing bigWigAverageOverBed step (getting info for all bins). Writing to ", bigWigAverageOverBedOutput, "....", sep=""))  
  bigWigAverageOverBed <- paste(pathToBigWigAverageOverBed, " ", bigWigFile, " ", pointBinnedWindowInFile, " ", bigWigAverageOverBedOutput, sep="")
  system(bigWigAverageOverBed)
  
  ## Return location to output file
  return(bigWigAverageOverBedOutput)
}
  
makeRowPile <- function(bigWigAverageOverBedOutput, flankSize, binSize){
  ##Read in the output of bigWigAverageOverBed
  ## out.tab file from bigWigAverageOverBed has following 6 columns: name, size of bin (bp), num bp in bin covered, sum, mean of entire bin, mean of only covered bases
  ## Will be interested in mean (col 5) over all bases (not just covered ones) -- alt. Sum (col4) over all bases 
  binnedSignalData <- read.table(bigWigAverageOverBedOutput)
  colnames(binnedSignalData) <- c("name", "size", "covered", "sum", "meanBin", "meanCovered")  
  ## bins were named pointi_j 
  ## put j in its own column (binNumber)
  binNumColumn <- as.integer(sub("point.+_", "", binnedSignalData$name))
  binnedSignalData <- cbind(binnedSignalData, binNumColumn)
  colnames(binnedSignalData) <- c("name", "size", "covered", "sum", "meanBin", "meanCovered", "binNum")
  ## -- shorted all pointnames to pointi (by removing _j)
  binnedSignalData$name <- sub("_.+", "", binnedSignalData$name)
  ## calculate number of bins each point should have
  numBinsPerPoint <- flankSize*2/binSize
  ## make matrix, where each row is the signal scores from left to right on each point's locus -- i.e. a pile of row vectors
  rowPile <- c()
  ## I think that below should be: dim(binnedSignalData)[1]/binSize or... dim(binnedSignalData)[1]/numBinsPerPoint
  for(i in 1:dim(binnedSignalData)[1]){
    pointName <- paste("point", i, sep="")
    ## continue if correct number of bins present
    print(pointName)
    if((length(binnedSignalData$name[binnedSignalData$name == pointName]) == numBinsPerPoint)){
      ## extract only data associated with current point
      pointOfInterest <- binnedSignalData[binnedSignalData$name == pointName, ]
      ## ensure bins are in order from 1:numPointsPerBin
      pointOfInterest <- pointOfInterest[order(pointOfInterest$binNum), ]
      ## Add pointOfInterest$meanBin vector to row pile
      rowPile <- rbind(rowPile, pointOfInterest$meanBin)
    }
  }
  return(rowPile)
}


#### overlap Analysis ############################################################################
basicBinomialTestForOverlapSignificance <- function(obsNumOverlap, numPeaksA, meanPeakWidthA, numPeaksB, meanPeakWidthB, genomeSize=2835679040, numConnectedComponents=257, minOverlap=1){
  totalPositions <- genomeSize - (numConnectedComponents*(meanPeakWidthA-1))
  numSuccessfulPositions <- (meanPeakWidthA + meanPeakWidthB - minOverlap)*numPeaksB
  numSuccessfulPositions <- min(numSuccessfulPositions, totalPositions)
  obsProportion <- obsNumOverlap/numPeaksA
  expProportion <- numSuccessfulPositions/totalPositions
#   expNumPeaks <- expProportion*numPeaksA
  pValue <- pbinom(q=obsNumOverlap, size=numPeaksA, prob=expProportion, lower=FALSE) 
  return(list(obsProp=obsProportion, expProp=expProportion, pvalue=pValue, obsToExpRatio=obsProportion/expProportion))
}

basicBinomialTestWithFilePaths <- function(obsNumOverlap, bedFileA, bedFileB, genomeSize=2835679040, numConnectedComponents=257, minOverlap=1){
  bedA <- read.table(bedFileA)[,2:3]; colnames(bedA) <- c("start", "end")
  bedB <- read.table(bedFileB)[,2:3]; colnames(bedB) <- c("start", "end")
  numPeaksA <- nrow(bedA); meanPeakWidthA <- mean(bedA$end-bedA$start)
  numPeaksB <- nrow(bedB); meanPeakWidthB <- mean(bedB$end-bedB$start)  
  print(c(numPeaksA, meanPeakWidthA, numPeaksB, meanPeakWidthB))
  basicBinomialTestForOverlapSignificance(obsNumOverlap=obsNumOverlap, numPeaksA=numPeaksA, meanPeakWidthA=meanPeakWidthA, numPeaksB=numPeaksB, meanPeakWidthB=meanPeakWidthB, genomeSize=genomeSize, numConnectedComponents=numConnectedComponents, minOverlap=minOverlap)
}

readNamesFile <- function(namesTxtFile){
  namesVector <- read.table(namesTxtFile)
  return(namesVector)
}

makeNamesDataFrame <- function(namesVector){
  ## where namesVector is a character vector -- e.g. c("A", "B", "C")
  return(as.data.frame(namesVector))
}

readOverlapTable <- function(overlapTableTxtFile, rowNames, colNames){
  ## names variables are Nx1 data frames (as output by readNamesFile())  
  overlapTable <- read.table(overlapTableTxtFile)
  overlapTable <- sapply(overlapTable, as.numeric)
  rownames(overlapTable) <- rowNames[,1]
  colnames(overlapTable) <- colNames[,1]
  return(overlapTable)
}



binomialOverLapPval <- function(obsNumPeaksOverlapTable, bedFileListFile, genomeSize=2835679040, numConnectedComponents=257, minOverlap=1){
  ## where:
  ##    obsNumPeaksOverlapTable is output of readOverlapTable() -- already processed with row/col names
  ##    --Assumes square table where rows == cols
  ##    --Assumes bedFileList is in same order as rows and cols -- i.e. order(rows) == order(cols) == order(bedfilelist)  
  ## ensure numeric parameters are numeric
  ## default genome size is female nuclear mappable genome size calculated for Hg19 by deleting total length of gaps from total length genome
  ##         nuclear female mappable genome size is 2835679040; 
  ##         Male nuclear+mito = 2861349177
  ## default numConnectedComponents estimated by how many regions are separated by gaps -- i.e. use complementBed on the hg19 gaps bed file
  ##         nuclear female numConComp = 257
  ##         nuclear+mito male numConComp = 275
  genomeSize <- as.numeric(genomeSize)
  numComponents <- as.numeric(numConnectedComponents)
  minOverlap <- as.numeric(minOverlap)
  
  ## Access BED files used
  BEDfiles <- readFilePaths(bedFileListFile)
  
  ## GET NUMPEAKS AND AVG PEAK WIDTHS FROM EACH BED FILE
  print("Getting numPeak and avgPeakWidth information from all BED files.....")
  numPeaks <- c()
  meanPeakWidth <- c()
  for(bedfile in BEDfiles){
    BED <- read.table(bedfile)
    ENDS <- as.integer(BED[,3])
    STARTS <- as.integer(BED[,2])
    WIDTHS <- ENDS-STARTS
    MEANWIDTH <- mean(WIDTHS)
    NUMPEAKS <- length(WIDTHS)
    meanPeakWidth <- c(meanPeakWidth, MEANWIDTH)
    numPeaks <- c(numPeaks, NUMPEAKS)
  }
  
  ## INITIALIZING OUTPUT TABLES
  ##obsNumPeaksOverlapTable is used to initilize new tables. This will be overwritten. It just has correct dimensions.
  obsProportionTable <- obsNumPeaksOverlapTable
  expProportionTable <- obsNumPeaksOverlapTable ## a.k.a. table of success probabilities
  expNumPeaksTable <- obsNumPeaksOverlapTable
  pValueTable <- obsNumPeaksOverlapTable
  setSizeRatioTable <- obsNumPeaksOverlapTable ## This will be |A|/|B| 
  
  ## FILLING IN OUTPUT TABLES
  print("Filling out output tables....")
  for (A in 1:length(meanPeakWidth)){
    for (B in 1:length(numPeaks)){
      ## calculate total positions
      totalPositions <- genomeSize - (numComponents*(meanPeakWidth[A]-1)) ## this was numPeaks[A] until june 24, 2014 when I changed it back to what it should have been --> meanWidthA-1 (got it right in matlab script - well I left -1 off since it doesnt make a diff with the genome size)
      ## calculate number of successful positions allowing redundancy (removed next step)
      numSuccessfulPositions <- (meanPeakWidth[A] + meanPeakWidth[B] - minOverlap)*numPeaks[B]
      ## remove redundantly counted successful positions -- where redundancy begins once all possible positions are considered
      ##  i.e. ensure that there cannot be more successful positions than total positions
      numSuccessfulPositions <- min(numSuccessfulPositions, totalPositions)
      ## calculate observed proportion -- numOverlaps/numPeaks
      obsProportionTable[A,B] <- obsNumPeaksOverlapTable[A,B]/numPeaks[A]
      ## calculate expected proportion - i.e. Pr(successfulOverlap)
      expProportionTable[A,B] <- numSuccessfulPositions/totalPositions
      ## calculate expected number of overlaps
      expNumPeaksTable[A,B] <- expProportionTable[A,B]*numPeaks[A]
      ## calculate pvalue for observed number of A peaks that overlap B peaks given num peaks in A (trials) and pr(success)
      pValueTable[A,B] <- pbinom(obsNumPeaksOverlapTable[A,B], numPeaks[A], expProportionTable[A,B], lower=FALSE) ##Alt. Hyp is more overlaps than expected by chance
      setSizeRatioTable[A,B] <- numPeaks[A]/numPeaks[B]
    }
  }
  allTables <- list(obsNumPeaksOverlapTable = obsNumPeaksOverlapTable, expNumPeaksTable = expNumPeaksTable, obsProportionTable = obsProportionTable, expProportionTable = expProportionTable, pValueTable = pValueTable, setSizeRatioTable = setSizeRatioTable)
  return(allTables)
}

writeAllTables <- function(allTables, outputFilename){
  i=0
  for(e in allTables){
    ## Add 1 to index i
    i = i+1
    ##Table name
    write.table(names(allTables)[i], file=outputFilename, append=TRUE, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
    ##Table
    write.table(e, file=outputFilename, append=TRUE, quote=FALSE, sep="\t", row.names=TRUE, col.names=NA)
    ## Space
    write.table(" ", file=outputFilename, append=TRUE, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  }
}

pairwiseBinomialInfo <- function(allTables){
  numRows <- dim(allTables$obsNumPeaksOverlapTable)[1]
  numCols <- dim(allTables$obsNumPeaksOverlapTable)[2]
  rowNames <- rownames(allTables$obsNumPeaksOverlapTable)
  colNames <- colnames(allTables$obsNumPeaksOverlapTable)
  pairwiseTable <- data.frame(fileA=c(), fileB=c(), expNum=c(), obsNum=c(), expProportion=c(), obsProportion=c(), pVal=c(), obsToExpRatio=c())
  #pairwiseTable <- data.frame()
  for(A in 1:numRows){
    for(B in 1:numCols){
      newRow <- data.frame(fileA=rowNames[A], fileB=colNames[B], expNum=allTables$expNumPeaksTable[A,B], obsNum=allTables$obsNumPeaksOverlapTable[A,B], expProportion=allTables$expProportionTable[A,B], obsProportion=allTables$obsProportionTable[A,B], pVal=allTables$pValueTable[A,B], obsToExpRatio=allTables$obsProportionTable[A,B]/allTables$expProportionTable[A,B], setSizeAtoBRatio=allTables$setSizeRatioTable[A,B], approxMaxProportion=min(1/allTables$setSizeRatioTable[A,B],1), percentOfMax=allTables$obsProportionTable[A,B]/min(1/allTables$setSizeRatioTable[A,B],1))
      pairwiseTable <- rbind(pairwiseTable, newRow)
      #return(pairwiseTable)
    }
  }
  return(pairwiseTable)
}

writePairwiseBinomialInfo <- function(pairwiseTable, outputFilename){
  write.table(pairwiseTable, file=outputFilename, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
}
  
plotExpVsObsOverlap <- function(pairwiseTable, noSelf=TRUE, legendX="topright", cex.names=0.7, ylim=c(0,1.3), yaxt="n", ...){
  ## If noSelf=FALSE, it will show the expected and observed of self vs. self (where observed is always 100%)
  numObservations <- dim(pairwiseTable)[1]
  if(noSelf){
    pairwiseTable <- pairwiseTable[pairwiseTable$fileA != pairwiseTable$fileB, ]
    numObservations <- dim(pairwiseTable)[1]
    }

  comparisons <- c()
  for(i in 1:numObservations){
    comparisons <- c(comparisons, paste0(pairwiseTable$fileA[i], " & ", pairwiseTable$fileB[i]))
  }
  barplot(rbind(pairwiseTable$expProportion,pairwiseTable$obsProportion), col=c("red","blue"), ylim=ylim, yaxt=yaxt, beside=TRUE, las=2, ylab="Proportion", names=comparisons, cex.names=0.6, horiz=F, ...)  
  legend(x=legendX, legend=c("Expected", "Observed"), fill=c("red","blue"), bty="n")
  axis(2, at=seq(0,1,0.1), labels=seq(0,1,0.1), las=1)
}

plotObsExpOverlapRatio <- function(pairwiseTable, noSelf=TRUE, cex.names=0.6, ...){
  numObservations <- dim(pairwiseTable)[1]
  if(noSelf){
    pairwiseTable <- pairwiseTable[pairwiseTable$fileA != pairwiseTable$fileB, ]
    numObservations <- dim(pairwiseTable)[1]
  }
  comparisons <- c()
  for(i in 1:numObservations){
    comparisons <- c(comparisons, paste0(pairwiseTable$fileA[i], " & ", pairwiseTable$fileB[i]))
  }
  barplot(pairwiseTable$obsToExpRatio, names=comparisons, las=2, cex.names=0.6, ylab="Observed:Expected", ...)
}

plotSetSizeRatios <- function(pairwiseTable, noSelf=TRUE, cex.names=0.7, ...){
  numObservations <- dim(pairwiseTable)[1]
  if(noSelf){
    pairwiseTable <- pairwiseTable[pairwiseTable$fileA != pairwiseTable$fileB, ]
    numObservations <- dim(pairwiseTable)[1]
  }
  comparisons <- c()
  for(i in 1:numObservations){
    comparisons <- c(comparisons, paste0(pairwiseTable$fileA[i], " & ", pairwiseTable$fileB[i]))
  }
  barplot(pairwiseTable$setSizeAtoBRatio, beside=TRUE, names=comparisons, las=2, cex.names=0.6, ylab="numPeaks(A):numPeaks(B)", ...)
}

plotSetSizeRatiosWithApproxMaxPercentMax <- function(pairwiseTable, noSelf=TRUE, ...){
  numObservations <- dim(pairwiseTable)[1]
  if(noSelf){
    pairwiseTable <- pairwiseTable[pairwiseTable$fileA != pairwiseTable$fileB, ]
    numObservations <- dim(pairwiseTable)[1]
  }
  comparisons <- c()
  for(i in 1:numObservations){
    comparisons <- c(comparisons, paste0(pairwiseTable$fileA[i], " & ", pairwiseTable$fileB[i]))
  }
  OneHundredPercent <- max(pairwiseTable$setSizeAtoBRatio)
  adjusted_approxMaxProportion <- OneHundredPercent*pairwiseTable$approxMaxProportion
  adjusted_percentOfMax <- OneHundredPercent*pairwiseTable$percentOfMax
  barplot(rbind(pairwiseTable$setSizeAtoBRatio, adjusted_approxMaxProportion, adjusted_percentOfMax), beside=TRUE, names=comparisons, cex.names=0.6, las=2, ylab="numPeaks(A):numPeaks(B)", ...)
  axis(4, at=OneHundredPercent*seq(0,1,0.1), labels=seq(0,1,0.1), las=1)
  ## TO FINISH --- needs legend... see plotExpVsObsOverlap for help
  return("TO FINISH: needs legend")
}

plotApproxMaxWithPercentMax <- function(pairwiseTable, noSelf=TRUE, legendX=1.3, ylim=c(0,1.3), col=c("light grey", "black"), ...){
  numObservations <- dim(pairwiseTable)[1]
  if(noSelf){
    pairwiseTable <- pairwiseTable[pairwiseTable$fileA != pairwiseTable$fileB, ]
    numObservations <- dim(pairwiseTable)[1]
  }
  comparisons <- c()
  for(i in 1:numObservations){
    comparisons <- c(comparisons, paste0(pairwiseTable$fileA[i], " & ", pairwiseTable$fileB[i]))
  }
  barplot(rbind(pairwiseTable$approxMaxProportion, pairwiseTable$percentOfMax), col=col, ylim=ylim, yaxt="n", beside=TRUE, names=comparisons, cex.names=0.6, las=2, ylab="Proportion", ...)
  legend(x=legendX, legend=c("approxMaxProportion", "proportionOfMax"), fill=col, bty="n")
  axis(2, at=seq(0,1,0.1), labels=seq(0,1,0.1), las=1)
  ## TO FINISH ---needs legend... see plotExpVsObsOverlap for help
  return("TO FINISH: needs legend")
}

plotlog10BinomialpVal <- function(pairwiseTable, noSelf=TRUE, cex.names=0.7, ...){
  numObservations <- dim(pairwiseTable)[1]
  if(noSelf){
    pairwiseTable <- pairwiseTable[pairwiseTable$fileA != pairwiseTable$fileB, ]
    numObservations <- dim(pairwiseTable)[1]
  }
  comparisons <- c()
  for(i in 1:numObservations){
    comparisons <- c(comparisons, paste0(pairwiseTable$fileA[i], " & ", pairwiseTable$fileB[i]))
  }
  ##pseudocount of 4.940656e-324 is added to all p-val (4.940656e-324 is ~smallest number R can handle)
  ### This makes all zeros cap off at -log10(p) = 323.3062
  ymax <- max(-log10(pairwiseTable$pVal+4.940656e-324))+49
  barplot(-log10(pairwiseTable$pVal+4.940656e-324), beside=TRUE, names=comparisons, ylim=c(0,ymax), las=2, cex.names=0.6, ylab="-log10(pValue)", ...)
}

getAvgOverDiagonalMatrix <- function(scoreMatrix){
  ## where 'scoreMatrix' is, for example, matrix of p-values (output of binomial overlap analysis)
  ## This was made to make pValue matrices (output by overlapTableBinomialPvalue.R) that are identical about the diagonal 
  ## by averaging the 2 pValues for A,B comparisons (A-->B, and B--A)
  ## The average pValues can then be used to make dendrograms and clustered heatmaps of p-value matrices
  ## This combines averages the pvalues of how well A^B AND B^A
  ## But it is nothing special. It is just averages row i and col j for anything.
  
  avgAB_BA <- scoreMatrix ## to give same dims
  for(A in 1:(dim(scoreMatrix)[1])){
    for(B in 1:(dim(scoreMatrix)[2])){
      avgAB_BA[A,B] <- mean(c(scoreMatrix[A,B], scoreMatrix[B,A])) #(avgAB_BA[A,B] + avgAB_BA[B,A])/2
    }
  }
  return(avgAB_BA)
}



matrixPlot <- function(scoreMatrix, title="Title", ylab="-log10(pValue) for the number of peaks in this set...",xlab="...that overlap peaks in this set",...){
  #return(scoreMatrix)
  levelplot(t(scoreMatrix)[, dim(scoreMatrix)[1]:1], pretty=TRUE, main=title, col.regions=heat.colors, ylab=ylab, xlab=xlab, ...)
}

matrixCluster <- function(scoreMatrix){
  ## *usually pValue matrix .... could do averageOverDiagonal of pvalues
  ## This was made to make clustered heatmaps of p-value matrices output by overlapTableBinomialPvalue.R
  ## But it is nothing special. It is just an R heatmap for anything.
  ##  
  dataDist <- dist(scoreMatrix)
  dataClust <- hclust(dataDist)
  plot(dataClust, xlab="xlab")
}


########## genome Feature Density Correlations Suite ##############################################################
###### These were made in another file some time in July 2013.
###### I want/need to update them to utilize the larger suite of functions where applicable as they have redundant code for now
## Some of the terms need to be generalized-- e.g. counts and percents are both just scores associated with bins
## A bash function (called bedGraphLoop) that automates the counts file production is usually used before this analysis
pathToBedGraphLoop <- "/Users/johnurban/searchPaths/shell_scripts_Mac/bedGraphLoop.sh"
bedGraphLoopWrapper <- function(){return()}
##TODO
## This may require that bedGraphLoop uses specific paths to BEDtools functions 

countMatrix <- function(countFiles){
  ##countFiles is a txt file of all files you want R to access for this
  ## The type of files it points to are counts (or some type of score such as GC%) for each bin across a genome
  ##    This usually is from the 4th column of a bedGraph type file
  ##can also be 'percentFiles.txt'
  files <- read.table(countFiles)
  files <- as.character(as.matrix(files))
  ##Create matrix of all counts (from all files)
  countMatrix <- c()
  for(i in 1:length(files)){
    data <- read.table(files[i]); 
    data <- as.numeric(as.matrix(data))
    countMatrix <- cbind(countMatrix, data)
  }
  colNames <- sub(".txt", "", files)
  colNames <- gsub("../","", colNames)
  colNames <- sub("counts_", "", colNames)
  colnames(countMatrix) <- colNames
  return(countMatrix)
}

corrMat_counts <- function(countMatrix, method){
  corMat <- matrix(nrow=dim(countMatrix)[2], ncol=dim(countMatrix)[2])
  for(i in 1:(dim(countMatrix)[2]-1)){
    for(j in (i+1):dim(countMatrix)[2]){
      feat1 <- countMatrix[,i]
      feat2 <- countMatrix[,j]
      corMat[i, j] <- cor(feat1, feat2, method=method)
      
    }
  }
  colRowNames <- colnames(countMatrix)
  colRowNames <- sub(".txt", "", colRowNames)
  colRowNames <- sub("counts_", "", colRowNames)
  colnames(corMat) <- colRowNames
  rownames(corMat) <- colRowNames
  return(corMat)    
}

corrMat_percentNT <- function(countMatrix, percentMatrix, method){
  corMat <- matrix(nrow=dim(percentMatrix)[2], ncol=dim(countMatrix)[2])
  for(i in 1:(dim(percentMatrix)[2])){
    percents <- percentMatrix[,i]
    for(j in 1:dim(countMatrix)[2]){
      features <- countMatrix[,j]
      corMat[i, j] <- cor(percents, features, method=method)
      
    }
  }
  colNames <- colnames(countMatrix)
  colNames <- sub(".txt", "", colNames)
  colRowNames <- sub("counts_", "", colNames)
  colnames(corMat) <- colRowNames
  rownames(corMat) <- colnames(percentMatrix)
  return(corMat)    
}

writeOutCorMat <- function(corMat, outFileName){
  write.table(corMat, file=outFileName, quote=FALSE, sep="\t", na="", row.names=TRUE, col.names=TRUE)
}

corLoopWrapper_allPairwiseSingleFile <- function(countFiles, outFileName, method){
  writeOutCorMat(corrMat_counts(countMatrix(countFiles), method=method), outFileName)
}

corLoopWrapper_allPairwiseTwoFiles <- function(countFiles, percentFiles, outFileName, method){
  writeOutCorMat(corrMat_percentNT(countMatrix(countFiles), countMatrix(percentFiles), method=method), outFileName)
}

corLoopWrapper_allPairwiseTwoFiles_transpose <- function(countFiles, percentFiles, outFileName, method){
  writeOutCorMat(t(corrMat_percentNT(countMatrix(countFiles), countMatrix(percentFiles), method=method)), outFileName)
}























###### Future or discontinued #######

intersectBedLoopWrapper <- function(filePathsTxtFile, outputFileName, rowNames=NULL, colNames=NULL, counts=TRUE, returnOverlapTable=TRUE, pathToIntersectBedToolsLoop="/Users/johnurban/searchPaths/shell_scripts_Mac/intersectBedLoop.sh"){
  ## Began 9/4/13
  ## Depends on intersectBedLoop.sh, which depends on BEDtools
  ## Where pathToIntersectBedLoop is the absolute path to this shell script
  ##    default on MAC: "/Users/johnurban/searchPaths/shell_scripts_Mac/intersectBedLoop.sh"
  ## Where 'filePathsTxtFile' is a text file with the paths (in order) to each BED file to compare to all others
  ## This is a wrapper just so everything can be done in R if one wants -- or this can be part of a larger pipeline
  ##    Otherwise - no need to do it in R per se and one can run "intersectBedLoop.sh" from the commandline
  ## This, by default, will also read the table into the R workspace (the default return of fxn) -- can turn this off and read in manually if want
  
  ### 9/4/13 This is broken still b/c 'intersectBed' path is not specified inside the shell script...
  ### My advice is to discontinue efforts here and just do it at the bash commandline
  if(counts){
    style="counts"
  } else {style="percents"}
  intersectBedLoop = paste0(pathToIntersectBedToolsLoop, " ", filePathsTxtFile, " ", style, " > ", outputFileName)
  system(intersectBedLoop)
  
  ###TODO
  if(returnOverlapTable){
    if(is.null(rowNames) | is.null(colNames)){return("Must provide rowNames and/or colNames")}
    return("Passed")
  } 
}