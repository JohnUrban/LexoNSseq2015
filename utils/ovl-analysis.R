args <- commandArgs(trailingOnly=TRUE)

narrowpeak <- args[1]
source(narrowpeak)

###genome <- args[2]
genomeFileVar <- readGenomeFile(args[2])
gapFileVar <- readGapFile(args[3])
numConnComp <- args[4]
ovlptable <- args[5]
namesfile <- args[6]
pathsfile <- args[7]
minovlp <- args[8]

mapGenSize <- mappableGenomeSize(genomeFileVar, gapFileVar)

ObsOverlaps <- readOverlapTable(overlapTableTxtFile=ovlptable, rowNames=readNamesFile(namesTxtFile=namesfile), colNames=readNamesFile(namesTxtFile=namesfile))

print("Tables...")
Tables <- binomialOverLapPval(obsNumPeaksOverlapTable=ObsOverlaps, bedFileListFile=pathsfile, genomeSize=mapGenSize, numConnectedComponents=numConnComp, minOverlap=minovlp)
suppressWarnings( writeAllTables(Tables, "Tables.txt") )

print("Pairwise...")
Pairwise <- pairwiseBinomialInfo(Tables)
writePairwiseBinomialInfo(pairwiseTable=Pairwise, outputFilename="PairwiseBinomialInfo.txt")


##
pdf("ExpObsChart.pdf")
par(mfrow=c(2,2), mar=c(8,3,2,2))

## fig1: exp vs. obs
plotExpVsObsOverlap(Pairwise, noSelf=TRUE, ylim=c(0,1.3), legendX=c(1,1.3), yaxt="n")
axis(2, at=seq(0,1,0.1), labels=seq(0,1,0.1), las=1)

## fig2: obs:exp ratio
plotObsExpOverlapRatio(pairwiseTable=Pairwise)

## fig3: pval
plotlog10BinomialpVal(pairwiseTable=Pairwise)

## fig 4: maxPoss, percentMaxPoss obtained
plotApproxMaxWithPercentMax(pairwiseTable=Pairwise, col=c("dark red","dark blue"))
garbage <- dev.off()


### pvalue matrix
pdf("logpValmatrix.pdf")
smallestRnum <- 4.940656e-324
matrixPlot(scoreMatrix=-log10(Tables$pValueTable+smallestRnum), title="")
garbage <- dev.off()

