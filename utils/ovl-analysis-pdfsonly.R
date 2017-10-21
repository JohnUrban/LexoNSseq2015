args <- commandArgs(trailingOnly=TRUE)

narrowpeak <- args[1]
source(narrowpeak)

Pairwise <- read.table(args[2], header=TRUE)

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


## Need to read in Tables for this...
### pvalue matrix
##pdf("logpValmatrix.pdf")
##smallestRnum <- 4.940656e-324
##matrixPlot(scoreMatrix=-log10(Tables$pValueTable+smallestRnum), title="")
##garbage <- dev.off()

