# FUNCTIONS LIST:
# For all peaks:
# - generateNsomePlot() - the main function used to generate our nsome figures for all peaks
#                             -- defaults: plotnsome=T, plotg4=F, g4Spacing=T, Nspacing=F, g4SpacingLinesType=5, window=1000, sample="nsg", cell.line="K562", nsomeScore="mean", nsomeLoessSpan=0.1, g4color=rgb(1,0,1,0.2), segCol="red", nsomeCol=NA, nsomeSummitParams=c(80,30,80,1.05), g4SummitParams=c(80,70,150,1.5)
# - exampleNsomePlotterUsage() -- quick, albeit slightly outdated (b/c I developed it more) tutorial on using generateNsomePlot()
#                             -- no parameters, just execute
# - exampleG4SpacingLineTypes() -- a wrapper over generateNsomePlot() that somewhat simplifies execution. This was actually used to generate most figures.
#                             -- if i=NA, it will go through all spacing patterns with other given params (e.g. sample and cell.line)
#                             -- defaults: i=NA, sample='nsl', cell.line="mean", g4color=rgb(1,0,0,0.5), segCol="red", nsomeCol=NA
#                             -- hardCoded defaults into generateNsomePlot(): plotnsome=T, plotg4=F, g4Spacing=T, Nspacing=F, g4SpacingLinesType=i, window=1000, sample=sample, cell.line=cell.line, nsomeScore="mean", nsomeLoessSpan=0.05, g4color=g4color, segCol=segCol, nsomeCol=nsomeCol, nsomeSummitParams=c(80,30,80,1.05), g4SummitParams=c(80,70,150,1.5)
#                             --> also any options not set default to generateNSomePlots()s defaults
# - generateNsomePlotAlt() - a slightly outdated and deviated version of generateNsomePlot() as the name implies
#                           - this will let you plot Nsome signal around summits and shuffled summits instead of just the FE from them
#                           - FE is optional and here it is FE of pre-smoothed signal, whereas in the original it smooths FE signal

#
readNsomeSummaryFile <-function(filename){
  data <- read.table(filename)
  colnames(data) <- c("pos", "sum", "mean", "numDataPoints")
  return(data)
}

generateNsomePlot <- function(plotnsome=T, plotg4=F, g4Spacing=T, Nspacing=F, g4SpacingLinesType=5, window=1000, sample="nsg", cell.line="K562", nsomeScore="mean", nsomeLoessSpan=0.1, g4color=rgb(1,0,1,0.2), segCol="red", nsomeCol=NA, nsomeSummitParams=c(80,30,80,1.05), g4SummitParams=c(80,70,150,1.5)){
  ## cell.line can be "K562" or "Gm12878", "both", "mean"
  ## sample can be "nsg", "nsl", "lg0"
  ## g4SpacingLinesType can be 1,2,3,4,6, 7 or 8, which tries to represent G4 summit height by :
  ##      lineWidthOnly, colorOnly, height only, lineWidthAndColor, lineWidthAndHeight, lineWidthAndHeightAndColor, polygon of entire G4 signal transformed same way as height lines, polygon+height lines
  ##    color of lines is defaulted as red unless color scale used, then blue to red.
  par(mar=c(5,5,5,5))
  if (sample == "nsg"){
    textcex <- 1
    if(is.na(nsomeCol)){nsomeCol <- "dark blue"}
    samplebase <- 10
    g4subtract <- 0.66
    g4linesdecrement <- 0.01
  } else if (sample == "nsl"){
    textcex <- 0.7
    if(is.na(nsomeCol)){nsomeCol <- "dark green"}
    samplebase <- 10
    g4subtract <- 1.35
    g4linesdecrement <- 0.01
  } else if (sample == "lg0"){
    textcex <- 1
    if(is.na(nsomeCol)){nsomeCol <- "dark cyan"}
    samplebase <- 500
    g4subtract <- 0.29
    g4linesdecrement <- 0.002
  }
  if (cell.line %in% c("K562", "Gm12878")){
    signal1 <- nsomeData[[sample]][[cell.line]][["Peaks"]][[nsomeScore]]/nsomeData[[sample]][[cell.line]][["shufPeaks"]][[nsomeScore]]
    pos1 <- nsomeData[[sample]][[cell.line]][["Peaks"]]$pos
    signal2 <- 0 ## zero so it does not affect ylim decision later
  }
  else if(cell.line %in% c("both", "mean")){
    signal1 <- nsomeData[[sample]][["K562"]][["Peaks"]][[nsomeScore]]/nsomeData[[sample]][["K562"]][["shufPeaks"]][[nsomeScore]]
    signal2 <- nsomeData[[sample]][["Gm12878"]][["Peaks"]][[nsomeScore]]/nsomeData[[sample]][["Gm12878"]][["shufPeaks"]][[nsomeScore]]
    pos1 <- nsomeData[[sample]][["K562"]][["Peaks"]]$pos
    pos2 <- nsomeData[[sample]][["Gm12878"]][["Peaks"]]$pos
    meanSignal <- (signal1+signal2)/2
    if (sum(pos1 == pos2) != length(pos1)){ return("Positions Error...see code")}
    if (cell.line == "mean") { 
      signal1 <- (signal1+signal2)/2 
      signal2 <- 0 ## zeroing this out so it does not affect ylim decision later
    }
  } 
  if(plotg4){
    lg4 <- loess(g4Data[[sample]][["G4"]]$count/g4Data[[sample]][["shufG4"]]$count ~ g4Data[[sample]][["G4"]]$mid, span=0.1)
  }
  ylim <- c(min(c(signal1[signal1 > 0], signal2[signal2 > 0])), max(c(signal1, signal2)))
  basicLwd <- 2
  
  ##PLOTTING
  loessSamp1 <- loess(signal1 ~ pos1, span=nsomeLoessSpan)
  datapoints <- loessSamp1$x >= -window & loessSamp1$x <= window
  if (plotnsome & plotg4){
    return("In Devo")
  }
  else if (plotnsome){
    plot(pos1, signal1, xaxt="n",type="n", ylim=ylim, ylab="Average Nucleosome Signal", xlab="Relative Position From Summit (bp)", las=2)
    lines(loessSamp1$x[datapoints], loessSamp1$fitted[datapoints], col=nsomeCol, lwd=basicLwd)
    if (cell.line == "both"){
      loessSamp2 <- loess(signal2 ~ pos2, span=nsomeLoessSpan)
      lines(loessSamp2$x[datapoints], loessSamp2$fitted[datapoints], col=nsomeCol, lwd=basicLwd)
    }
  }
  else if (plotg4){plot(lg4$x, lg4$fitted, ylim=c(0, max(lg4$fitted)), col=nsomeCol, ylab="G4 Enrichment", xlab="Relative Position (bp)", las=2, type="l", xaxt="n")}
  axis(side=1,at=seq(-1000,1000,200), labels=seq(-1000,1000,200), las=1, cex.axis=0.8)
  abline(v=0, lty=2)

  ## spacing
  if (Nspacing){
    if (cell.line == "both"){loessSamp1 <- loess(meanSignal ~ pos1, span=nsomeLoessSpan)}
    summits <- loessSamp1$x[allsummits5(vector=loessSamp1$fitted, controlVector=NA, vectorAlreadyIsFE=T, maxWindowSearch=nsomeSummitParams[1], minWindowSearch=nsomeSummitParams[2], totalWinLenCutOff=nsomeSummitParams[3], FEcut=nsomeSummitParams[4])] ## allsummits5 c(80,30,80,1.05)
    InterLens <- summits[2:length(summits)] - summits[1:(length(summits)-1)]
    abline(v=summits, col="black", lwd=3, lty=1)
    yArrow <- quantile(x=loessSamp1$fitted, probs=0.98)
    yText <- quantile(x=loessSamp1$fitted, probs=1)
    for(i in 1:(length(summits)-1)){arrows(x0=summits[i],y0=yArrow,x1=summits[i+1],y1=yArrow, length=0.1, lwd=2)}
    for(i in 1:(length(summits)-1)){arrows(x0=summits[i+1],y0=yArrow,x1=summits[i],y1=yArrow, length=0.1, lwd=2)}
    for(i in 1:(length(summits)-1)){text(x=mean(c(summits[i],summits[i+1])),y=yText,labels=InterLens[i], cex=textcex)}
  }
  # G4s
  if (g4Spacing){
    lg4Samp1 <- loess(g4Data[[sample]][["G4"]]$counts ~ g4Data[[sample]][["G4"]]$mids, span=0.1)
    lg4FE1 <- loess(g4Data[[sample]][["G4"]]$counts/g4Data[[sample]][["shufG4"]]$counts ~ g4Data[[sample]][["G4"]]$mids, span=0.1)
    ## summits obtained from counts data, FE plotted
    summits <- lg4Samp1$x[allsummits5(vector=lg4Samp1$fitted, controlVector=g4Data[[sample]][["shufG4"]]$counts, vectorAlreadyIsFE=F, maxWindowSearch=g4SummitParams[1], minWindowSearch=g4SummitParams[2], totalWinLenCutOff=g4SummitParams[3], FEcut=g4SummitParams[4])] ## allsummits5, c(80,70,150,1.5)
    ########### above could maybe also index loessSamp1$x...
    InterLens <- summits[2:length(summits)] - summits[1:(length(summits)-1)]
#     lwdsDenominator <- sum(lg4FE1$fitted[which(lg4FE1$x %in% summits)])
    lwdsDenominator <- max(lg4FE1$fitted[which(lg4FE1$x %in% summits)])
    lwds <- lg4FE1$fitted[which(lg4FE1$x %in% summits)]#/lwdsDenominator
    g4ColorsFxn <- colorRamp(colors=c("blue","grey", "red"))
    g4Colors <- g4ColorsFxn((lwds-min(lwds))/max(lwds-min(lwds)))/255
    ranks <- (1:length(lwds))[order(lwds)]
    rankColors <- g4ColorsFxn(ranks/max(ranks))/255
#     points(500:600, rep(1.5, 101), col=rgb(g4ColorsFxn(seq(0, 1, 0.01))/255), pch=19, cex=0.1)
#     lineWidthOnly, colorOnly, height only, lineWidthAndColor, lineWidthAndHeight, lineWidthAndHeightAndColor, polygon of entire G4 signal transformed same way as height lines
    lineWidths <- lwds*0.75
    if (g4SpacingLinesType == 1){abline(v=summits, col="red", lwd=lineWidths, lty=1)}
    else if (g4SpacingLinesType == 2){abline(v=summits, col=rgb(g4Colors), lwd=2, lty=1)}
    else if (g4SpacingLinesType == 4){abline(v=summits, col=rgb(g4Colors), lwd=lineWidths, lty=1)}
    else if (g4SpacingLinesType %in% c(3,5,6)){
      for (i in 1:length(summits)){
          if (g4SpacingLinesType == 6){
            segments(x0=summits[i], y0=0, x1=summits[i], y1=log(lwds[i],base=samplebase)+ylim[2]-ylim[1]*g4subtract, col=rgb(g4Colors[i,1], g4Colors[i,2], g4Colors[i,3]), lwd=lineWidths[i])
          }
          else if (g4SpacingLinesType == 5){
            segments(x0=summits[i], y0=0, x1=summits[i], y1=log(lwds[i],base=samplebase)+ylim[2]-ylim[1]*g4subtract, col="red", lwd=lineWidths[i])
          }
          else {
            segments(x0=summits[i], y0=0, x1=summits[i], y1=log(lwds[i],base=samplebase)+ylim[2]-ylim[1]*g4subtract, col="red", lwd=2)
          }
      }
    }
    else if (g4SpacingLinesType %in% c(7,8,9,10)){
      ## polygons
#       allheights <- c(lg4FE1$fitted/sum(lg4FE1$fitted))
      allheights <- lg4FE1$fitted#/max(lg4FE1$fitted)
      pos <- c(lg4FE1$x[1]-1, lg4FE1$x, lg4FE1$x[length(lg4FE1$x)]+1)
      #heights <- c(0, log10((10^ylim[2])*allheights)+(ylim[2]-max(log10((10^ylim[2])*allheights))), 0)
      #heights2 <- c(0, log10((10^ylim[2])*lwds[i])+(ylim[2]-max(log10((10^ylim[2])*lwds)))
#       if (sample %in% c("lg0")){heights <- c(0, log10((10^ylim[2])*allheights)+(ylim[2]-max(log10((10^ylim[2])*allheights))), 0)}
#       else if (sample %in% c("nsg")){heights <- c(0, log2((10^ylim[2])*allheights)+(ylim[2]-max(log2((10^ylim[2])*allheights)))-ylim[1]*0.1, 0)}
#       else if (sample %in% c("nsl")){heights <- c(0, log10((10^ylim[2])*allheights)+(ylim[2]-max(log10((10^ylim[2])*allheights)))-ylim[1]*0.5, 0)}
#       if (sample %in% c("lg0")){heights <- c(0, log(allheights/max(allheights), base=samplebase)+ylim[2]-ylim[1]/20, 0)}
#       else if (sample %in% c("nsg")){heights <- c(0, log(allheights/max(allheights), base=samplebase)+ylim[2]-ylim[1]*0.1, 0)}
#       else if (sample %in% c("nsl")){heights <- c(0, log(allheights/max(allheights), base=samplebase)+ylim[2]-ylim[1]*0.5, 0)}
      if (sample %in% c("lg0")){heights <- c(0, log(allheights, base=samplebase)+ylim[2]-ylim[1]*0.29, 0)}
      else if (sample %in% c("nsg")){heights <- c(0, log(allheights, base=samplebase)+ylim[2]-ylim[1]*0.66, 0)}
      else if (sample %in% c("nsl")){heights <- c(0, log(allheights, base=samplebase)+ylim[2]-ylim[1]*1.35, 0)}
#       heights <- c(0, log10(allheights), 0)*ylim[2]/max(log10(lwds))
#       polygon(x=pos, y=heights, col=g4color, border=NA)
      if (g4SpacingLinesType %in% c(7,8)){polygon(x=pos, y=heights, col=g4color, border="black", lwd=0.5)}
      else if (g4SpacingLinesType %in% c(9,10)){
        lines(x=pos, y=heights, col=g4color)
        for (i in seq(0, max(heights), g4linesdecrement)){lines(x=pos, y=heights-i, col=g4color)}
      }
#       polygon(x=pos, y=heights, col=rgb(0.3,0.1,180/255,0.3), border=NA)
#       polygon(x=pos, y=heights, col=rgb(0.6627451,0.6627451,0.6627451,0.3), border=NA)
      ## redraw nsome lines over polygon
      if (g4SpacingLinesType %in% c(7,8,9)){lines(loessSamp1$x[datapoints], loessSamp1$fitted[datapoints], col=nsomeCol, lwd=basicLwd)}
      if (cell.line == "both"){lines(loessSamp2$x[datapoints], loessSamp2$fitted[datapoints], col=nsomeCol, lwd=basicLwd)}
      ## if 8, add height lines
      if (g4SpacingLinesType %in% c(8,10)){
        if(g4SpacingLinesType == 8){segwd <- 1}
        else if(g4SpacingLinesType == 10){segwd <- 2}
        for (i in 1:length(summits)){
          if (sample %in% c("lg0")) {
            segments(x0=summits[i], y0=0, x1=summits[i], y1=log(lwds[i], base=samplebase)+ylim[2]-ylim[1]*0.29, col=segCol, lwd=segwd, lty=1)
            ## taking max of lwds is same as max of allheights b/c lwds are the summits inside allheights
#             segments(x0=summits[i], y0=0, x1=summits[i], y1=log10((10^ylim[2])*lwds[i])+(ylim[2]-max(log10((10^ylim[2])*allheights))), col="red", lwd=lineWidths)
          }
          else if (sample %in% c("nsg")) {
            segments(x0=summits[i], y0=0, x1=summits[i], y1=log(lwds[i], base=samplebase)+ylim[2]-ylim[1]*0.66, col=segCol, lwd=segwd,lty=1)
#             segments(x0=summits[i], y0=0, x1=summits[i], y1=log10((10^ylim[2])*lwds[i])+(ylim[2]-max(log10((10^ylim[2])*allheights)))-ylim[1]*0.1, col="red", lwd=lineWidths)
          }
          else if (sample %in% c("nsl")) {
            segments(x0=summits[i], y0=0, x1=summits[i], y1=log(lwds[i], base=samplebase)+ylim[2]-ylim[1]*1.35, col=segCol, lwd=segwd, lty=1)
#             segments(x0=summits[i], y0=0, x1=summits[i], y1=log10((10^ylim[2])*lwds[i])+(ylim[2]-max(log10((10^ylim[2])*allheights)))-ylim[1]*0.5, col="red", lwd=lineWidths)
          }
          
        }
      }      
    }
    if (g4SpacingLinesType %in% c(7,8,9)){lines(loessSamp1$x[datapoints], loessSamp1$fitted[datapoints], col=nsomeCol, lwd=basicLwd)}
    
    ## redraw nsome lines over styles other than polygone looked like shit
    ## Add side 4 axis
    if (g4SpacingLinesType %in% c(3,5,6,7,8)){
      if (sample == "lg0"){
        f <- function(x){samplebase^(x-ylim[2]+ylim[1]*0.29)}
        g <- function(x){log(x, base=samplebase)+ylim[2]-ylim[1]*0.29}
        at <- g(c(0.5,0.7,1,1.5,2,2.5,3:6, 8, 10))
        labels <- round(f(at), digits=1)
        at <- g(labels)
#         abline(h=g(max(lwds))); print(max(lwds))
      } 
      else if (sample == "nsg"){
        f <- function(x){samplebase^(x-ylim[2]+ylim[1]*0.66)}
        g <- function(x){log(x, base=samplebase)+ylim[2]-ylim[1]*0.66}
        at <- g(c(0.5,0.7,1,1.5,2,2.5,3:6, 8, 10))
        labels <- round(f(at), digits=1)
        at <- g(labels)
#         abline(h=g(max(lwds))); print(max(lwds))
      } 
      else if (sample == "nsl"){
        f <- function(x){samplebase^(x-ylim[2]+ylim[1]*1.35)}
        g <- function(x){log(x, base=samplebase)+ylim[2]-ylim[1]*1.35}
        at <- g(c(0.5,0.7,1,1.5,2,2.5,3:6, 8, 10))
        labels <- round(f(at), digits=1)
        at <- g(labels)
#         abline(h=g(max(lwds))); print(max(lwds))
      }
      mtext("G4 Fold Enrichment", side=4, line=3)
      axis(side=4, at=at, labels=labels, las=2)
    }
    ## spacing info
    if (g4SpacingLinesType != 7){
      yvals <- seq(ylim[1],ylim[2], 0.01)
      yArrow <- quantile(x=yvals,probs=0)
      yText <- quantile(x=yvals,probs=0.03)
      textcex<-1
      for(i in 1:(length(summits)-1)){arrows(x0=summits[i],y0=yArrow,x1=summits[i+1],y1=yArrow, length=0.1, lwd=2)}
      for(i in 1:(length(summits)-1)){arrows(x0=summits[i+1],y0=yArrow,x1=summits[i],y1=yArrow, length=0.1, lwd=2)}
      for(i in 1:(length(summits)-1)){text(x=mean(c(summits[i],summits[i+1])),y=yText,labels=InterLens[i], cex=textcex)}
    }
  }
}




exampleNsomePlotterUsage <- function(){
  print("Examples using lg0")
  print("note: sample and cell.line values are actually strings though they look like variables here")
  par(ask=T)
  generateNsomePlot(g4Spacing=T, Nspacing=F, window=1000, sample="lg0", cell.line="K562", nsomeLoessSpan=0.05, nsomeSummitParams=c(80,30,80,1.05), g4SummitParams=c(80,70,150,1.5))
  print("K562 with G4 distances")
  print("generateNsomePlot(g4Spacing=T, Nspacing=F, window=1000, sample=lg0, cell.line=K562, nsomeLoessSpan=0.05, nsomeSummitParams=c(80,30,80,1.05), g4SummitParams=c(80,70,150,1.5))")
  
  par(ask=T)
  generateNsomePlot(g4Spacing=T, Nspacing=F, window=1000, sample="lg0", cell.line="Gm12878", nsomeLoessSpan=0.05, nsomeSummitParams=c(80,30,80,1.05), g4SummitParams=c(80,70,150,1.5))
  print("Gm12878 with G4 distances")
  print("generateNsomePlot(g4Spacing=T, Nspacing=F, window=1000, sample=\"lg0\", cell.line=\"Gm12878\", nsomeLoessSpan=0.05, nsomeSummitParams=c(80,30,80,1.05), g4SummitParams=c(80,70,150,1.5))")

  par(ask=T)
  generateNsomePlot(g4Spacing=T, Nspacing=F, window=1000, sample="lg0", cell.line="both", nsomeLoessSpan=0.05, nsomeSummitParams=c(80,30,80,1.05), g4SummitParams=c(80,70,150,1.5))
  print("both cell lines")
  print("generateNsomePlot(g4Spacing=T, Nspacing=F, window=1000, sample=lg0, cell.line=both, nsomeLoessSpan=0.05, nsomeSummitParams=c(80,30,80,1.05), g4SummitParams=c(80,70,150,1.5))")
  
  par(ask=T)
  generateNsomePlot(g4Spacing=T, Nspacing=F, window=1000, sample="lg0", cell.line="mean", nsomeLoessSpan=0.05, nsomeSummitParams=c(80,30,80,1.05), g4SummitParams=c(80,70,150,1.5))
  print("mean from both cell lines")
  print("generateNsomePlot(g4Spacing=T, Nspacing=F, window=1000, sample=lg0, cell.line=mean, nsomeLoessSpan=0.05, nsomeSummitParams=c(80,30,80,1.05), g4SummitParams=c(80,70,150,1.5))")
  
  par(ask=T)
  generateNsomePlot(g4Spacing=F, Nspacing=T, window=1000, sample="lg0", cell.line="both", nsomeLoessSpan=0.05, nsomeSummitParams=c(80,30,80,1.05), g4SummitParams=c(80,70,150,1.5))
  print("Measuring nsome spacing -- if 'both' given for cell.line, it measures spacing of mean")
  print("...See spacing params at the end of command -- these are fed into allsummits5")
  print("both")
  print("generateNsomePlot(g4Spacing=F, Nspacing=T, window=1000, sample=lg0, cell.line=Gm12878, nsomeLoessSpan=0.05, nsomeSummitParams=c(80,30,80,1.05), g4SummitParams=c(80,70,150,1.5))")
  
  par(ask=T)
  generateNsomePlot(g4Spacing=F, Nspacing=T, window=1000, sample="lg0", cell.line="mean", nsomeLoessSpan=0.05, nsomeSummitParams=c(80,30,80,1.05), g4SummitParams=c(80,70,150,1.5))
  print("mean")
  print("generateNsomePlot(g4Spacing=F, Nspacing=T, window=1000, sample=lg0, cell.line=mean, nsomeLoessSpan=0.05, nsomeSummitParams=c(80,30,80,1.05), g4SummitParams=c(80,70,150,1.5))")
  
  par(ask=T)
  generateNsomePlot(g4Spacing=T, Nspacing=T, window=1000, sample="lg0", cell.line="mean", nsomeLoessSpan=0.05, nsomeSummitParams=c(80,30,80,1.05), g4SummitParams=c(80,70,150,1.5))
  print("Both spacings can be on plot by setting: g4Spacing=T, Nspacing=T")
  
  par(ask=T)
  generateNsomePlot(g4Spacing=T, Nspacing=F, window=1000, sample="lg0", cell.line="mean", nsomeLoessSpan=0.05, nsomeSummitParams=c(80,30,80,1.05), g4SummitParams=c(80,70,150,1.5))
  print("Messing around with windows of data with window=x")
  print("Window=1000 -- what we have been seeing")
  
  par(ask=T)
  generateNsomePlot(g4Spacing=T, Nspacing=F, window=500, sample="lg0", cell.line="mean", nsomeLoessSpan=0.05, nsomeSummitParams=c(80,30,80,1.05), g4SummitParams=c(80,70,150,1.5))
  print("window=500 -- can pick any number <=1000")  
  
  par(ask=F)
}
  
exampleG4SpacingLineTypes <- function(i=NA, sample='nsl', cell.line="mean", g4color=rgb(1,0,0,0.5), segCol="red", nsomeCol=NA){
  print(sample)
  g4SpacingLinesTypeMsg <- c("1=lineWidthOnly", "2=colorOnly", "3=height only", "4=lineWidthAndColor", "5=lineWidthAndHeight", "6=lineWidthAndHeightAndColor", "7=polygon of entire G4 signal transformed same way as height lines -- spacing info not included (use option 8 for that to be added)", "8=polygon of entire G4 signal transformed same way as height lines as well as the heightOnly lines and spacing information")
  if (is.na(i)){
    par(ask=T)
    for(i in 1:8){
      generateNsomePlot(plotnsome=T, plotg4=F, g4Spacing=T, Nspacing=F, g4SpacingLinesType=i, window=1000, sample=sample, cell.line=cell.line, nsomeScore="mean", nsomeLoessSpan=0.05, g4color=g4color, segCol=segCol, nsomeCol=nsomeCol, nsomeSummitParams=c(80,30,80,1.05), g4SummitParams=c(80,70,150,1.5))
      print(g4SpacingLinesTypeMsg[i])
    }
    par(ask=F)
  } else {
    generateNsomePlot(plotnsome=T, plotg4=F, g4Spacing=T, Nspacing=F, g4SpacingLinesType=i, window=1000, sample=sample, cell.line=cell.line, nsomeScore="mean", nsomeLoessSpan=0.05, g4color=g4color, segCol=segCol, nsomeCol=nsomeCol, nsomeSummitParams=c(80,30,80,1.05), g4SummitParams=c(80,70,150,1.5))
    print(g4SpacingLinesTypeMsg[i])
  }
}

## for plotting the treatment and random together (not FE) by default
## if FE=TRUE, it plots FE of pre-smoothed treat and pre-smoothed shuf 
##  the other function does smoothing of FEd treat/shuf
generateNsomePlotAlt <- function(g4Spacing=T, Nspacing=F, FE=F,  window=500, sample="nsg", cell.line="K562", nsomeScore="mean", nsomeLoessSpan=0.1, nsomeSummitParams=c(80,30,80,1.05), g4SummitParams=c(80,70,150,1.5)){
  ## cell.line can be "K562" or "Gm12878", "both", "mean"
  ## sample can be "nsg", "nsl", "lg0"
  ## nsomeScore can be "mean" or "sum"
  if (sample == "nsg"){
    textcex <- 1
    nsomeCol <- "dark blue"
  } else if (sample == "nsl"){
    textcex <- 0.7
    nsomeCol <- "dark green"
  } else if (sample == "lg0"){
    textcex <- 1
    nsomeCol <- "dark cyan"
  }
  if (cell.line %in% c("K562", "Gm12878")){
    signal1 <- nsomeData[[sample]][[cell.line]][["Peaks"]][[nsomeScore]]
    shuf1 <- nsomeData[[sample]][[cell.line]][["shufPeaks"]][[nsomeScore]]
    pos1 <- nsomeData[[sample]][[cell.line]][["Peaks"]]$pos
    signal2 <- 0 ## zero so it does not affect ylim decision later
  }
  else if(cell.line %in% c("both", "mean")){
    signal1 <- nsomeData[[sample]][["K562"]][["Peaks"]][[nsomeScore]]
    shuf1 <- nsomeData[[sample]][["K562"]][["shufPeaks"]][[nsomeScore]]
    signal2 <- nsomeData[[sample]][["Gm12878"]][["Peaks"]][[nsomeScore]]
    shuf2 <- nsomeData[[sample]][["Gm12878"]][["shufPeaks"]][[nsomeScore]]
    pos1 <- nsomeData[[sample]][["K562"]][["Peaks"]]$pos
    pos2 <- nsomeData[[sample]][["Gm12878"]][["Peaks"]]$pos
    meanSignal <- (signal1+signal2)/2
    meanShuf <- (shuf1+shuf2)/2
    if (sum(pos1 == pos2) != length(pos1)){ return("Positions Error...see code")}
    if (cell.line == "mean") { 
      signal1 <- (signal1+signal2)/2 
      shuf1 <- (shuf1+shuf2)/2
      signal2 <- 0 ## zeroing this out so it does not affect ylim decision later
    }
  } 

  ylim <- c(0, max(c(signal1, signal2)))

  
  ##PLOTTING
  loessSamp1 <- loess(signal1 ~ pos1, span=nsomeLoessSpan)
  datapoints <- loessSamp1$x >= -window & loessSamp1$x <= window
  loessShuf1 <- loess(shuf1 ~ pos1, span=nsomeLoessSpan)
  plot(pos1, signal1, xaxt="n",type="n", ylim=ylim, ylab="Average Nucleosome Signal", xlab="Relative Position From Summit (bp)", las=2)
  if(FE){
    fe1 <- loessSamp1$fitted[datapoints]/loessShuf1$fitted[datapoints]
    lines(loessSamp1$x[datapoints], fe1, col=nsomeCol, lwd=2)
  } else {
    lines(loessSamp1$x[datapoints], loessSamp1$fitted[datapoints], col=nsomeCol, lwd=2)
    lines(loessShuf1$x[datapoints], loessShuf1$fitted[datapoints], col="dark grey", lwd=2) 
  }
  if (cell.line == "both"){
    loessSamp2 <- loess(signal2 ~ pos2, span=nsomeLoessSpan)
    loessShuf2 <- loess(shuf2 ~ pos2, span=nsomeLoessSpan)
    if(FE){
      fe2 <- loessSamp2$fitted[datapoints]/loessShuf2$fitted[datapoints]
      lines(loessSamp2$x[datapoints], fe2, col=nsomeCol, lwd=2)
    } else {
      lines(loessSamp2$x[datapoints], loessSamp2$fitted[datapoints], col=nsomeCol, lwd=2)
      lines(loessShuf2$x[datapoints], loessShuf2$fitted[datapoints], col="dark grey", lwd=2)
    }
  }
  axis(side=1,at=seq(-1000,1000,200), labels=seq(-1000,1000,200), las=1, cex.axis=0.8)
  abline(v=0, lty=2)
  
  ## spacing
  if (Nspacing){
    if (cell.line == "both"){loessSamp1 <- loess(meanSignal ~ pos1, span=nsomeLoessSpan)}
    summits <- loessSamp1$x[allsummits5(vector=loessSamp1$fitted, controlVector=NA, vectorAlreadyIsFE=T, maxWindowSearch=nsomeSummitParams[1], minWindowSearch=nsomeSummitParams[2], totalWinLenCutOff=nsomeSummitParams[3], FEcut=nsomeSummitParams[4])] ## allsummits5 c(80,30,80,1.05)
    InterLens <- summits[2:length(summits)] - summits[1:(length(summits)-1)]
    abline(v=summits, col="black", lwd=3, lty=1)
    yArrow <- quantile(x=loessSamp1$fitted, probs=0.98)
    yText <- quantile(x=loessSamp1$fitted, probs=1)
    for(i in 1:(length(summits)-1)){arrows(x0=summits[i],y0=yArrow,x1=summits[i+1],y1=yArrow, length=0.1, lwd=2)}
    for(i in 1:(length(summits)-1)){arrows(x0=summits[i+1],y0=yArrow,x1=summits[i],y1=yArrow, length=0.1, lwd=2)}
    for(i in 1:(length(summits)-1)){text(x=mean(c(summits[i],summits[i+1])),y=yText,labels=InterLens[i], cex=textcex)}
  }
  # G4s
  if (g4Spacing){
    lg4Samp1 <- loess(g4Data[[sample]][["G4"]]$counts ~ g4Data[[sample]][["G4"]]$mids, span=0.1)
    summits <- lg4Samp1$x[allsummits5(vector=lg4Samp1$fitted, controlVector=g4Data[[sample]][["shufG4"]]$counts, vectorAlreadyIsFE=F, maxWindowSearch=g4SummitParams[1], minWindowSearch=g4SummitParams[2], totalWinLenCutOff=g4SummitParams[3], FEcut=g4SummitParams[4])] ## allsummits5, c(80,70,150,1.5)
    ########### above could maybe also index loessSamp1$x...
    InterLens <- summits[2:length(summits)] - summits[1:(length(summits)-1)]
    lwds <- lg4Samp1$fitted[which(lg4Samp1$x %in% summits)]/sum(lg4Samp1$fitted[which(lg4Samp1$x %in% summits)])
    g4ColorsFxn <- colorRamp(colors=c("blue","grey", "red"))
    g4Colors <- g4ColorsFxn((lwds-min(lwds))/max(lwds-min(lwds)))/255
    #     points(500:600, rep(1.5, 101), col=rgb(g4ColorsFxn(seq(0, 1, 0.01))/255), pch=19, cex=0.1)
    #     abline(v=summits, col="blue", lwd=lwds*15, lty=1)
    abline(v=summits, col=rgb(g4Colors), lwd=lwds*20, lty=1)
    #     abline(v=summits, col=rgb(g4Colors), lwd=2, lty=1)
    yArrow <- 1.0
    yText <- 1.05
    for(i in 1:(length(summits)-1)){arrows(x0=summits[i],y0=yArrow,x1=summits[i+1],y1=yArrow, length=0.1, lwd=2)}
    for(i in 1:(length(summits)-1)){arrows(x0=summits[i+1],y0=yArrow,x1=summits[i],y1=yArrow, length=0.1, lwd=2)}
    for(i in 1:(length(summits)-1)){text(x=mean(c(summits[i],summits[i+1])),y=yText,labels=InterLens[i], cex=textcex)}
  }
}

