


# FUNCTIONS LIST:
# see each sourced file for FXN descriptions
#
# For new analysis focusing on partition peaks

prelimWrapper <- function(sample="nsl", withG4s=TRUE, cell.line="mean", addNsomeSpacing=TRUE, yArrow=1.2,yText=1.25, textcex=0.75){
  
}
gen3Stack.K562 <- function(){
  plot(ylab="Nucleosome Signal", xlab="Relative Position from Summit (bp)", yaxt="n",las=1,k_nsl_wG4$pos, loess(k_nsl_wG4$mean/k_nslshuf$mean ~ k_nsl_wG4$pos, span=0.075)$fitted, type="n", ylim=c(0,1.2))
  abline(h=c(0.4,0.8))
  axis(side=2, at=c(c(0.8, 0.9, 1)-0.7,c(0.8, 0.9, 1)-0.3,c(0.9, 1, 1.1,1.2)), labels=c(c(0.8, 0.9, 1),c(0.8, 0.9, 1),c(0.9, 1, 1.1,1.2)), las=1)
  # abline(h=c(1,1-0.3,1-0.7), lty=3)
  lines(k_nsl_wG4$pos, loess(k_nsl_wG4$mean/k_nslshuf$mean ~ k_nsl_wG4$pos, span=0.075)$fitted-0.7, col="dark green", lwd=3)
  lines(k_nsg_wG4$pos, loess(k_nsg_wG4$mean/k_nsgshuf$mean ~ k_nsg_wG4$pos, span=0.075)$fitted-0.3, col="dark blue", lwd=2)
  lines(k_lg0_wG4$pos, loess(k_lg0_wG4$mean/k_nsgshuf$mean ~ k_lg0_wG4$pos, span=0.075)$fitted, col="dark cyan", lwd=2)
  
}

gen3Stack.Gm12878 <- function(){
  plot(ylab="Nucleosome Signal", xlab="Relative Position from Summit (bp)", yaxt="n",las=1,k_nsl_wG4$pos, loess(k_nsl_wG4$mean/k_nslshuf$mean ~ k_nsl_wG4$pos, span=0.075)$fitted, type="n", ylim=c(0,1.35))
  abline(h=c(0.45,0.9))
  axis(side=2, at=c(c(0.8, 0.9, 1,1.1)-0.75,c(1, 1.1, 1.2)-0.4,c(1.1,1.2,1.3,1.4)-0.05), labels=c(c(0.8, 0.9, 1,1.1),c(1, 1.1, 1.2),c(1.1,1.2,1.3,1.4)), las=1)
  # abline(h=c(1,1-0.3,1-0.7), lty=3)
  lines(g_nsl_wG4$pos, loess(g_nsl_wG4$mean/g_nslshuf$mean ~ g_nsl_wG4$pos, span=0.075)$fitted-0.75, col="dark green", lwd=3)
  lines(g_nsg_wG4$pos, loess(g_nsg_wG4$mean/g_nsgshuf$mean ~ g_nsg_wG4$pos, span=0.075)$fitted-0.4, col="dark blue", lwd=2)
  lines(g_lg0_wG4$pos, loess(g_lg0_wG4$mean/g_nsgshuf$mean ~ g_lg0_wG4$pos, span=0.075)$fitted-0.05, col="dark cyan", lwd=2)
  
}

gen3Stack.mean <- function(){
  plot(ylab="Nucleosome Signal", xlab="Relative Position from Summit (bp)", yaxt="n",las=1,k_nsl_wG4$pos, loess(k_nsl_wG4$mean/k_nslshuf$mean ~ k_nsl_wG4$pos, span=0.075)$fitted, type="n", ylim=c(0,1.35))
  abline(h=c(0.45,0.9))
  axis(side=2, at=c(c(0.8, 0.9, 1, 1.1)-0.75,c(0.9,1, 1.1, 1.2)-0.35,c(1.0, 1.1,1.2,1.3)), labels=c(c(0.8, 0.9, 1,1.1),c(0.9, 1, 1.1, 1.2),c(1.0,1.1,1.2,1.3)), las=1)
  # abline(h=c(1,1-0.3,1-0.7), lty=3)
  nsl <- ((g_nsl_wG4$mean/g_nslshuf$mean) + (k_nsl_wG4$mean/k_nslshuf$mean))/2
  nsg <- ((g_nsg_wG4$mean/g_nsgshuf$mean) + (k_nsg_wG4$mean/k_nsgshuf$mean))/2
  lg0 <- ((g_lg0_wG4$mean/g_lg0shuf$mean) + (k_lg0_wG4$mean/k_lg0shuf$mean))/2
  lines(g_nsl_wG4$pos, loess(nsl ~ g_nsl_wG4$pos, span=0.075)$fitted-0.75, col="dark green", lwd=3)
  lines(g_nsg_wG4$pos, loess(nsg ~ g_nsg_wG4$pos, span=0.075)$fitted-0.35, col="dark blue", lwd=2)
  lines(g_lg0_wG4$pos, loess(lg0 ~ g_lg0_wG4$pos, span=0.075)$fitted, col="dark cyan", lwd=2)
  
}


gen3Stack.all <- function(g4=FALSE, g4score="counts", spacinglinesonlybehind=F, spacinglinesonlyfront=F){
  par(mar=c(5,5,5,5))
  #g4Score can be counts or FE
  plot(ylab="Nucleosome Signal", xlab="Relative Position from Summit (bp)", yaxt="n",las=1,k_nsl_wG4$pos, loess(k_nsl_wG4$mean/k_nslshuf$mean ~ k_nsl_wG4$pos, span=0.075)$fitted, type="n", ylim=c(0,1.35))
  abline(h=c(0.45,0.9))
  axis(side=2, at=c(c(0.8, 0.9, 1, 1.1)-0.79,c(0.9,1, 1.1, 1.2)-0.35,c(1.0, 1.1,1.2,1.3)), labels=c(c(0.8, 0.9, 1,1.1),c(0.9, 1, 1.1, 1.2),c(1.0,1.1,1.2,1.3)), las=1)
  # abline(h=c(1,1-0.3,1-0.7), lty=3)
  if(spacinglinesonlybehind){
    gen3stack.addg4spacing()
    gen3stack.addnsomespacing(plotarrows=F, plotlengths=F)
  }
  if(g4){gen3stack.addG4signal(plotaxis=T, g4score=g4score)}
  nsl <- ((g_nsl_wG4$mean/g_nslshuf$mean) + (k_nsl_wG4$mean/k_nslshuf$mean))/2
  nsg <- ((g_nsg_wG4$mean/g_nsgshuf$mean) + (k_nsg_wG4$mean/k_nsgshuf$mean))/2
  lg0 <- ((g_lg0_wG4$mean/g_lg0shuf$mean) + (k_lg0_wG4$mean/k_lg0shuf$mean))/2
  
  lines(g_nsl_wG4$pos, loess(nsl ~ g_nsl_wG4$pos, span=0.075)$fitted-0.79, col="black", lwd=3)
  lines(k_nsl_wG4$pos, loess(k_nsl_wG4$mean/k_nslshuf$mean ~ k_nsl_wG4$pos, span=0.075)$fitted-0.79, col="dark green", lwd=3)
  lines(g_nsl_wG4$pos, loess(g_nsl_wG4$mean/g_nslshuf$mean ~ g_nsl_wG4$pos, span=0.075)$fitted-0.79, col="dark green", lwd=3)
  
  lines(g_nsg_wG4$pos, loess(nsg ~ g_nsg_wG4$pos, span=0.075)$fitted-0.35, col="black", lwd=2)
  lines(k_nsg_wG4$pos, loess(k_nsg_wG4$mean/k_nsgshuf$mean ~ k_nsg_wG4$pos, span=0.075)$fitted-0.35, col="dark blue", lwd=2)
  lines(g_nsg_wG4$pos, loess(g_nsg_wG4$mean/g_nsgshuf$mean ~ g_nsg_wG4$pos, span=0.075)$fitted-0.35, col="dark blue", lwd=2)
  
  lines(g_lg0_wG4$pos, loess(lg0 ~ g_lg0_wG4$pos, span=0.075)$fitted, col="black", lwd=2)
  lines(k_lg0_wG4$pos, loess(k_lg0_wG4$mean/k_nsgshuf$mean ~ k_lg0_wG4$pos, span=0.075)$fitted, col="dark cyan", lwd=2)
  lines(g_lg0_wG4$pos, loess(g_lg0_wG4$mean/g_nsgshuf$mean ~ g_lg0_wG4$pos, span=0.075)$fitted, col="dark cyan", lwd=2)
  
  if(spacinglinesonlyfront){
    gen3stack.addg4spacing()
    gen3stack.addnsomespacing(plotarrows=F, plotlengths=F)
  }
}

gen3stack.addnsomespacing <- function(plotarrows=TRUE, plotlengths=TRUE){
  lg0 <- addNsomeLines("lg0", add=F)
  print(lg0)
  nsg <- addNsomeLines("nsg", add=F)
  print(nsg)
  nsl <- addNsomeLines("nsl", add=F)
  print(nsl)
#   abline(v=summits, lty=3, lwd=1.5)
  gen3stack.addSpacingLinesFromSpacingLinesObject(lg0, y0=0.9, y1=1.55, yArrow=1.3, yText=1.35, textcex=1, plotarrows=plotarrows, plotlengths=plotlengths)
  gen3stack.addSpacingLinesFromSpacingLinesObject(nsg, y0=0.45, y1=0.9, yArrow=0.8, yText=0.85, textcex=1, plotarrows=plotarrows, plotlengths=plotlengths)
  gen3stack.addSpacingLinesFromSpacingLinesObject(nsl, y0=-1, y1=0.45, yArrow=0.35, yText=0.4, textcex=1, plotarrows=plotarrows, plotlengths=plotlengths)
}

gen3stack.addSpacingLinesFromSpacingLinesObject <- function(spacingobject, y0, y1, yArrow, yText, textcex, plotlines=TRUE, plotarrows=TRUE, plotlengths=TRUE, linecol="black", lwd=0.75,lty=1, bold=F){
  font <- 1
  if(bold){font <- 2}
  ## nsomelinesobject is output of addNsomeLines() with add=F
  data <- spacingobject
  #add vertical lines at apexes
  if(plotlines){
    for(i in 1:length(data$apexes)){segments(x0=data$apexes[i], y0=y0, x1=data$apexes[i], y1=y1, col=linecol, lty=lty, lwd=lwd)}
  }
  ## construct arrows
  if(plotarrows){
    for(i in 1:(length(data$apexes)-1)){arrows(x0=data$apexes[i],y0=yArrow,x1=data$apexes[i+1],y1=yArrow, length=0.1, lwd=0.75)}
    for(i in 1:(length(data$apexes)-1)){arrows(x0=data$apexes[i+1],y0=yArrow,x1=data$apexes[i],y1=yArrow, length=0.1, lwd=0.75)}
  }
  ## add lengths above arrows
  if(plotlengths){
    if(length(textcex) != length(data$apexes)-1){textcex <- rep(textcex,times=(length(data$apexes)-1))}
    for(i in 1:(length(data$apexes)-1)){text(x=mean(c(data$apexes[i],data$apexes[i+1])),y=yText,labels=data$interApexDistances[i], cex=textcex[i], font=font)
#     print(c(data$apexes[i],data$apexes[i+1], data$interApexDistances[i]))}
    }
  }
}

gen3stack.addg4spacing <- function(){
  ## add just bars with relative heights OR just lines
  print("hi")
  data <- loess(g4Data[["lg0"]][["G4"]]$counts ~ g4Data[["lg0"]][["G4"]]$mids, span=0.1)
  lg0 <- findG4lines(data$fitted, hlg0Shufg4$counts, data$x)
  print(lg0)
  #nsg
  data <- loess(g4Data[["nsg"]][["G4"]]$counts ~ g4Data[["nsg"]][["G4"]]$mids, span=0.1)
  nsg <- findG4lines(data$fitted, hnsgShufg4$counts, data$x)
  print(nsg)
  #nsl
  data <- loess(g4Data[["nsl"]][["G4"]]$counts ~ g4Data[["nsl"]][["G4"]]$mids, span=0.1)
  nsl <- findG4lines(data$fitted, hnslShufg4$counts, data$x)
  print(nsl)
  gen3stack.addSpacingLinesFromSpacingLinesObject(lg0, y0=0.9, y1=1.55, yArrow=1.3, yText=1.35, textcex=1, plotarrows=FALSE, plotlengths=FALSE, linecol="red")
  gen3stack.addSpacingLinesFromSpacingLinesObject(nsg, y0=0.45, y1=0.9, yArrow=0.8, yText=0.85, textcex=1, plotarrows=FALSE, plotlengths=FALSE, linecol="red")
  gen3stack.addSpacingLinesFromSpacingLinesObject(nsl, y0=-1, y1=0.45, yArrow=0.35, yText=0.4, textcex=1, plotarrows=FALSE, plotlengths=FALSE, linecol="red")
  
}

gen3stack.addnsomeAndG4distances <- function(textcexlg0=0.5, textcexnsg=0.5,textcexnsl=0.5, yTextlg0=1.35, yTextnsg=0.85, yTextnsl=0.4, bold=T){
#   textcexlg0=c(1,1,0.5,1,0.5,0.5,0.5,0.5,1,0.5,1,1)
#   textcexnsg=c(1,0.5,1,1,1,0.5,1,0.5,0.5,1,0.5,1,1,1,0.5,1)
#   textcexnsl=c(1,0.5,1,0.5,1,0.5,1,0.5,0.5,1,0.5,1,0.5,1,0.5,1)
  data <- loess(g4Data[["lg0"]][["G4"]]$counts ~ g4Data[["lg0"]][["G4"]]$mids, span=0.1)
  lg0g4 <- findG4lines(data$fitted, hlg0Shufg4$counts, data$x)
  lg0nsome <- addNsomeLines("lg0", add=F)
  apexes <- sort(c(round(lg0g4$apexes), lg0nsome$apexes))
  distances <- apexes[2:length(apexes)]-apexes[1:(length(apexes)-1)]
  object <- list(apexes=apexes, interApexDistances=distances)
  gen3stack.addSpacingLinesFromSpacingLinesObject(object, yText=yTextlg0, textcex=textcexlg0, plotlines=F, plotarrows=F, plotlengths=TRUE, bold=bold)

  data <- loess(g4Data[["nsg"]][["G4"]]$counts ~ g4Data[["nsg"]][["G4"]]$mids, span=0.1)
  nsgg4 <- findG4lines(data$fitted, hnsgShufg4$counts, data$x)
  nsgnsome <- addNsomeLines("nsg", add=F)
  apexes <- sort(c(round(nsgg4$apexes), nsgnsome$apexes))
  distances <- apexes[2:length(apexes)]-apexes[1:(length(apexes)-1)]
  object <- list(apexes=apexes, interApexDistances=distances)
  gen3stack.addSpacingLinesFromSpacingLinesObject(object, yText=yTextnsg, textcex=textcexnsg, plotlines=F, plotarrows=F, plotlengths=TRUE, bold=bold)
  
  data <- loess(g4Data[["nsl"]][["G4"]]$counts ~ g4Data[["nsl"]][["G4"]]$mids, span=0.1)
  nslg4 <- findG4lines(data$fitted, hnslShufg4$counts, data$x)
  nslnsome <- addNsomeLines("nsl", add=F)
  apexes <- sort(c(round(nslg4$apexes), nslnsome$apexes))
  distances <- apexes[2:length(apexes)]-apexes[1:(length(apexes)-1)]
  object <- list(apexes=apexes, interApexDistances=distances)
  gen3stack.addSpacingLinesFromSpacingLinesObject(object, yText=yTextnsl, textcex=textcexnsl, plotlines=F, plotarrows=F, plotlengths=TRUE, bold=bold)
  
}

gen3stack.addnsomeAndG4lines <- function(plotarrows=TRUE, plotlengths=TRUE){
  gen3stack.addnsomespacing(plotarrows=F, plotlengths=F)
  gen3stack.addg4spacing()
}

gen3stack.addABC <- function(){
  text(x=-1000,y=1.35,labels="A", cex=1.5)
  text(x=-1000,y=0.85,labels="B", cex=1.5)
  text(x=-1000,y=0.4,labels="C", cex=1.5)
}

gen3stack.findApexDistances <- function(){
  ## add just bars with relative heights OR just lines
  data <- loess(g4Data[["lg0"]][["G4"]]$counts ~ g4Data[["lg0"]][["G4"]]$mids, span=0.1)
  lg0.g4 <- findG4lines(data$fitted, hlg0Shufg4$counts, data$x)
  lg0.nsome <- addNsomeLines("lg0", add=F)
  lg0 <- gen3stack.apexPositionDistanceAnalysis(g4=lg0.g4$apexes, nsome=lg0.nsome$apexes)
  print(list(g4Apexes=lg0.g4$apexes, g4interApexDistances=lg0.g4$interApexDistances, nsomeApexes=lg0.nsome$apexes, nsomeinterApexDistances=lg0.nsome$interApexDistances, nsomeToG4Dist=lg0))
  
  #nsg
  data <- loess(g4Data[["nsg"]][["G4"]]$counts ~ g4Data[["nsg"]][["G4"]]$mids, span=0.1)
  nsg.g4 <- findG4lines(data$fitted, hnsgShufg4$counts, data$x)
  nsg.nsome <- addNsomeLines("nsg", add=F)
  nsg <- gen3stack.apexPositionDistanceAnalysis(g4=nsg.g4$apexes, nsome=nsg.nsome$apexes)
  print(list(g4Apexes=nsg.g4$apexes, g4interApexDistances=nsg.g4$interApexDistances, nsomeApexes=nsg.nsome$apexes, nsomeinterApexDistances=nsg.nsome$interApexDistances, nsomeToG4Dist=nsg))
  #nsl
  data <- loess(g4Data[["nsl"]][["G4"]]$counts ~ g4Data[["nsl"]][["G4"]]$mids, span=0.1)
  nsl.g4 <- findG4lines(data$fitted, hnslShufg4$counts, data$x)
  nsl.nsome <- addNsomeLines("nsl", add=F)
  nsl <- gen3stack.apexPositionDistanceAnalysis(g4=nsl.g4$apexes, nsome=nsl.nsome$apexes)
  print(list(g4Apexes=nsl.g4$apexes, g4interApexDistances=nsl.g4$interApexDistances, nsomeApexes=nsl.nsome$apexes, nsomeinterApexDistances=nsl.nsome$interApexDistances, nsomeToG4Dist=nsl))
  return(list(lg0=lg0, nsg=nsg, nsl=nsl))
}

gen3stack.apexPositionDistanceAnalysis <- function(g4, nsome){
  whichG4 <- c()
  whichNsome <- c()
  minDist <- c()
  if(length(g4) < length(nsome)){
    for (e in g4){
      whichG4 <- c(whichG4, e)
      whichNsome <- c(whichNsome, nsome[which.min(abs(e-nsome))])
      minDist <- c(minDist, min(abs(e-nsome)))    
    }
  } else {
    whichG4 <- c(whichG4, nsome[which.min(abs(e-g4))])
    whichNsome <- c(whichNsome, e)
    minDist <- c(minDist, min(abs(e-nsome)))
  }
  distances <- data.frame(g4Apex=whichG4, nsomeApex=whichNsome, distance=minDist)
  return(distances)
}

findG4lines <- function(vector, controlVector, indexes){
  apexes <- indexes[allsummits5(vector=vector, controlVector=controlVector, vectorAlreadyIsFE=F, maxWindowSearch=80, minWindowSearch=70, totalWinLenCutOff=150, FEcut=1.5)] # 80,70,150,1.5 
  interApexDistances <- apexes[2:length(apexes)]-apexes[1:(length(apexes)-1)]
  return(list(apexes=apexes, interApexDistances=interApexDistances))
}

gen3stack.addG4signal <- function(plotaxis=F, g4score="counts"){
  #lg0
  if(g4score == "counts"){counts=T; FE=F} 
  else if (g4score == "FE"){FE=T; counts=F}
  if(counts){
    g4 <- loess(g4Data[["lg0"]][["G4"]]$counts ~ g4Data[["lg0"]][["G4"]]$mids, span=0.1)
    gen3stack.getG4signal(g4=g4, sampbase=10, sampadd=-0.8, plotaxis=plotaxis)
    #nsg
    g4 <- loess(g4Data[["nsg"]][["G4"]]$counts ~ g4Data[["nsg"]][["G4"]]$mids, span=0.1)
    gen3stack.getG4signal(g4=g4, sampbase=10, sampadd=-1.2, plotaxis=plotaxis)
    #nsl
    g4 <- loess(g4Data[["nsl"]][["G4"]]$counts ~ g4Data[["nsl"]][["G4"]]$mids, span=0.1)
    gen3stack.getG4signal(g4=g4, sampbase=100, sampadd=-0.4, plotaxis=plotaxis)
    #
    gen3stack.addG4axislabel("counts")
  } else if(FE){
    g4 <- loess(g4Data[["lg0"]][["G4"]]$counts/g4Data[["lg0"]][["shufG4"]]$counts ~ g4Data[["lg0"]][["G4"]]$mids, span=0.1)
    gen3stack.getG4signal(g4=g4, sampbase=10, sampadd=0.65, plotaxis=plotaxis)
    #nsg
    g4 <- loess(g4Data[["nsg"]][["G4"]]$counts/g4Data[["nsg"]][["shufG4"]]$counts ~ g4Data[["nsg"]][["G4"]]$mids, span=0.1)
    gen3stack.getG4signal(g4=g4, sampbase=10, sampadd=0.17, plotaxis=plotaxis)
    #nsl
    g4 <- loess(g4Data[["nsl"]][["G4"]]$counts/g4Data[["nsl"]][["shufG4"]]$counts ~ g4Data[["nsl"]][["G4"]]$mids, span=0.1)
    gen3stack.getG4signal(g4=g4, sampbase=100, sampadd=0.02, plotaxis=plotaxis)
    #
    gen3stack.addG4axislabel("FE")
  }
}


gen3stack.getG4signal <- function(g4, sampbase, sampadd, plotaxis=F){
  #g4 is loess smoothed line
  f <- function(x){log(x, base=sampbase)+sampadd}
  g <- function(y){sampbase^(y-sampadd)}
#   lines(g4$x, log(g4$fitted, base=sampbase)+sampadd, col="grey", lwd=2)
  lines(g4$x, f(g4$fitted), col="grey", lwd=2)
  print(c("G4 min/max:", min(g4$fitted), max(g4$fitted)))
  if (plotaxis){
    r <- range(f(g4$fitted))
    d <- 0.1*(r[2]-r[1])
    at <- seq(r[1]+d, r[2]-d, length.out=4)
    labels <- round(g(at), digits=0)
    at <- f(labels)
    axis(side=4, at=at, labels=labels, las=2)
  }
}

gen3stack.addG4axislabel <- function(countsOrFE="counts"){
  if(countsOrFE == "counts"){
    mtext("G4 Counts", side=4, line=3) 
  } else if(countsOrFE == "FE"){
    mtext("G4 Fold Enrichment", side=4, line=3) 
  }
}
gen3stack.NsomeG4Correlations <- function(ALL=F, LG0=F, NSG=F, NSL=F, SUMSQRERR=T){
  ## sample = all, lg0, nsg, nsl
  if(ALL){LG0=T; NSG=T; NSL=T}
  if(LG0){
    print("LG0")
    lg0.k.raw <- k_lg0_wG4$mean/k_lg0shuf$mean
    lg0.g.raw <- g_lg0_wG4$mean/g_lg0shuf$mean
    lg0 <- (lg0.k.raw + lg0.g.raw)/2
    lg0.mean <- loess(lg0 ~ g_lg0_wG4$pos, span=0.075)
    lg0.k <- loess(lg0.k.raw ~ k_lg0_wG4$pos, span=0.075)
    lg0.g <- loess(lg0.g.raw ~ g_lg0_wG4$pos, span=0.075)
    printNsomeG4Correlations(k=lg0.k, g=lg0.g, cell.line.mean=lg0.mean, g4=llg0g4, k.raw=lg0.k.raw, g.raw=lg0.g.raw)
    if(SUMSQRERR){
      k.error <- sum((lg0.k$fitted - lg0.mean$fitted)^2)
      g.error <- sum((lg0.g$fitted - lg0.mean$fitted)^2)
      total.error <- k.error+g.error
      print(c("Sum of Squared Devs from Means for K562", k.error))
      print(c("Sum of Squared Devs from Means for GM12878", g.error))
      print(c("Sum of Squared Devs from Means for Total", total.error))
      print(c("Sum of Squared Diffs is", sum((lg0.k$fitted-lg0.g$fitted)^2)))
      print(c("Area between cell lines (sum of abs diffs)", sum(abs(lg0.k$fitted-lg0.g$fitted))))
    }
  }
  print("---")
  if(NSG){
    print("NSG")
    nsg.k.raw <- k_nsg_wG4$mean/k_nsgshuf$mean
    nsg.g.raw <- g_nsg_wG4$mean/g_nsgshuf$mean
    nsg <- (nsg.k.raw + nsg.g.raw)/2
    nsg.mean <- loess(nsg ~ g_nsg_wG4$pos, span=0.075)
    nsg.k <- loess(nsg.k.raw ~ k_nsg_wG4$pos, span=0.075)
    nsg.g <- loess(nsg.g.raw ~ g_nsg_wG4$pos, span=0.075)
    printNsomeG4Correlations(k=nsg.k, g=nsg.g, cell.line.mean=nsg.mean, g4=lnsgg4, k.raw=nsg.k.raw, g.raw=nsg.g.raw)
    if(SUMSQRERR){
      k.error <- sum((nsg.k$fitted - nsg.mean$fitted)^2)
      g.error <- sum((nsg.g$fitted - nsg.mean$fitted)^2)
      total.error <- k.error+g.error
      print(c("Sum of Squared Devs from Means for K562", k.error))
      print(c("Sum of Squared Devs from Means for GM12878", g.error))
      print(c("Sum of Squared Devs from Means for Total", total.error))
      print(c("Sum of Squared Diffs is", sum((nsg.k$fitted-nsg.g$fitted)^2)))
      print(c("Area between cell lines (sum of abs diffs)", sum(abs(nsg.k$fitted-nsg.g$fitted))))
    }
  }
  print("---")
  if(NSL){
    print("NSL")
    nsl.k.raw <- k_nsl_wG4$mean/k_nslshuf$mean 
    nsl.g.raw <- g_nsl_wG4$mean/g_nslshuf$mean
    nsl <- (nsl.k.raw + nsl.g.raw)/2
    nsl.mean <- loess(nsl ~ g_nsl_wG4$pos, span=0.075)
    nsl.k <- loess(nsl.k.raw ~ k_nsl_wG4$pos, span=0.075)
    nsl.g <- loess(nsl.g.raw ~ g_nsl_wG4$pos, span=0.075)
    printNsomeG4Correlations(k=nsl.k, g=nsl.g, cell.line.mean=nsl.mean, g4=lnslg4, k.raw=nsl.k.raw, g.raw=nsl.g.raw)
    if(SUMSQRERR){
      k.error <- sum((nsl.k$fitted - nsl.mean$fitted)^2)
      g.error <- sum((nsl.g$fitted - nsl.mean$fitted)^2)
      total.error <- k.error+g.error
      print(c("Sum of Squared Devs from Means for K562", k.error))
      print(c("Sum of Squared Devs from Means for GM12878", g.error))
      print(c("Sum of Squared Devs from Means for Total", total.error))
      print(c("Sum of Squared Diffs is", sum((nsl.k$fitted-nsl.g$fitted)^2)))
      print(c("Area between cell lines (sum of abs diffs)", sum(abs(nsl.k$fitted-nsl.g$fitted))))
    }
  }
}

printNsomeG4Correlations <- function(k, g, cell.line.mean, g4, k.raw, g.raw){
  print(c("mean v. g4 pearson/spearman:", cor(x=g4$fitted, cell.line.mean$fitted[1:2000]), cor(x=g4$fitted, cell.line.mean$fitted[1:2000], method="spearman")))
  print(c("k562 v. g4 pearson/spearman:", cor(x=g4$fitted, k$fitted[1:2000]), cor(x=g4$fitted, k$fitted[1:2000], method="spearman")))
  print(c("gm12878 v. g4 pearson/spearman:", cor(x=g4$fitted, g$fitted[1:2000]), cor(x=g4$fitted, g$fitted[1:2000], method="spearman")))
  print(c("k562 v. gm12878 smooth pear/spear:", cor(x=k$fitted, g$fitted), cor(x=k$fitted, g$fitted, method="spearman")))
  print(c("k562 v. mean pearson/spearman:", cor(x=k$fitted, cell.line.mean$fitted), cor(x=k$fitted, cell.line.mean$fitted, method="spearman")))
  print(c("gm12878 v. mean pearson/spearman:", cor(x=g$fitted, cell.line.mean$fitted), cor(x=g$fitted, cell.line.mean$fitted, method="spearman")))
  print(c("k562 v. gm12878 raw pear/spear:", cor(x=k.raw, g.raw), cor(x=k.raw, g.raw, method="spearman")))
  
#   print(c("mean vs. g4 pearson/spearman:", cor(x=g4$fitted, cell.line.mean$fitted[1:2000]), cor(x=g4$fitted, cell.line.mean$fitted[1:2000], method="spearman")))
#   print(c("k562 vs. gm12878 pearson:", cor(x=k$fitted, g$fitted), cor(x=k$fitted, g$fitted, method="spearman")))
#   print(c("k562 vs. gm12878 pearson/spearman:", cor(x=k$fitted, g$fitted), cor(x=k$fitted, g$fitted, method="spearman"))) 
}



neighboringSinWaveModel  <- function(){
  plot(seq(0,40,0.1),sin(seq(0,40,0.1)), type='l', ylab="Amplitude", xlab="distance", lwd=2, ylim=c(-4,4), las=1, yaxt='n')
  axis(side=2, at=(-4:4), las=1)
  lines(seq(0,40,0.1), sin(seq(0,40,0.1)+1.6), col="blue")
  print(c("Pearson/Spearman Cor =", cor(sin(seq(0,40,0.1)), sin(seq(0,40,0.1)+1.6)), cor(sin(seq(0,40,0.1)), sin(seq(0,40,0.1)+1.6), method="spearman")))
}

preliminaryNsomesAroundSummitsWithG4Analysis <- function(sample="lg0", K562=FALSE, GM12878=FALSE, MEAN=TRUE, plotG4=TRUE,colstyle=1, FEoverGenomicAvg=TRUE){
  ## colstyle can be 1,2,3
  ## if FEoverGenomicAvg=TRUE, then cell lines nsome scores at each position are normalized by means over each position from all shuffled peaks
  ##  The mean is then (normK+normG)/2 rather than (K+G)/2
  ##  The normalization happens before smoothing
  ##  The mean also happens before smoothing
  span=0.075
  if(sample == "lg0"){
    k562 <- k_lg0_wG4
    gm12878 <- g_lg0_wG4
    if (FEoverGenomicAvg){
      k562$mean <- k562$mean/k_lg0shuf$mean
      gm12878$mean <- gm12878$mean/g_lg0shuf$mean
    }
    l <- loess((k562$mean+gm12878$mean)/2 ~ gm12878$pos, span=span)
    color <- "dark cyan"
#     ylim <- c(0.9,1.3) ## style1
#     ylim <- c(0.95,1.25)
    g4 <- llg0g4
    sampbase <- 100
    sampadd <- 0.85
    if(MEAN & !(K562 | GM12878)){ylim <- c(1,1.3); sampbase <- 100; sampadd <- 0.85}
    else if((MEAN & K562 & GM12878) | (K562 & GM12878)){ylim <- c(0.9,1.3); sampbase <- 10; sampadd <- 0.55}
    else if(K562 & !(MEAN | GM12878)){ylim <- c(0.9,1.15); sampbase <- 100; sampadd <- 0.85}
    else if(GM12878 & !(MEAN | K562)){ylim <- c(1.1,1.3); sampbase <- 100; sampadd <- 0.85}
    else if((K562 & MEAN) & !(GM12878)){ylim <- c(0.9,1.2); sampbase <- 100; sampadd <- 0.85}
    else if((GM12878 & MEAN) & !(K562)){ylim <- c(1.05,1.3); sampbase <- 100; sampadd <- 0.85}
  } else if (sample == "nsg"){
    k562 <- k_nsg_wG4
    gm12878 <- g_nsg_wG4
    if (FEoverGenomicAvg){
      k562$mean <- k562$mean/k_nsgshuf$mean
      gm12878$mean <- gm12878$mean/g_nsgshuf$mean
    }
    l <- loess((k562$mean+gm12878$mean)/2 ~ gm12878$pos, span=span)
    color <- "dark blue"
#     ylim <- c(0.93,1.15) ## style 1
    ylim <- c(0.93,1.05)
    ylim <- c(0.93,1.1)
    g4 <- lnsgg4
    sampbase <- 10000
    sampadd <- 0.85
    if(MEAN & !(K562 | GM12878)){ylim <- c(0.93,1.1); sampbase <- 10000; sampadd <- 0.85}
    else if((MEAN & K562 & GM12878) | (K562 & GM12878)){ylim <- c(0.8,1.2); sampbase <- 10000; sampadd <- 0.85}
    else if(K562 & !(MEAN | GM12878)){ylim <- c(0.8,1.15); sampbase <- 10000; sampadd <- 0.85}
    else if(GM12878 & !(MEAN | K562)){ylim <- c(1.0,1.2); sampbase <- 10000; sampadd <- 0.85}
    else if((K562 & MEAN) & !(GM12878)){ylim <- c(0.8,1.2); sampbase <- 10000; sampadd <- 0.85}
    else if((GM12878 & MEAN) & !(K562)){ylim <- c(0.93,1.2); sampbase <- 10000; sampadd <- 0.85}
  } else if (sample == "nsl"){
    k562 <- k_nsl_wG4
    gm12878 <- g_nsl_wG4
    if (FEoverGenomicAvg){
      k562$mean <- k562$mean/k_nslshuf$mean
      gm12878$mean <- gm12878$mean/g_nslshuf$mean
    }
    l <- loess((k562$mean+gm12878$mean)/2 ~ gm12878$pos, span=span)
    color <- "dark green"
#     ylim <- c(0.6,1.3) ## style1 with sampbase=500,sampadd=0.75, ylim=c(0.6,1.3)
    ylim <- c(0.6,1.15)
    ylim <- c(0.75,1.2)
    ylim <- c(0.6,1.25)
    if(MEAN & !(K562 | GM12878)){ylim <- c(0.6,1.25); ylim <- c(0,1.5);}
    else if((MEAN & K562 & GM12878) | (K562 & GM12878)){ylim <- c(0.6,1.25)}
    else if(K562 & !(MEAN | GM12878)){ylim <- c(0.6,1.25)}
    else if(GM12878 & !(MEAN | K562)){ylim <- c(0.6,1.25)}
    else if((K562 & MEAN) & !(GM12878)){ylim <- c(0.6,1.25)}
    else if((GM12878 & MEAN) & !(K562)){ylim <- c(0.6,1.25)}
    g4 <- lnslg4
    sampbase <- 100
    sampadd <- 0.7
  }
#   g4color <- color
#   color <- "black"
  if (colstyle == 1){
    g4color <- "red"
    if(MEAN & (K562 | GM12878)){meancolor <- "black"}
    else{meancolor <- color}
  } else if (colstyle == 2){
    g4color <- "grey"
    meancolor <- color
  } else if (colstyle == 3){
    g4color <- "red"
    meancolor <- "black"
    } else if (colstyle == 4){
      g4color <- color
      color <- "black"
      meancolor <- color
  }
  ###### PLOTTING!! ####
  par(mar=c(5,5,5,5))
  plot(k_lg0_wG4$pos, k_lg0_wG4$mean, col="dark cyan", las=1, ylim=ylim, type='n', ylab="Average Nucleosome Signal", xlab="Relative Position (bp)")
#   abline(v=0, lty=2)
  if(GM12878){lines(gm12878$pos, loess(gm12878$mean ~ gm12878$pos, span=span)$fitted, col=color, lwd=2)} 
  if (K562){lines(k562$pos, loess(k562$mean ~ k562$pos, span=span)$fitted, col=color, lwd=2)}
  if (MEAN){lines(l$x, l$fitted, col=meancolor, lwd=2)}
  if(plotG4){
#     lines(g4$x, 0.3+g4$fitted/max(g4$fitted), col="grey", lwd=2)
    lines(g4$x, log(g4$fitted, base=sampbase)+sampadd, col=g4color, lwd=2)
    f <- function(x){log(x, base=sampbase)+sampadd}
    g <- function(y){sampbase^(y-sampadd)}
#     print(max(g4$fitted)); print(min(g4$fitted))
    at <- f(c(0.5,0.7,1,1.5,2,2.5,3:6, 8, 10))
#     print(at)
    labels <- round(g(at), digits=1)
#     print(labels)
    at <- f(labels)
#     print(at)
    mtext("G4 Fold Enrichment", side=4, line=3)
    axis(side=4, at=at, labels=labels, las=2)
  }
  k <- loess(k562$mean ~ k562$pos, span=span)
  g <- loess(gm12878$mean ~ gm12878$pos, span=span)
  print(c("mean v. g4 pearson/spearman:", cor(x=g4$fitted, l$fitted[1:2000]), cor(x=g4$fitted, l$fitted[1:2000], method="spearman")))
  print(c("k562 v. gm12878 smooth pear/spear:", cor(x=k$fitted, g$fitted), cor(x=k$fitted, g$fitted, method="spearman")))
  print(c("k562 v. mean pearson/spearman:", cor(x=k$fitted, l$fitted), cor(x=k$fitted, l$fitted, method="spearman")))
  print(c("gm12878 v. mean pearson/spearman:", cor(x=g$fitted, l$fitted), cor(x=g$fitted, l$fitted, method="spearman")))
  print(c("k562 v. gm12878 raw pear/spear:", cor(x=k562$mean, gm12878$mean), cor(x=k562$mean, gm12878$mean, method="spearman")))
  
}


preliminaryNsomesAroundSummitsWithOutG4Analysis <- function(sample="lg0", K562=FALSE, GM12878=FALSE, MEAN=TRUE, plotG4=FALSE){
  span=0.075
  if(sample == "lg0"){
    k562 <- k_lg0_woG4
    gm12878 <- g_lg0_woG4
    l <- loess((k_lg0_woG4$mean+g_lg0_woG4$mean)/2 ~ g_lg0_woG4$pos, span=0.05)
    color <- "dark cyan"
    ylim <- c(1.1,1.9)
    g4 <- llg0g4
    sampbase <- 10
    sampadd <- 0.9
  } else if (sample == "nsg"){
    k562 <- k_nsg_woG4
    gm12878 <- g_nsg_woG4
    l <- loess((k_nsg_woG4$mean+g_nsg_woG4$mean)/2 ~ g_nsg_woG4$pos, span=0.05)
    color <- "dark blue"
    ylim <- c(0.9,1.9)
    g4 <- lnsgg4
    sampbase <- 10
    sampadd <- 0.8
  } else if (sample == "nsl"){
    k562 <- k_nsl_woG4
    gm12878 <- g_nsl_woG4
    l <- loess((k_nsl_woG4$mean+g_nsl_woG4$mean)/2 ~ g_nsl_woG4$pos, span=0.05)
    color <- "dark green"
    ylim <- c(0.9,2.9)
    g4 <- lnslg4
    sampbase <- 100
    sampadd <- 1.2
  }
  meancolor <- color
  plot(k_lg0_woG4$pos, k_lg0_woG4$mean, col="dark cyan", ylim=ylim, type='n', ylab="Average Nucleosome Signal", xlab="Relative Position (bp)")
#   abline(v=0, lty=2)
  if(GM12878){lines(gm12878$pos, gm12878$mean, col=color, lwd=2)} 
  if (K562){lines(k562$pos, k562$mean, col=color, lwd=2)}
  if (MEAN){lines(l$x, l$fitted, col=meancolor, lwd=2)}
  if(plotG4){lines(g4$x, log(g4$fitted, base=sampbase)+sampadd, col="red", lwd=2)}
  k <- loess(k562$mean ~ k562$pos, span=span)
  g <- loess(gm12878$mean ~ gm12878$pos, span=span)
  print(c("mean vs. g4 pearson/spearman:", cor(x=g4$fitted, l$fitted[1:2000]), cor(x=g4$fitted, l$fitted[1:2000], method="spearman")))
  print(c("k562 vs. gm12878 pearson:", cor(x=k$fitted, g$fitted), cor(x=k$fitted, g$fitted, method="spearman")))
  print(c("k562 vs. gm12878 pearson/spearman:", cor(x=k$fitted, g$fitted), cor(x=k$fitted, g$fitted, method="spearman")))
  
}




addNsomeLines <- function(sample="nsl",yArrow=NA,yText=NA, textcex=NA, add=TRUE){
#   print("What do nsomes look like from local shuffle?")
#   print("What do nsomes look like when select G4 mids randomly?")
  
  if (sample == "nsl"){
    k <- k_nsl_wG4$mean/k_nslshuf$mean
    g <- g_nsl_wG4$mean/g_nslshuf$mean
    x <- g_nsl_wG4
    if(is.na(yArrow)){yArrow <- 1.15}
    if(is.na(yText)){yText <- 1.18}
    if(is.na(textcex)){textcex <- 0.75}
  } else if(sample == "nsg"){
    k <- k_nsg_wG4$mean/k_nsgshuf$mean
    g <- g_nsg_wG4$mean/g_nsgshuf$mean
    x <- g_nsg_wG4
    if(is.na(yArrow)){yArrow <- 1.08}
    if(is.na(yText)){yText <- 1.09}
    if(is.na(textcex)){textcex <- 0.75}
  } else if(sample == "lg0"){
    k <- k_lg0_wG4$mean/k_lg0shuf$mean
    g <- g_lg0_wG4$mean/g_lg0shuf$mean
    x <- g_lg0_wG4
    if(is.na(yArrow)){yArrow <- 1.25}
    if(is.na(yText)){yText <- 1.26}
    if(is.na(textcex)){textcex <- 0.75}
  }  
  rawmean <- (k+g)/2
  l <- loess(rawmean ~ x$pos, span=0.075)
  summits <- l$x[allsummits3(l$fitted, maxWindowSearch=80, minWindowSearch=70, totalWinLenCutOff=130)]
  summits <- l$x[allapexes(l$fitted, rawmean, maxWindowSearch=80, minWindowSearch=50, totalWinLenCutOff=130, pValCut=0.1)]
  InterLens <- summits[2:length(summits)]-summits[1:(length(summits)-1)]
  if(add){
    abline(v=summits, lty=3, lwd=1.5)
    for(i in 1:(length(summits)-1)){arrows(x0=summits[i],y0=yArrow,x1=summits[i+1],y1=yArrow, length=0.1, lwd=2)}
    for(i in 1:(length(summits)-1)){arrows(x0=summits[i+1],y0=yArrow,x1=summits[i],y1=yArrow, length=0.1, lwd=2)}
    for(i in 1:(length(summits)-1)){text(x=mean(c(summits[i],summits[i+1])),y=yText,labels=InterLens[i], cex=textcex)}
  }
  return(list(apexes=summits, interApexDistances=InterLens))
}
#   } else if(sample == "nsg"){
#     
#     summits <- l$x[allsummits3(l$fitted, maxWindowSearch=80, minWindowSearch=30, totalWinLenCutOff=130)]
#     print(summits)
#     abline(v=summits)
#     InterLens <- summits[2:length(summits)]-summits[1:(length(summits)-1)]
#     for(i in 1:(length(summits)-1)){arrows(x0=summits[i],y0=yArrow,x1=summits[i+1],y1=yArrow, length=0.1, lwd=2)}
#     for(i in 1:(length(summits)-1)){arrows(x0=summits[i+1],y0=yArrow,x1=summits[i],y1=yArrow, length=0.1, lwd=2)}
#     for(i in 1:(length(summits)-1)){text(x=mean(c(summits[i],summits[i+1])),y=yText,labels=InterLens[i], cex=textcex)}
#   } else if(sample == "lg0"){
#     l <- loess((k_lg0_wG4$mean+g_lg0_wG4$mean)/2 ~ g_lg0_wG4$pos, span=0.05)
#     summits <- l$x[allsummits3(l$fitted, maxWindowSearch=80, minWindowSearch=30, totalWinLenCutOff=130)]
#     print(summits)
#     abline(v=summits)
#     InterLens <- summits[2:length(summits)]-summits[1:(length(summits)-1)]
#     for(i in 1:(length(summits)-1)){arrows(x0=summits[i],y0=yArrow,x1=summits[i+1],y1=yArrow, length=0.1, lwd=2)}
#     for(i in 1:(length(summits)-1)){arrows(x0=summits[i+1],y0=yArrow,x1=summits[i],y1=yArrow, length=0.1, lwd=2)}
#     for(i in 1:(length(summits)-1)){text(x=mean(c(summits[i],summits[i+1])),y=yText,labels=InterLens[i], cex=textcex)}
#   }
# }






 ### OBSOLETE -- older versions
# preliminaryNsomesAroundSummitsWithG4Analysis <- function(sample="lg0", bothAndMean=FALSE){
#   if(sample == "lg0"){
#     plot(k_lg0_wG4$pos, k_lg0_wG4$mean, col="dark cyan", ylim=c(0.9,1.3), type='n', ylab="Average Nucleosome Signal", xlab="Relative Position (bp)")
#     if(bothAndMean){
#       lines(g_lg0_wG4$pos, g_lg0_wG4$mean, col="dark cyan")
#       lines(k_lg0_wG4$pos, k_lg0_wG4$mean, col="dark cyan")
#     }
#     l <- loess((k_lg0_wG4$mean+g_lg0_wG4$mean)/2 ~ g_lg0_wG4$pos, span=0.05)
#     lines(l$x, l$fitted, col="black")
#     lines(llg0g4$x, log(llg0g4$fitted, base=100)+0.85, col="red")
#     print(c("pearson:", cor(x=llg0g4$fitted, l$fitted[1:2000])))
#     print(c("spearman:", cor(x=llg0g4$fitted, l$fitted[1:2000], method="spearman")))
#   } else if (sample == "nsg"){
#     plot(k_nsg_wG4$pos, k_nsg_wG4$mean, col="dark blue", ylim=c(0.8,1.3), type='n', ylab="Average Nucleosome Signal", xlab="Relative Position (bp)")
#     if(bothAndMean){
#       lines(g_nsg_wG4$pos, g_nsg_wG4$mean, col="dark blue")
#       lines(k_nsg_wG4$pos, k_nsg_wG4$mean, col="dark blue")
#     }
#     l <- loess((k_nsg_wG4$mean+g_nsg_wG4$mean)/2 ~ g_nsg_wG4$pos, span=0.05)
#     lines(l$x, l$fitted, col="black")
#     lines(lnsgg4$x, log(lnsgg4$fitted, base=100)+0.7, col="red")
#     print(c("pearson:", cor(x=lnsgg4$fitted, y=l$fitted[1:2000])))
#     print(c("spearman:", cor(x=lnsgg4$fitted, y=l$fitted[1:2000], method="spearman")))
#   } else if (sample == "nsl"){
#     plot(k_nsl_wG4$pos, k_nsl_wG4$mean, col="dark green", ylim=c(0.6,1.3), type='n', ylab="Average Nucleosome Signal", xlab="Relative Position (bp)")
#     if(bothAndMean){
#       lines(g_nsl_wG4$pos, g_nsl_wG4$mean, col="dark green")
#       lines(k_nsl_wG4$pos, k_nsl_wG4$mean, col="dark green")
#     }
#     l <- loess((k_nsl_wG4$mean+g_nsl_wG4$mean)/2 ~ g_nsl_wG4$pos, span=0.05)
#     lines(l$x, l$fitted, col="black")
#     lines(lnslg4$x, log(lnslg4$fitted, base=500)+0.75, col="red")
#     print(c("pearson", cor(x=lnslg4$fitted, y=l$fitted[1:2000])))
#     print(c("spearman", cor(x=lnslg4$fitted, y=l$fitted[1:2000], method="spearman")))
#   }
# }
# 
# 
# preliminaryNsomesAroundSummitsWithOutG4Analysis <- function(sample="lg0", bothAndMean=FALSE){
# if(sample == "lg0"){
#   plot(k_lg0_woG4$pos, k_lg0_woG4$mean, col="dark cyan", ylim=c(1.1,1.9), type='n', ylab="Average Nucleosome Signal", xlab="Relative Position (bp)")
#   if(bothAndMean){
#     lines(g_lg0_woG4$pos, g_lg0_woG4$mean, col="dark cyan")
#     lines(k_lg0_woG4$pos, k_lg0_woG4$mean, col="dark cyan")
#   }
#   l <- loess((k_lg0_woG4$mean+g_lg0_woG4$mean)/2 ~ g_lg0_woG4$pos, span=0.05)
#   lines(l$x, l$fitted, col="black")
#   lines(llg0g4$x, log(llg0g4$fitted, base=10)+0.9, col="red")
#   abline(v=0,lty=2)
#   print(c("pearson:", cor(x=llg0g4$fitted, l$fitted[1:2000])))
#   print(c("spearman:", cor(x=llg0g4$fitted, l$fitted[1:2000], method="spearman")))
# } else if (sample == "nsg"){
#   plot(k_nsg_woG4$pos, k_nsg_woG4$mean, col="dark blue", ylim=c(0.9,1.9), type='n', ylab="Average Nucleosome Signal", xlab="Relative Position (bp)")
#   if(bothAndMean){
#     lines(g_nsg_woG4$pos, g_nsg_woG4$mean, col="dark blue")
#     lines(k_nsg_woG4$pos, k_nsg_woG4$mean, col="dark blue")
#   }
#   l <- loess((k_nsg_woG4$mean+g_nsg_woG4$mean)/2 ~ g_nsg_woG4$pos, span=0.05)
#   lines(l$x, l$fitted, col="black")
#   abline(v=0,lty=2)
#   lines(lnsgg4$x, log(lnsgg4$fitted, base=10)+0.8, col="red")
#   print(c("pearson:", cor(x=lnsgg4$fitted, y=l$fitted[1:2000])))
#   print(c("spearman:", cor(x=lnsgg4$fitted, y=l$fitted[1:2000], method="spearman")))
# } else if (sample == "nsl"){
#   plot(k_nsl_woG4$pos, k_nsl_woG4$mean, col="dark green", ylim=c(0.9,2.9), type='n', ylab="Average Nucleosome Signal", xlab="Relative Position (bp)")
#   if(bothAndMean){
#     lines(k_nsl_woG4$pos, k_nsl_woG4$mean, col="dark green")
#     lines(g_nsl_woG4$pos, g_nsl_woG4$mean, col="dark green")
#   }
#   l <- loess((k_nsl_woG4$mean+g_nsl_woG4$mean)/2 ~ g_nsl_woG4$pos, span=0.05)
#   lines(l$x, l$fitted, col="black")
#   lines(lnslg4$x, log(lnslg4$fitted, base=100)+1.2, col="red")
#   print(c("pearson", cor(x=lnslg4$fitted, y=l$fitted[1:2000])))
#   print(c("spearman", cor(x=lnslg4$fitted, y=l$fitted[1:2000], method="spearman")))
#   abline(v=0,lty=2)
# }
# }

# 
# addNsomeLines <- function(sample="nsl",yArrow=1.2,yText=1.25, textcex=0.75){
#   print("What do nsomes look like from local shuffle?")
#   print("What do nsomes look like when select G4 mids randomly?")
#   if (sample == "nsl"){
#     l <- loess((k_nsl_wG4$mean+g_nsl_wG4$mean)/2 ~ g_nsl_wG4$pos, span=0.05)
#     summits <- l$x[allsummits3(l$fitted, maxWindowSearch=80, minWindowSearch=30, totalWinLenCutOff=130)]
#     print(summits)
#     abline(v=summits, lty=3, lwd=1.5)
#     InterLens <- summits[2:length(summits)]-summits[1:(length(summits)-1)]
#     for(i in 1:(length(summits)-1)){arrows(x0=summits[i],y0=yArrow,x1=summits[i+1],y1=yArrow, length=0.1, lwd=2)}
#     for(i in 1:(length(summits)-1)){arrows(x0=summits[i+1],y0=yArrow,x1=summits[i],y1=yArrow, length=0.1, lwd=2)}
#     for(i in 1:(length(summits)-1)){text(x=mean(c(summits[i],summits[i+1])),y=yText,labels=InterLens[i], cex=textcex)}
#     
#   } else if(sample == "nsg"){
#     l <- loess((k_nsg_wG4$mean+g_nsg_wG4$mean)/2 ~ g_nsg_wG4$pos, span=0.05)
#     summits <- l$x[allsummits3(l$fitted, maxWindowSearch=80, minWindowSearch=30, totalWinLenCutOff=130)]
#     print(summits)
#     abline(v=summits)
#     InterLens <- summits[2:length(summits)]-summits[1:(length(summits)-1)]
#     for(i in 1:(length(summits)-1)){arrows(x0=summits[i],y0=yArrow,x1=summits[i+1],y1=yArrow, length=0.1, lwd=2)}
#     for(i in 1:(length(summits)-1)){arrows(x0=summits[i+1],y0=yArrow,x1=summits[i],y1=yArrow, length=0.1, lwd=2)}
#     for(i in 1:(length(summits)-1)){text(x=mean(c(summits[i],summits[i+1])),y=yText,labels=InterLens[i], cex=textcex)}
#   } else if(sample == "lg0"){
#     l <- loess((k_lg0_wG4$mean+g_lg0_wG4$mean)/2 ~ g_lg0_wG4$pos, span=0.05)
#     summits <- l$x[allsummits3(l$fitted, maxWindowSearch=80, minWindowSearch=30, totalWinLenCutOff=130)]
#     print(summits)
#     abline(v=summits)
#     InterLens <- summits[2:length(summits)]-summits[1:(length(summits)-1)]
#     for(i in 1:(length(summits)-1)){arrows(x0=summits[i],y0=yArrow,x1=summits[i+1],y1=yArrow, length=0.1, lwd=2)}
#     for(i in 1:(length(summits)-1)){arrows(x0=summits[i+1],y0=yArrow,x1=summits[i],y1=yArrow, length=0.1, lwd=2)}
#     for(i in 1:(length(summits)-1)){text(x=mean(c(summits[i],summits[i+1])),y=yText,labels=InterLens[i], cex=textcex)}
#   }
# }


##CSHL2015
gen3Stack.cshl2015 <- function(g4=FALSE, g4score="counts", spacinglinesonlybehind=F, spacinglinesonlyfront=F){
  par(mar=c(5,5,5,5))
  #g4Score can be counts or FE
  plot(ylab="Nucleosome Signal", xlab="Relative Position from Summit (bp)", yaxt="n",las=1,k_nsl_wG4$pos, loess(k_nsl_wG4$mean/k_nslshuf$mean ~ k_nsl_wG4$pos, span=0.075)$fitted, type="n", ylim=c(0,1.35))
  abline(h=c(0.45,0.9))
  axis(side=2, at=c(c(0.8, 0.9, 1, 1.1)-0.79,c(0.9,1, 1.1, 1.2)-0.35,c(1.0, 1.1,1.2,1.3)), labels=c(c(0.8, 0.9, 1,1.1),c(0.9, 1, 1.1, 1.2),c(1.0,1.1,1.2,1.3)), las=1)
  # abline(h=c(1,1-0.3,1-0.7), lty=3)
  if(spacinglinesonlybehind){
    gen3stack.addg4spacing()
    gen3stack.addnsomespacing(plotarrows=F, plotlengths=F)
  }
  if(g4){gen3stack.addG4signal(plotaxis=T, g4score=g4score)}
  nsl <- ((g_nsl_wG4$mean/g_nslshuf$mean) + (k_nsl_wG4$mean/k_nslshuf$mean))/2
  nsg <- ((g_nsg_wG4$mean/g_nsgshuf$mean) + (k_nsg_wG4$mean/k_nsgshuf$mean))/2
  lg0 <- ((g_lg0_wG4$mean/g_lg0shuf$mean) + (k_lg0_wG4$mean/k_lg0shuf$mean))/2
  
  lines(g_nsl_wG4$pos, loess(nsl ~ g_nsl_wG4$pos, span=0.075)$fitted-0.79, col="dark green", lwd=3)
  
  lines(g_nsg_wG4$pos, loess(nsg ~ g_nsg_wG4$pos, span=0.075)$fitted-0.35, col="dark blue", lwd=2)
  
  lines(g_lg0_wG4$pos, loess(lg0 ~ g_lg0_wG4$pos, span=0.075)$fitted, col="dark cyan", lwd=2)
  
  if(spacinglinesonlyfront){
    gen3stack.addg4spacing()
    gen3stack.addnsomespacing(plotarrows=F, plotlengths=F)
  }
}
