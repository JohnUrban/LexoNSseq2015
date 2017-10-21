args <- commandArgs(trailingOnly=TRUE)
## 1=input 
## 2=method
## 3=outputpre

library(lattice)

df <- read.table(args[1], header=TRUE)

z<-cor(df, method=args[2]) 

write.table(x=z, file=paste0(args[3],'.csv'), quote=FALSE, na='', row.names=TRUE, col.names=NA, sep=',')

pdf(paste0(args[3],'.pdf'))
levelplot(z[,dim(z)[1]:1], scales=list(x=list(rot=0), tck = c(1,0)), pretty=TRUE, 
          col.regions=colorRampPalette(c("blue", "black","red")), at=seq(-1.0, 1.0, 0.0001), xlab="", ylab="",
          colorkey=list(space="right", col=colorRampPalette(c("blue", "black","red")), 
                        at=seq(-1.0, 1.0, 0.0001), raster=TRUE, tck=2))
garbage <- dev.off()



#levelplot(z[,dim(z)[1]:1], scales=list(x=list(rot=0), tck = c(1,0)), pretty=TRUE, 
#          col.regions=colorRampPalette(c("blue", "black","red")), at=seq(round(min(z)-0.01,digits = 3), 1.0, 0.0001), xlab="", ylab="",
#          colorkey=list(space="right", col=colorRampPalette(c("blue", "black","red")), 
#                        at=seq(round(min(z)-0.01,digits = 3), 1.0, 0.0001), raster=TRUE, tck=2))
