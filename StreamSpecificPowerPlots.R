setwd("V:/Analysis/5_Coastwide/Multispecies/Alaska Hatchery Research Program/Parentage Simulations/Simulations based on Christie et al. 2015 Box 2")

source("H:/R Source Scripts/Functions.GCL_KS.R")

#### Target Dir
# setwd("StreamSpecific")  # var.n = var.h
setwd("StreamSpecific_ConstSize")  # size.n = size.h
# setwd("StreamSpecific_ConstSize_RRS_0.8")  # size.n = size.h

stream.dat <- read.table(file = "V:/Presentations/Public/AHRP/March 2016 Science Panel Meeting/KyleStreamData.txt", header = TRUE, stringsAsFactors = FALSE)
stream.mat <- data.matrix(stream.dat[, 2:5])
rownames(stream.mat) <- stream.dat$Stream


#stream <- c("Spring", "Stockdale", "Gilmour", "HoganBay", "Erb", "Paddy", "Fish", "Admiralty", "Prospect", "Sawmill")[1]
stream <- c("Erb", "HoganBay", "Paddy", "Spring", "Stockdale", "Gilmour", "Admiralty", "Fish", "Prospect", "Sawmill")[6]
year <- 2013
props <- c(0.05, 0.1, 0.167, 0.333, 0.5, 0.667, 0.833, 1)
kSize <- 10
# for(kSize in c(1, 2, 5, 10)){
for(stream in c("Erb", "HoganBay", "Paddy", "Spring", "Stockdale")) {
  
  OUT <- setNames(object = sapply(props, function(prop) {read.table(file = grep(pattern = paste(year, "_prop_", prop, "_", sep=""), x = list.files(path=getwd(), pattern=stream, full.names=TRUE), value = TRUE))}, simplify = FALSE), nm = props)
  
  ## Isolate a specific "size.n"
  # different size.n values
  unique(OUT[[1]][, "size.n"])  # if NA shows up, there is a problem
  
  # Check for NAs
  lapply(OUT, function (prop) {prop[which(is.na(prop[, "size.n"])), ]})
  #OUT[["0.05"]][4, 1:3] <- c(0.5, 0.25, 10)  # Fix NAs
  
  newlst <- lapply(OUT, function(prop) {prop[which(prop[,"size.n"] == kSize), ]})
  
  newmat <- sapply(newlst, function(prop) {prop[, "Perm.Power"]})
  rownames(newmat) <- unique(newlst[[1]][,"mu.n"])
  
  par(mar = c(5.1, 6.1, 5.1, 6.1))
  par(mfrow = c(1, 1))
  par(bg = "white")
  #   par(bg = colors()[c(356)])
  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  zmin <- 0
  zmax <- 1
  
  filled.contour3(x = as.numeric(rownames(newmat)),y = as.numeric(colnames(newmat)), z = newmat, color = jet.colors,
                  xlim= c (0, 5), ylim = c(0, max(as.numeric(colnames(newmat)))), zlim = c(zmin,zmax),
                  xlab = "", ylab = "", cex.lab = 2, axes = FALSE)
  
  mtext(expression(paste("F"[1], " sampling proportion", sep = "")), side = 2, line = 3.5, cex = 2, las = 0)
  mtext(side = 1, "Mean RS Natural-origin", cex = 2, line = 3)
  mtext(side=3, paste(stream, " ", year, sep = ""), cex = 2, line = 3)
  mtext(side=3, paste("Dispersion = ", kSize, "; N = ", OUT[[1]][1, "N.w.parents"], "; H = ", OUT[[1]][1, "N.h.parents"], sep = ""), cex = 2, line = 1)
  plot.axes = { axis(1, at=seq(0,5,by=1), cex.axis=1.5);
                axis(2, at=c(0,as.numeric(colnames(newmat))),labels=c("0","1/20","1/10","1/6","1/3","1/2","2/3","5/6","1"),las=1,cex.axis=1.5) }
  
  
  abline(v = stream.mat[stream, c(2, 4)], h = stream.mat[stream, c(1, 3)], lwd = 3, lty = 2)
  polygon(x = rep(range(stream.mat[stream, c(2, 4)]), each = 2), y =  c(range(stream.mat[stream, c(1, 3)]), rev(range(stream.mat[stream, c(1, 3)]))), col = rgb(1, 1, 1, 0.5))
  # abline(v = c(0.82, 1.87), h = c(0.45, 0.9), lwd = 3, lty = 2)
  # polygon(x = c(0.82, 0.82, 1.87, 1.87), y = c(0.45, 0.9, 0.9, 0.45), col = rgb(1, 1, 1, 0.5))  # , lwd = 5, lty = 2
  
  #   points(x = rep(x = as.numeric(rownames(newmat)), each = ncol(newmat)), y = rep(x = as.numeric(colnames(newmat)), times = nrow(newmat)), cex = 2, pch = 16)  # add points to show how interpolations are done
  contour(x = as.numeric(rownames(newmat)), y = as.numeric(colnames(newmat)), z = newmat, levels = 0.8, add = TRUE, lwd = 5, labcex = 2.5, method = "edge")
  #   
  par(new = "TRUE",plt = c(0.87,0.92,0.20,0.80),las = 1,cex.axis = 1.5)
  filled.legend(as.numeric(rownames(newmat)),y=as.numeric(colnames(newmat)),newmat,color=jet.colors,xlab="",ylab="",xlim=c(min(xintercepts),max(xintercepts)),ylim=c(min(slopes),max(slopes)),zlim=c(zmin,zmax))
  filled.legend(as.numeric(rownames(newmat)),y=as.numeric(colnames(newmat)),newmat,color=jet.colors,xlab="",ylab="",xlim=c(min(xintercepts),max(xintercepts)),ylim=c(min(slopes),max(slopes)),zlim=c(zmin,zmax))
  mtext("Power",cex=2,padj=-0.3, adj = 0.2)
}