impScores.plot <-
function(impScores.table, yaxis.max=10, lastchr=22, color="red"){

######################################################################################################
#
# function that plots importance scores
#
######################################################################################################

# internal function needed for plotting
rangebpchr <- function(chr, data){
datchr <- data[data[,1]==chr,]
en <- nrow(datchr)
rangebp <- datchr[en,2]
return(rangebp/en)
}

# blow up to total numer of features on the array
impScores.table.extended <- numeric()
for (r in 1:dim(impScores.table)[1]){
junk <- matrix(rep(impScores.table[r,-c(4:5)], impScores.table[r,5]), ncol=4, byrow=TRUE)
step.size <- as.numeric(round((impScores.table[r,4] - impScores.table[r,3]) / (impScores.table[r,5]-1), digits=0))
grid <- impScores.table[r,3] + step.size * c(0:(impScores.table[r,5]-1))
junk[,3] <- grid
impScores.table.extended <- rbind(impScores.table.extended, junk)
}

# extract relevant info from the extended table
bppos <- impScores.table.extended[,3]
chrom <- impScores.table.extended[,2]
chrbp <- cbind(chrom, bppos) 
impScores <- impScores.table.extended[,4]

# dit is de schaling. scaleall is een vector ter lengte aantal datapunten; de waardes zijn constant per 
# chromosoom. In de barplot functie zie je dat de wijdte van de bars (1 bar per clone) wordt 
# opgeschaald met scaleall. Zo zorg ik ervoor dat de lengte van een chromosoom klopt. De precieze
# x-as positie klopt niet met de bp-positie maar dat is ook niet zo belangrijk
rangeallchr <- sapply(unique(chrom),rangebpchr,data=chrbp)
sca <- rangeallchr/max(rangeallchr)
scaleall <- sapply(chrom,function(x)sca[x])
chrom.labels <- as.character(unique(chrom))

# set plot border
par(mar = c(5, 4, 4, 4) + 0.2)
# plot
G <- barplot(rep(0, length(impScores)), width=scaleall, border="black",space=0,col=FALSE, las=1, cex.axis=1, cex.lab=1, ylim=c(0, min(max(impScores)+1, yaxis.max)), beside=TRUE, ylab="importance score", xlab="chromosome", main="Importance score plot")
barplot(rep(0, length(impScores)), width=scaleall, border="black",space=0,col=FALSE, las=1, cex.axis=1, cex.lab=1, ylim=c(0, min(max(impScores)+1, yaxis.max)), beside=TRUE, ylab="importance score", xlab="chromosome", main="Importance score plot")
lines(y=t(impScores),x=G, type="l", col="red", lwd=1.5)
lim <- par("usr")
lim[3:4] <- c(-5,5)
par(usr=lim)
# add vert lines at chromosome ends
chromnew <- table(chrom) * sca
for (iii in 1:length(cumsum(chromnew))) {
segments(cumsum(chromnew)[[iii]],-5,cumsum(chromnew)[[iii]],5,lty=2)
}
box()
# add chromomsome labels
ax <- (cumsum(chromnew) + c(0,cumsum(chromnew)[-length(cumsum(chromnew))])) / 2
axis(side=1, at=ax, labels=chrom.labels, cex=.2, lwd=.25, las=1, cex.axis=1, cex.lab=1)

}
