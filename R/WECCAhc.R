WECCAhc <-
function(cghdata.regioned, dist.measure="agree", linkage="ward", weight.type="all.equal"){

hierarchy2dendrogram <- function(hierarchy,order){  
##############################################################################################
# Function that transforms a hierarchy into a dendogram.
##############################################################################################
     tree <- list(merge = hierarchy[,1:2], height= hierarchy[,3], method=NULL, call = match.call(), order = order, dist.method = "whatever")
class(tree) <- "hclust" 
     return(tree) 
} 

merge.agree.avlink.weight <- function(X, a, corr){
#############################################################################################
# Function that constructs the heights and merging pattern,
# with the weighted agreement similarity and average linkage.
##############################################################################################
clusters <- cluster.list(dim(X)[2]-1)
D0 <- agree.comlink.weight.dist(X,clusters,a)
if (corr){ D0 <- matrix(as.numeric(lapply(as.numeric(D0), function(x) max((x-1/a)/(1-1/a),0))),ncol=dim(D0)[1]) }
        D0[upper.tri(D0)] <- t(D0)[upper.tri(t(D0))]
diag(D0) <- 1
hc <- hclust(as.dist(1-D0),method="average")
hc
}

merge.agree.comtradlink.weight <- function(X,a,corr){
#############################################################################################
# Function that constructs the heights and merging pattern, with the weighted agreement 
# similarity and complete linkage in the traditional sense.
##############################################################################################
clusters <- cluster.list(dim(X)[2]-1)
D0 <- agree.comlink.weight.dist(X,clusters,a)
if (corr){ D0 <- matrix(as.numeric(lapply(as.numeric(D0), function(x) max((x-1/a)/(1-1/a),0))),ncol=dim(D0)[1]) }
        D0[upper.tri(D0)] <- t(D0)[upper.tri(t(D0))]
diag(D0) <- 1
hc <- hclust(as.dist(1-D0),method="complete")
return(hc)
}

merge.agree.wardlink.weight <- function(X, a, corr){
#############################################################################################
# Function that constructs the heights and merging pattern, with the weighted agreement 
# similarity and complete linkage in the traditional sense.
##############################################################################################
clusters <- cluster.list(dim(X)[2]-1)
D0 <- agree.comlink.weight.dist(X,clusters,a)
if (corr){ D0 <- matrix(as.numeric(lapply(as.numeric(D0), function(x) max((x-1/a)/(1-1/a),0))),ncol=dim(D0)[1]) }
        D0[upper.tri(D0)] <- t(D0)[upper.tri(t(D0))]
diag(D0) <- 1
hc <- hclust(as.dist(1-D0), method="ward")
return(hc)
}



merge.conc.comtradlink.weight <- function(X,a,corr,strict){
#############################################################################################
# Function that constructs the heights and merging pattern, with the weighted agreement 
# similarity and complete linkage in the traditional sense.
##############################################################################################
clusters <- cluster.list(dim(X)[2]-1)
      D0 <- conc.comlink.weight.dist(X,clusters,a,strict)
      if (!strict){ 
if (corr){ D0 <- matrix(as.numeric(lapply(as.numeric(D0), function(x) max((x-((a^2-2*a+1)/(2*a^2) + 1/a^3))/(1-((a^2-2*a+1)/(2*a^2) + 1/a^3)),0))),ncol=dim(D0)[1]) }
} else {
if (corr){ D0 <- matrix(as.numeric(lapply(as.numeric(D0), function(x) max((x-(a^2-2*a+1)/(2*a^2))/(1-(a^2-2*a+1)/(2*a^2)),0))),ncol=dim(D0)[1]) }
}
D0[upper.tri(D0)] <- t(D0)[upper.tri(t(D0))]
diag(D0) <- 1
hc <- hclust(as.dist(1-D0),method="complete")
return(hc)
}


merge.conc.avlink.weight <- function(X,a,corr,strict){
#############################################################################################
# Function that constructs the heights and merging pattern,
# with the weighted agreement similarity and average linkage.
##############################################################################################
clusters <- cluster.list(dim(X)[2]-1)
        D0 <- conc.comlink.weight.dist(X,clusters,a,strict)
      if (!strict){ 
 if (corr){ D0 <- matrix(as.numeric(lapply(as.numeric(D0), function(x) max((x-((a^2-2*a+1)/(2*a^2) + 1/a^3))/(1-((a^2-2*a+1)/(2*a^2) + 1/a^3)),0))),ncol=dim(D0)[1]) }
} else {
if (corr){ D0 <- matrix(as.numeric(lapply(as.numeric(D0), function(x) max((x-(a^2-2*a+1)/(2*a^2))/(1-(a^2-2*a+1)/(2*a^2)),0))),ncol=dim(D0)[1]) }
}
D0[upper.tri(D0)] <- t(D0)[upper.tri(t(D0))]
diag(D0) <- 1
hc <- hclust(as.dist(1-D0),method="average")
return(hc)
}



merge.conc.wardlink.weight <- function(X,a,corr,strict){
#############################################################################################
# Function that constructs the heights and merging pattern,
# with the weighted agreement similarity and average linkage.
##############################################################################################
clusters <- cluster.list(dim(X)[2]-1)
        D0 <- conc.comlink.weight.dist(X,clusters,a,strict)
      if (!strict){ 
 if (corr){ D0 <- matrix(as.numeric(lapply(as.numeric(D0), function(x) max((x-((a^2-2*a+1)/(2*a^2) + 1/a^3))/(1-((a^2-2*a+1)/(2*a^2) + 1/a^3)),0))),ncol=dim(D0)[1]) }
} else {
if (corr){ D0 <- matrix(as.numeric(lapply(as.numeric(D0), function(x) max((x-(a^2-2*a+1)/(2*a^2))/(1-(a^2-2*a+1)/(2*a^2)),0))),ncol=dim(D0)[1]) }
}
D0[upper.tri(D0)] <- t(D0)[upper.tri(t(D0))]
diag(D0) <- 1
hc <- hclust(as.dist(1-D0),method="ward")
return(hc)
}


merge.agree.comlink.weight <- function(X, a, corr){
#############################################################################################
# Function that constructs the heights and merging pattern,
# with the weighted agreement similarity and complete linkage.
##############################################################################################
clusters <- cluster.list(dim(X)[2]-1)
D0 <- agree.comlink.weight.dist(X,clusters,a)
if (corr){ D0 <- matrix(as.numeric(lapply(as.numeric(D0), function(x) max((x-1/a)/(1-1/a),0))),ncol=dim(D0)[1]) }
height <- NULL
merge <- NULL
for (j in 1:(dim(D0)[1]-1)){
height <- c(height,max(D0[lower.tri(D0)]))
if (height[length(height)] != 0){
   ind <- which(D0 == max(D0[lower.tri(D0)]),arr.ind=TRUE)[1,]
   new.cluster <- sort(c(clusters[[ind[1]]],clusters[[ind[2]]]))
} else {
               ind <- which(lapply(1:length(clusters), function(i, clusters){ is.na(max(match(clusters[[i]], clusters[[1]]))) }, clusters )==TRUE)[1]
   new.cluster <- sort(c(clusters[[1]],clusters[[ind]]))
}
for (i in 1:length(new.cluster)){ clusters[[new.cluster[i]]] <- new.cluster }
when <- rep(0,dim(D0)[1])
when[new.cluster] <- j
merge <- rbind(merge,when)
dist.vec <- NULL
for (i in 1:dim(D0)[1]){
    if (corr){ chance <- 1/a^(length(new.cluster)+length(clusters[i])) } else { chance <- 0 }
          dist.vec <- c(dist.vec, max((prob.agree.comlink.weight(X,1,c(new.cluster),c(clusters[[i]]),a)-chance)/(1-chance),0))
}
dist.vec[new.cluster] <- 0
for (i in 1:length(new.cluster)){ D0[new.cluster[i], ] <- dist.vec }
for (i in 1:length(new.cluster)){ D0[, new.cluster[i]] <- dist.vec }
D0[upper.tri(D0,diag=TRUE)] <- 0
}
return(list(merge=merge,height=height))
}



merge.conc.comlink.weight <- function(X,a,corr,strict){
#############################################################################################
# Function that constructs the heights and merging pattern,
# with the weighted concordance similarity and complete linkage.
##############################################################################################
clusters <- cluster.list(dim(X)[2]-1)
D0 <- conc.comlink.weight.dist(X,clusters,a,strict)
      if (!strict){ 
if (corr){ D0 <- matrix(as.numeric(lapply(as.numeric(D0), function(x) max((x-((a^2-2*a+1)/(2*a^2) + 1/a^3))/(1-((a^2-2*a+1)/(2*a^2) + 1/a^3)),0))),ncol=dim(D0)[1]) }
} else {
if (corr){ D0 <- matrix(as.numeric(lapply(as.numeric(D0), function(x) max((x-(a^2-2*a+1)/(2*a^2))/(1-(a^2-2*a+1)/(2*a^2)),0))),ncol=dim(D0)[1]) }
}
height <- NULL
merge <- NULL
for (j in 1:(dim(D0)[1]-1)){
height <- c(height,max(D0[lower.tri(D0)]))
if (height[length(height)] != 0){
ind <- which(D0 == max(D0[lower.tri(D0)]),arr.ind=TRUE)[1,]
new.cluster <- sort(c(clusters[[ind[1]]],clusters[[ind[2]]]))
} else {
ind <- which(lapply(1:length(clusters), function(i, clusters){ is.na(max(match(clusters[[i]], clusters[[1]]))) }, clusters )==TRUE)[1]
new.cluster <- sort(c(clusters[[1]],clusters[[ind]]))
}
for (i in 1:length(new.cluster)){ clusters[[new.cluster[i]]] <- new.cluster }
when <- rep(0,dim(D0)[1])
when[new.cluster] <- j
merge <- rbind(merge,when)
dist.vec <- NULL
for (i in 1:dim(D0)[1]){
if (!strict){ 
if (corr){ chance <- 2*((a-1)/(2*a))^(length(new.cluster)+length(clusters[i])) + 1/a^(2*(length(new.cluster)+length(clusters[i]))-1) } else { chance <- 0 }
} else {
if (corr){ chance <- 2*((a-1)/(2*a))^(length(new.cluster)+length(clusters[i])) } else { chance <- 0 }
}
dist.vec <- c(dist.vec, max((prob.conc.comlink.weight(X,1,c(new.cluster),c(clusters[[i]]),a,strict)-chance)/(1-chance),0))
}
dist.vec[new.cluster] <- 0
for (i in 1:length(new.cluster)){ D0[new.cluster[i], ] <- dist.vec }
for (i in 1:length(new.cluster)){ D0[, new.cluster[i]] <- dist.vec }
D0[upper.tri(D0,diag=TRUE)] <- 0
}
return(list(merge=merge,height=height))
}



build.dendro <- function(merge,height){
#############################################################################################
# Function that transforms the merging patterns into an efficient one and the permutation 
# order for efficient plotting, and turns this into a dendogram.
##############################################################################################
H1 <- apply(merge,2,cummax)
H <- NULL
for (j in 1:dim(merge)[2]){
H <- c(H,which(H1[,j]==min(H1[(H1[,j]>0),j]),arr.ind=TRUE)[1])
}
M <- (H==1)*c(1:dim(merge)[2])
M <- -1*M[(M!=0)]
for (j in 2:dim(merge)[1]){
id <- (merge[j,]!=0)*c(1:dim(merge)[2])
id <- id[(id!=0)]
H1a <- H1[j-1,id]
for (k1 in 1:length(id)){ 
if(H1a[k1]==0){ H1a[k1] <- -id[k1] }
}
M <- rbind(M,as.numeric(names(table(H1a))))
}
colnames(merge) <- c(1:dim(merge)[2])
for (j in 1:(dim(merge)[1]-1)){
merge <- merge[,order(merge[j,])]
}
order <- as.numeric(colnames(merge))
      mergePheight <- cbind(M,1-height)
colnames(mergePheight) <- NULL
rownames(mergePheight) <- NULL
hc <- hierarchy2dendrogram(mergePheight,order)
return(hc)
}

cluster.list <- function(n){
###################################################################################
# Function that constructs a list for n items equal to 1,...,n.
###################################################################################
c.null <- list()
for (i in 1:n){ c.null[[i]] <- c(i) }
return(c.null)
}


conc.comlink.weight.dist <- function(X,clusters,a,strict){
##############################################################################################
# Function that constructs the concordance similarity matrix with complete linkage.
##############################################################################################
dist.mat <- matrix(0,nrow=length(clusters),ncol=length(clusters))
for (i1 in 1:(length(clusters)-1)){
for (i2 in (i1+1):length(clusters)){
dist.mat[i2,i1] <- prob.conc.comlink.weight(X,1,as.numeric(clusters[i1]),as.numeric(clusters[i2]),a,strict)
}
}
return(dist.mat)
}



agree.comlink.weight.dist <- function(X,clusters,a){
##############################################################################################
# Function that constructs the agreement similarity matrix with complete linkage.
##############################################################################################
dist.mat <- matrix(0,nrow=length(clusters),ncol=length(clusters))
for (i1 in 1:(length(clusters)-1)){
for (i2 in (i1+1):length(clusters)){
dist.mat[i2,i1] <- prob.agree.comlink.weight(X,1,as.numeric(clusters[i1]),as.numeric(clusters[i2]),a)
}
}
return(dist.mat)
}



prob.agree.comlink.weight <- function(X,w.col,c1.cols,c2.cols,a){
##############################################################################################
# Function that calculates the weighted probability of agreement with complete linkage.
##############################################################################################
X1 <- matrix(X[,c1.cols+1],ncol=length(c1.cols),nrow=dim(X)[1])
X2 <- matrix(X[,c2.cols+1],ncol=length(c2.cols),nrow=dim(X)[1])
Xt <- cbind(X[,w.col],X1,X2)
index <- as.numeric(factor(apply(Xt[,2:dim(Xt)[2]], 1, paste, collapse=":"))) 
Xt <- cbind(Xt,index)
Xt <- Xt[order(Xt[,dim(Xt)[2]]),]
counter <- NULL
for (i in 1:length(table(index))){
counter <- c(counter, cumsum(Xt[(Xt[,dim(Xt)[2]]==i),dim(Xt)[2]])/i)
}
Xt <- cbind(Xt,counter,1)
Xc <- NULL
for (i in 1:length(table(index))){
Xc <- rbind(Xc,c(i,sum(Xt[(Xt[,dim(Xt)[2]-2]==i),dim(Xt)[2]]),Xt[((Xt[,dim(Xt)[2]-2]==i)&(Xt[,dim(Xt)[2]-1]==1)),2:(dim(Xt)[2]-3)],sum(Xt[(Xt[,dim(Xt)[2]-2]==i),1]),sum((Xt[(Xt[,dim(Xt)[2]-2]==i),1])^2)))
}
Pa <- 0
for (i in 1:a){
K <- rep(i,length(c1.cols)+length(c2.cols))
P <- matrix((Xc[,3:(length(c1.cols)+length(c2.cols)+2)]==K),ncol=(length(c1.cols)+length(c2.cols)+2)-(3-1))
H <- apply(P,1,prod)*c(1:dim(P)[1])
if (sum(H)>0){ Pa <- Pa + Xc[H,(length(c1.cols)+length(c2.cols)+3)] }
}
Pa <- Pa/sum(X[,w.col])
return(Pa)
}


prob.conc.comlink.weight <- function(X, w.col, c1.cols, c2.cols, a, strict){
##############################################################################################
# Function that calculates the weighted probability of concordance with complete linkage.
##############################################################################################

# select relevant data
X <- X[, c(w.col, c1.cols+1, c2.cols+1)]

# add column with indicator for each unique pattern
# X <- X[order(apply(X[,-1], 1, sum)), ]
# X <- X[,c(1, order(apply(X[,-1], 2, sum)) + 1)]
X <- cbind(X, as.numeric(factor(apply(X[,2:dim(X)[2]], 1, paste, collapse=":"))) )

# sort by pattern indicator
X <- X[order(X[,dim(X)[2]]), ]

# create intermediate variables for counting the response pattern frequency later
counter <- NULL
for (i in 1:length(table(X[,dim(X)[2]]))){
counter <- c(counter, cumsum(X[(X[,dim(X)[2]]==i), dim(X)[2]])/i)
}
X <- cbind(X, counter, 1)

# make table, denoted Xc
# 1st column: response pattern indicator
# 2nd column: response pattern frequency
# next columns: response patterns
# last two columns summary statistic of the weights
Xc <- NULL
for (i in 1:length(table(X[,dim(X)[2]-2]))){
Xc <- rbind(Xc, c(i, sum(X[(X[,dim(X)[2]-2]==i), dim(X)[2]]), 
X[((X[,dim(X)[2]-2]==i) & (X[, dim(X)[2]-1]==1)), 2:(dim(X)[2]-3)], 
sum(X[(X[,dim(X)[2]-2]==i), 1]), sum((X[(X[,dim(X)[2]-2]==i), 1])^2)))
}

m <- Xc[,3:(dim(Xc)[2]-2), drop=FALSE]

# find regions with tied aberrations pattern
equal <- function(X, P){ all(X==P) }
welke <- NULL
if (is.matrix(m)){
for (i in 1:a){ welke <- cbind(welke, apply(m, 1, equal, i)) }
}
if (is.matrix(m)==FALSE){
for (i in 1:a){ welke <- cbind(welke, all(m==i)) }
}

# calculate tie contribution of "tied regions" to concordance probability
if (is.matrix(welke)){
Pc.e <- sum(Xc[apply(welke, 1, any), dim(Xc)[2]-1]^2) - sum(Xc[apply(welke, 1, any), dim(Xc)[2]])
}
if (is.numeric(welke)){
Pc.e <- sum(Xc[apply(welke, 1, any), dim(Xc)[2]-1]^2) - sum(Xc[apply(welke, 1, any), dim(Xc)[2]])
}

# calculate 
Pc.l <- 0
bigger <- function(X, P){ all(X > P) }
for (i in 1:(nrow(Xc))){
Ui <- apply(Xc[, 3:(ncol(Xc)-2)], 1, "bigger", Xc[i, 3:(ncol(Xc)-2)]) * c(1:nrow(Xc))
Ui <- Ui[(Ui!=0)]    
if (length(Ui) != 0){ Pc.l <- Pc.l + Xc[i, (ncol(Xc)-1)] * sum(Xc[Ui, (ncol(Xc)-1)]) }
}

# sort out which definition of the concordance probability to use
if (!strict){ Pc <- (2*Pc.l + Pc.e) / (sum(X[, w.col])^2 - sum(X[, w.col]^2)) }
if (strict){ Pc <- (2*Pc.l) / (sum(X[, w.col])^2 - sum(X[, w.col]^2)) }

return(Pc)
}

##############################################################################################
# Actual start of function WECCAhc.
##############################################################################################

# make object for compatability
corr <- FALSE
strict <- FALSE
a <- dim(cghdata.regioned$softcalls)[2] / dim(cghdata.regioned$hardcalls)[2]

# determine weights
cat("Weight calculation...")
if (weight.type == "all.equal"){
weights <- rep(1, dim(cghdata.regioned$hardcalls)[1])
}

# the orginal wecca works with different coding of aberration. transform if coding differs.
if (min(cghdata.regioned$hardcalls) == -1){
cghdata.regioned$hardcalls <- cghdata.regioned$hardcalls + 2
}

if (weight.type == "heterogeneity"){
minCall <- min(cghdata.regioned$hardcalls)
maxCall <- max(cghdata.regioned$hardcalls)
callLevels <- seq(from=minCall, to=maxCall, by=1)
sumData <- t(apply(cghdata.regioned$hardcalls, 1, function(Z, callLevels){ as.numeric(table(factor(Z, levels=callLevels))) }, callLevels=callLevels)) / dim(cghdata.regioned$hardcalls)[2]
weights.entropy <-  apply(sumData, 1, function(x){ x*log(x, 2) })
weights.entropy[is.na(weights.entropy)] <- 0
weights.entropy <-  -colSums(weights.entropy)
weights <- dim(cghdata.regioned$hardcalls)[1] * weights.entropy / sum(weights.entropy)
}
cat("...done.", "\n")
data <- cbind(weights, cghdata.regioned$hardcalls)

cat("Distance matrix calculation and dendrogram construction...", "\n")
if (dist.measure=="conc" && linkage=="total" && corr==FALSE && strict==FALSE){
cat("NOTE: total-linkage may take a while...", "\n")
hier <- merge.conc.comlink.weight(data,a,corr=FALSE,strict=FALSE)
dendro <- build.dendro(hier$merge,hier$height)
}
if (dist.measure=="conc" && linkage=="total" && corr==FALSE && strict==TRUE){
cat("NOTE: total-linkage may take a while...", "\n")
hier <- merge.conc.comlink.weight(data,a,corr=FALSE,strict=TRUE)
dendro <- build.dendro(hier$merge,hier$height)
}
if (dist.measure=="conc" && linkage=="total" && corr==TRUE && strict==FALSE){
cat("NOTE: total-linkage may take a while...", "\n")
hier <- merge.conc.comlink.weight(data,a,corr=TRUE,strict=FALSE)
dendro <- build.dendro(hier$merge,hier$height)
}
if (dist.measure=="conc" && linkage=="total" && corr==TRUE && strict==TRUE){
cat("NOTE: total-linkage may take a while...", "\n")
hier <- merge.conc.comlink.weight(data,a,corr=TRUE,strict=TRUE)
dendro <- build.dendro(hier$merge,hier$height)
}
if (dist.measure=="conc" && linkage=="complete" && corr==FALSE && strict==TRUE){
dendro <- merge.conc.comtradlink.weight(data,a,corr=FALSE,strict=TRUE)
}
if (dist.measure=="conc" && linkage=="complete" && corr==FALSE && strict==FALSE){
dendro <- merge.conc.comtradlink.weight(data,a,corr=FALSE,strict=FALSE)
}
if (dist.measure=="conc" && linkage=="complete" && corr==TRUE && strict==TRUE){
dendro <- merge.conc.comtradlink.weight(data,a,corr=TRUE,strict=TRUE)
}
if (dist.measure=="conc" && linkage=="complete" && corr==TRUE && strict==FALSE){
dendro <- merge.conc.comtradlink.weight(data,a,corr=TRUE,strict=FALSE)
}
if (dist.measure=="conc" && linkage=="average" && corr==FALSE && strict==FALSE){
dendro <- merge.conc.avlink.weight(data,a,corr=FALSE,strict=FALSE)
}
if (dist.measure=="conc" && linkage=="average" && corr==FALSE && strict==TRUE){
dendro <- merge.conc.avlink.weight(data,a,corr=FALSE,strict=TRUE)
}
if (dist.measure=="conc" && linkage=="average" && corr==TRUE && strict==FALSE){
dendro <- merge.conc.avlink.weight(data,a,corr=TRUE,strict=FALSE)
}
if (dist.measure=="conc" && linkage=="average" && corr==TRUE && strict==TRUE){
dendro <- merge.conc.avlink.weight(data,a,corr=TRUE,strict=TRUE)
}
if (dist.measure=="conc" && linkage=="ward" && corr==FALSE && strict==FALSE){
dendro <- merge.conc.wardlink.weight(data,a,corr=FALSE,strict=FALSE)
}
if (dist.measure=="conc" && linkage=="ward" && corr==FALSE && strict==TRUE){
dendro <- merge.conc.wardlink.weight(data,a,corr=FALSE,strict=TRUE)
}
if (dist.measure=="conc" && linkage=="ward" && corr==TRUE && strict==FALSE){
dendro <- merge.conc.wardlink.weight(data,a,corr=TRUE,strict=FALSE)
}
if (dist.measure=="conc" && linkage=="ward" && corr==TRUE && strict==TRUE){
dendro <- merge.conc.wardlink.weight(data,a,corr=TRUE,strict=TRUE)
}
if (dist.measure=="agree" && linkage=="total" && corr==FALSE){
cat("NOTE: total-linkage may take a while...", "\n")
hier <- merge.agree.comlink.weight(data,a,corr=FALSE)
dendro <- build.dendro(hier$merge,hier$height)
}
if (dist.measure=="agree" && linkage=="total" && corr==TRUE){
cat("NOTE: total-linkage may take a while...", "\n")
hier <- merge.agree.comlink.weight(data,a,corr=TRUE)
dendro <- build.dendro(hier$merge,hier$height)
}
if (dist.measure=="agree" && linkage=="complete" && corr==FALSE){
dendro <- merge.agree.comtradlink.weight(data,a,corr=FALSE)
}
if (dist.measure=="agree" && linkage=="complete" && corr==TRUE){
dendro <- merge.agree.comtradlink.weight(data,a,corr=TRUE)
}
if (dist.measure=="agree" && linkage=="average" && corr==FALSE){
dendro <- merge.agree.avlink.weight(data,a,corr=FALSE)
}
if (dist.measure=="agree" && linkage=="average" && corr==TRUE){
dendro <- merge.agree.avlink.weight(data,a,corr=TRUE)
}
if (dist.measure=="agree" && linkage=="ward" && corr==FALSE){
dendro <- merge.agree.wardlink.weight(data, a, corr=FALSE)
}
if (dist.measure=="agree" && linkage=="ward" && corr==TRUE){
dendro <- merge.agree.wardlink.weight(data, a, corr=TRUE)
}
cat("...done.", "\n")
return(dendro)
}
