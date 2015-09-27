WECCAsc <-
function(cghdata.regioned, dist.measure="KLdiv", linkage="ward", weight.type="all.equal"){

######################################################################################################
#
# WECCA for call probabilities
#
######################################################################################################

reg.dist.ordinal <- function(one.reg.probs.two.samples, nclass){

##############################################################################################
# 
# Function that calculates the difference between the cumulative call
# probability distributions of two regions.
#
##############################################################################################

return(sum(abs(cumsum(one.reg.probs.two.samples[1:nclass]) - cumsum(one.reg.probs.two.samples[-c(1:nclass)]))))
}

reg.dist.KLdiv <- function(one.reg.probs.two.samples, nclass, precision){

##############################################################################################
# 
# Function that calculates the difference between the cumulative call
# probability distributions of two regions.
#
##############################################################################################

# return(one.reg.probs.two.samples[1:nclass] * log(one.reg.probs.two.samples[-c(1:nclass)]) + one.reg.probs.two.samples[-c(1:nclass)] * log(one.reg.probs.two.samples[1:nclass]))
KL.matrix <- KLdiv(cbind(one.reg.probs.two.samples[1:nclass], one.reg.probs.two.samples[-c(1:nclass)]), overlap=FALSE, eps=max((min(one.reg.probs.two.samples))^2, precision))
return((KL.matrix[1,2] + KL.matrix[2,1])/2)
}

##############################################################################################
# 
# Function that constructs the agreement similarity matrix with complete linkage.
#
##############################################################################################

weight.dist.matrix.ordinal <- function(cgh.reg.probs, nclass, weights){

##############################################################################################
# 
# Function that constructs the distance matrix with the distance measure 
# that takes the ordinality of the calls into account.
# It calculates the distance between the cumulative distributions for each region.
#
##############################################################################################

dist.mat <- matrix(0, nrow=dim(cgh.reg.probs)[2]/nclass, ncol=dim(cgh.reg.probs)[2]/nclass)
for (i1 in 1:(dim(cgh.reg.probs)[2] / nclass - 1)){
for (i2 in (i1+1):(dim(cgh.reg.probs)[2] / nclass)){
corresponding.cols <- c(nclass * (i1-1) + c(1:nclass), nclass*(i2-1) + c(1:nclass))
dist.mat[i2,i1] <- sum(weights * apply(cgh.reg.probs[,corresponding.cols], 1, reg.dist.ordinal, nclass=nclass))
dist.mat[i1,i2] <- dist.mat[i2,i1]
}
}
return(dist.mat)
}

weight.dist.matrix.KLdivergence <- function(cgh.reg.probs, nclass, weights){

##############################################################################################
# 
# Function that constructs the distance matrix with the distance measure 
# that does not take the ordinality of the calls into account.
# It calculates the Kullback-Leibler divergence between the distributions for each region.
#
##############################################################################################

precision <- sort(unique(as.numeric(cgh.reg.probs)))[2:201]
precision <- min(precision[(precision!=0)], na.rm=TRUE)/10

dist.mat <- matrix(0, nrow=dim(cgh.reg.probs)[2]/nclass, ncol=dim(cgh.reg.probs)[2]/nclass)
for (i1 in 1:(dim(cgh.reg.probs)[2]/nclass - 1)){
for (i2 in (i1+1):(dim(cgh.reg.probs)[2] / nclass)){
corresponding.cols <- c(nclass * (i1-1) + c(1:nclass), nclass * (i2-1) + c(1:nclass))
dist.mat[i2,i1] <- sum(weights * apply(cgh.reg.probs[,corresponding.cols], 1, reg.dist.KLdiv, nclass=nclass, precision=precision))
dist.mat[i1,i2] <- dist.mat[i2,i1]
}
}
return(dist.mat)
}

# START OF WECCAsc FUNCTION

# extract cgh data from data object
cgh.reg.probs <- cghdata.regioned$softcalls

# calculated number of classes used
nclass <- dim(cghdata.regioned$softcalls)[2] / dim(cghdata.regioned$hardcalls)[2]

# calculate number of samples present
nosamp <- dim(cghdata.regioned$softcalls)[2]

# determine weights
cat("Weight calculation...")
if (weight.type == "all.equal"){
weights <- rep(1, dim(cgh.reg.probs)[1])
}
if (weight.type == "heterogeneity"){
total.call.probs <- numeric()
for (k in 1:nclass){
total.call.probs <- cbind(total.call.probs, rowSums(cgh.reg.probs[,nclass*(c(1:(dim(cgh.reg.probs)[2]/nclass)) - 1) + k])/(dim(cgh.reg.probs)[2]/nclass))
}
weights.entropy <-  apply(total.call.probs, 1, function(x){ x*log(x, 2) })
weights.entropy[is.na(weights.entropy)] <- 0
weights.entropy <-  -colSums(weights.entropy)
weights <- dim(cgh.reg.probs)[1] * weights.entropy / sum(weights.entropy)
}
cat("...done.", "\n")

# calculate the distance matrix
cat("Distance matrix calculation...")
if (dist.measure == "KLdiv"){
dist.mat <- weight.dist.matrix.KLdivergence(cgh.reg.probs, nclass, weights)
} else {
dist.mat <- weight.dist.matrix.ordinal(cgh.reg.probs, nclass, weights)
}
cat("...done.", "\n")

# construct dendrogram
cat("Dendrogram construction...")
dendro <- hclust(as.dist(dist.mat), method=linkage)
cat("...done.", "\n")

# return dendro
return(dendro)
}
