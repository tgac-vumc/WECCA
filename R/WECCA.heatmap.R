WECCA.heatmap <-
function (cghdata.regioned, dendrogram,...){
######################################################################################################
# heatmap with hard calls and chrosome indicator 
######################################################################################################

# calculated number of classes used
nclass <- length(unique(as.numeric(cghdata.regioned$hardcalls))) #28/7/2011: adjustment MvdW

# Preparation of the plotting of the heatmap
# Generate alternating colors for chromosomes.
chr.color <- rep("blue", dim(cghdata.regioned$hardcalls)[1])
ids <- ((cghdata.regioned$ann[, 1]%%2) == 0)
chr.color[ids] <- c("yellow")

# Generate labels for begin points of chromosomes.    
Y <- rep(FALSE, dim(cghdata.regioned$hardcalls)[1])
    for (i in 2:(dim(cghdata.regioned$ann)[1])) {
        if ((cghdata.regioned$ann[i - 1, 1] != cghdata.regioned$ann[i, 1])) {
            Y[i] <- TRUE
        }
}
Y[1] <- TRUE
begin.chr <- rep("", dim(cghdata.regioned$ann)[1])
begin.chr[Y] <- cghdata.regioned$ann[Y, 1]

# The heatmap is plotted.    
color.coding <- c("red", "black", "green", "white")[1:nclass]
heatmap(cghdata.regioned$hardcalls, Colv = as.dendrogram(dendrogram), 
        Rowv = NA, col = color.coding, labRow = begin.chr,RowSideColors = chr.color, 
        scale = "none",...) #28/7/2011: adjustment MvdW, added ... to allow for using heatmap options in the mother function
}
