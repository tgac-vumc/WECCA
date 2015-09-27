regioning <-
function (cghdata.called, threshold = 1e-05){

#cghdata.called<-WiltingCalled;threshold<-0.0001

    find.reg.modus <- function(x) {
        splitter <- list()
        splitter[[1]] <- c(1)
        index.temp <- 1
        j <- 1
        for (i in 1:(dim(x)[1] - 1)) {
            if (all(x[i, ] == x[i + 1, ])) {
                index.temp <- c(index.temp, i + 1)
                splitter[[j]] <- index.temp
            }
            else {
                index.temp <- i + 1
                j <- j + 1
                splitter[[j]] <- index.temp
            }
        }
        region.details <- NULL
        for (i in 1:length(splitter)) {
            region.details <- rbind(region.details, c(min(splitter[[i]]), 
                max(splitter[[i]])))
        }
        modus <- which.max(region.details[, 2] - region.details[, 
            1] + 1)
        return(x[region.details[modus[1], 1], ])
    }
    cat("CGHregions of hard call data...")
    cghdata.regions <- CGHregions(cghdata.called, averror = threshold)
    cat("...done", "\n")
    print(paste("threshold used:", threshold, sep = " "))
    calls.annotation <- pData(featureData(cghdata.called))
    regions.annotation <- pData(featureData(cghdata.regions))
    cat("Map regions to clones...")
    reg.to.clones <- list()
    counter <- 0
    for (chr in 1:max(calls.annotation[, 1])) {
        reg.ann.temp <- regions.annotation[regions.annotation[, 
            1] == chr, 1:4]
        for (r in 1:dim(reg.ann.temp)[1]) {
            counter <- counter + 1
            A1 <- which(calls.annotation[, 1] == chr)
            A2 <- which(calls.annotation[, 2] >= reg.ann.temp[r, 
                2])
            A3 <- which(calls.annotation[, 2] <= reg.ann.temp[r, 
                3])
            reg.to.clones[[counter]] <- intersect(intersect(A1, 
                A2), A3)
        }
    }
#    cat("...done", "\n")
#    cghdata.probs <- numeric()
#    for (i in 1:dim(calls(cghdata.called))[2]) {
#        cghdata.probs <- cbind(cghdata.probs, cbind(probloss(cghdata.called)[, 
#            i], probnorm(cghdata.called)[, i], probgain(cghdata.called)[, 
#            i], probamp(cghdata.called)[, i]))
#    }
#    cat("Calculate mode soft call signature for each region...")
#    cghdata.regprobs <- numeric()
#    for (i in 1:length(reg.to.clones)) {
#        cghdata.regprobs <- rbind(cghdata.regprobs, find.reg.modus(cghdata.probs[reg.to.clones[[i]], 
#            , drop = FALSE]))
#    } 
#    rtc <- 1:10;cghdata.call <- calls

    #MODIFIED BY MvdW 25/7/2012, FOR MEMORY EFFICIENCY
    cat("Calculate mode soft call signature for each region...")
    cghdata.regprobs <- numeric()
    for (i in 1:length(reg.to.clones)) {
        rtc <- reg.to.clones[[i]]
        cghdatrtc <- cghdata.called[rtc,,drop=FALSE]
        cghdata.probs <- numeric()
        for (j in 1:dim(calls(cghdata.called))[2]) {
        cghdata.probs <- cbind(cghdata.probs, cbind(probloss(cghdatrtc)[, 
            j], probnorm(cghdatrtc)[, j], probgain(cghdatrtc)[, 
            j], probamp(cghdatrtc)[, j]))
        }
        cghdata.regprobs <- rbind(cghdata.regprobs, find.reg.modus(cghdata.probs))
    }
    #END MODIFICATION

    cat("...done", "\n")
    softcalls.samplenames <- character()
    for (i in 1:dim(calls(cghdata.called))[2]) {
        if (dim(cghdata.regprobs)[2]/dim(calls(cghdata.called))[2] == 
            3) {
            softcalls.samplenames <- c(softcalls.samplenames, 
                paste(c("probloss_", "probnorm_", "probgain_"), 
                  colnames(regions(cghdata.regions))[i], sep = ""))
        }
        if (dim(cghdata.regprobs)[2]/dim(calls(cghdata.called))[2] == 
            4) {
            softcalls.samplenames <- c(softcalls.samplenames, 
                paste(c("probloss_", "probnorm_", "probgain_", 
                  "probamp_"), colnames(regions(cghdata.regions))[i], 
                  sep = ""))
        }
    }
    colnames(cghdata.regprobs) <- softcalls.samplenames
    rownames(cghdata.regprobs) <- rownames(regions(cghdata.regions))
    regdata <- list()
    regdata$ann <- regions.annotation
    regdata$hardcalls <- regions(cghdata.regions)
    regdata$softcalls <- cghdata.regprobs
    return(regdata)
}
