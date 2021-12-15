getRedInfRV <- function(tree, ytemp, mapDf=NULL, cores = 1) {
    if(length(tree$tip) != nrow(ytemp)) {
        if(is.null(mapDf))
            stop("MapDf cant be NULL")
        tinds <- as.numeric(tree$tip)+1 ##indexes in tree corresponding to equivalence class file, 1 for R based
        tnames <- rownames(mapDf)[tinds] ##transcripts in y
        ytemp <- y[tnames,]
    }
        
    ytemp <- prepSwish(tree, ytemp)
    ytemp <- computeInfRV(ytemp)
    
    innNodes <- length(tree$tip)+1:tree$Nnode
    netRedIRV <- 0
    childNodes <- Descendants(tree, innNodes, "child")
    if(class(childNodes)!="list")
        childNodes <- list(childNodes)
    
    netRedIRV <- sum(unlist(mclapply(seq_along(innNodes), function(n) mcols(ytemp)[innNodes[n], "meanInfRV"] - mean(mcols(ytemp)[childNodes[[n]], "meanInfRV"]), mc.cores = cores)))
    
    return(netRedIRV)
}

findBestTree <- function(tree, y, mapDf, type = "infRV", cores=1) {
    tinds <- as.numeric(tree$tip)+1 ##indexes in tree corresponding to equivalence class file, 1 for R based
    tnames <- rownames(mapDf)[tinds] ##transcripts in y
    ytemp <- y[tnames,]
    
    if(length(tree$tip) == 2) {
        redInfRV <- getRedInfRV(tree, ytemp, cores=cores)
        return(list(tree = write.tree(tree), "infRVRed" = redInfRV))
    }
    #tree <- unroot(tree)
    ts <- mclapply(length(tree$tip)+1:tree$Nnode, function(n) root(tree, node=n), mc.cores=2)
    
    if(type=="infRV") {
        redInfRV <- unlist(mclapply(ts, function(t) getRedInfRV(t, ytemp, cores=2), mc.cores=cores))
        return(list(tree = gsub(';',':0;',write.tree(ts[[which.min(redInfRV)]])), "infRVRed" = redInfRV))
    }
}

computeRF <- function(t1, t2, rooted = T, norm = F) {
    comm <- intersect(t1$tip, t2$tip)
    t1u <- drop.tip(t1, setdiff(t1$tip,comm))
    t2u <- drop.tip(t2, setdiff(t2$tip,comm))
    RF.dist(t1u, t2u, rooted = rooted, normalize = norm)
}

readInpTrees <- function(filepath) {
    con = file(filepath, "r")
    i=1
    trees <- list()
    group <- ""
    while (TRUE) {
        line = readLines(con, n = 1)
        if(length(line) == 0)
            break()
        if(grepl(">Group ", line)) {
            group <- strsplit(line, " ", fixed=T)[[1]][2]
            trees[[group]] <- list()
            j=1
        }
        else {
            tr <- read.tree(text=gsub(';',':0;', line))
                trees[[group]][[j]] <- tr
                j=j+1
                
        }
        #i = i+1
    }
    close(con)
    return(trees)
}
