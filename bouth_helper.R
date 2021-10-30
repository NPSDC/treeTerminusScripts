extractBouthNodes <- function(bouth_ob, se, type = "driver") {
    if(type=="driver")
        detNodes <- bouth_ob$results.by.node[bouth_ob$results.by.node$is.driver,]
    else
        detNodes <- bouth_ob$results.by.node[bouth_ob$results.by.node$is.detected,]
    levels <- ncol(detNodes) - 3
    det <- c()
    for(i in seq(nrow(detNodes)))
    {
        for(j in seq(levels))
        {
            if(detNodes[i,j]!=" " & detNodes[i,j+1]==" ")
            {
                node <- detNodes[i,j]
                if(node=="Root")
                    det <- nrow(y)+1
                else
                    det <- c(det, as.numeric(detNodes[i,j]))
                break()
            }

            else if(j==levels)
                det <- c(det, match(detNodes[i,j], rownames(y)))
        }    
    }
    return(det)
}