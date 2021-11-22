.bottom.up.weighted<-function(tree, pvalue, far = 0.10, tau = 0.3, mode = 'one.stage'){
    
    get_weight<-function(hash.id, tree){## calculate the least favorable weights
        
        n.node = nrow(tree@node)
        weight = array(1, dim = n.node)
        hash.xid = hash.id
        subtree<-vector("list",n.node)
        this = hash.id[1]
        
        while(this <= n.node){
            parent = tree@node$parent[this]
            if(parent <= n.node & !parent %in% hash.id){
                hash.id =  c(hash.id, parent)
            }
            hash.id = hash.id[-1]
            if(length(hash.id)==0){
                break
            }
            if(parent <= n.node){
                subtree[[parent]] = c(subtree[[parent]], this)
            }
            this = hash.id[1]
        }
        subtree = lapply(subtree, unique)
        
        ## compute least favorable weights
        LFW<-vector("list",n.node); LFW[hash.xid] = 1
        for(i in 1:n.node){
            if(length(subtree[[i]]) > 0){
                LFW[[i]] = sort(unlist(LFW[subtree[[i]]]), decreasing = TRUE)
                LFW[[i]][1] = LFW[[i]][1] + 1
            }
        }
        return(sort(LFW[[n.node]]))
    }
    
    
    n.node = nrow(tree@node)
    max.depth = max(tree@node$level)
    if(length(pvalue)!=sum(tree@node$level==1)){
        stop("Incorrect dimension of p-values!")
    }
    if(!is.null(names(pvalue))){
        level.one.id = tree@node$id[tree@node$level==1]
        if(setequal(level.one.id, names(pvalue))){
            pvalue = pvalue[level.one.id]
        }else{
            stop("Mis-matched symbols between p-values and tree")
        }
    }
    q.spent = array(0, max.depth)
    alpha.out = array(NA, max.depth)
    nl.out = array(NA, max.depth)
    
    if(mode == 'one.stage'){
        for(i in 1:max.depth){
            nl.out[i] = sum(tree@node$level == i)
            q.spent[i] = far/n.node*nl.out[i]
        }
    }
    if(mode == 'stageI.in.two.stage'){
        q.spent[1] = far
        nl.out = sapply(c(1:max.depth),function(x)return(sum(tree@node$level==x)))
    }
    if(mode == 'stageII.in.two.stage'){
        if(max.depth<2){stop('Two stage analysis requires at least 2 levels.')}
        n.node1 = n.node - length(which(tree@node$level == 1))
        for(i in 2:max.depth){
            q.spent[i] = far/n.node1*length(which(tree@node$level == i))
        }
    }
    
    cum.rej = 0; last.cut.off = 0
    sig.out = array(FALSE,n.node)
    pvalue.in = pvalue.out = dec.level = array(NA, n.node); pvalue.in[1:length(pvalue)]= pvalue
    
    
    for(i in 1:max.depth){
        this.level.set = which(tree@node$level == i)
        q.this.level = q.spent[i] ### by default q-spending by proportion
        this.level.tree = tree@edge[this.level.set]
        
        # removing significant children from previous levels
        if(i>1){
            
            ## nodes that are significant before testing this level
            pre.sig.nodes = which(sig.out[1:(min(this.level.set)-1)])
            ## nodes to be tested now
            this.level.tree = lapply(this.level.tree,function(x)return(setdiff(x,pre.sig.nodes)))
            ind = (sapply(this.level.tree,length) > 0)
            above.level.set = which(tree@node$level >= i)
            xind = (sapply(lapply(tree@edge[above.level.set],function(x)return(setdiff(x,pre.sig.nodes))),length)>0)
            xind = (xind | sig.out[above.level.set])
            na.nodes = above.level.set[!xind]
            sig.out[na.nodes] = TRUE
            pvalue.out[na.nodes] = NA
            
            
            if(length(na.nodes)>0){
                # find its least significant offspring on the i-1 th level
                last.level.set = intersect(which(tree@node$level == (i-1)),which(sig.out))
                for(j in 1:length(na.nodes)){
                    jind =  (tree@node$parent[last.level.set]==na.nodes[j])
                    pvalue.out[na.nodes[j]] = max(pvalue.out[last.level.set[jind]])
                }
                dec.level[na.nodes] = i-1
            }
            if(mode == 'stageI.in.two.stage') next;
            
            cum.rej = cum.rej + sum(!xind)
            this.level.set = this.level.set[ind]
            this.level.tree = this.level.tree[ind]
        }
        
        
        if(length(this.level.set)>0){
            if(i>1){
                this.level.p = sapply(this.level.tree,function(x)return(.StoufferCombination(pvalue.in[x])))
            }else{
                this.level.p = pvalue
            }
            
            weight = get_weight(this.level.set, tree)
            
            this.level.test.summary = .WeightedStepdownControl(this.level.p, weight, cum.rej, q.this.level, tau)
            this.level.detect = this.level.test.summary$output
            sig.out[this.level.set] = this.level.detect
            pvalue.out[this.level.set] = this.level.p
            dec.level[this.level.set] = i
            if(i ==16){
                print(this.level.p)
                print(cum.rej)
                print(this.level.test.summary)
            }
            
            
            cum.rej = cum.rej + sum(this.level.detect)
            last.cut.off = this.level.test.summary$cut.off
            alpha.out[i] = last.cut.off
            
            
            if(sum(!this.level.detect)>0){
                pvalue.in[this.level.set[!this.level.detect]]=
                    (this.level.p[!this.level.detect]-last.cut.off)/(1-last.cut.off)
            }
        }
        
    }
    
    tree@info$alpha = alpha.out
    tree@info$nl = nl.out
    tree@info$ql = q.spent
    
    tree@test$reject = sig.out; tree@test$pvalue = pvalue.out; tree@node$dec.level = dec.level
    if(mode=='stageI.in.two.stage'){#1001
        tree@test$num = last.cut.off
    }else{
        tree@test$num = far
    }
    return(tree)
}

.StoufferCombination<-function(pvalue){
    nan.pvalue = pvalue[!is.na(pvalue)]
    n = sum(nan.pvalue<=1)
    if(n>0){
        z.stat = sum(qnorm(1-nan.pvalue))/sqrt(n)
        return(1-pnorm(z.stat))
    }else{
        return(0)
    }
}

setClass("Tree", representation = list(node = "data.frame",
                                       edge = "list", test = "list", info = "list"))


setMethod("initialize", "Tree", function(.Object, ...) {
    value <- callNextMethod()
    value
})


.build.tree <-function(anno.table, na.symbol = 'unknown'){
    if(is.null(anno.table)) stop("Incorrect number of parameters.")
    n.rank = ncol(anno.table)
    num = id = xparent = parent=level= NULL
    count = 0
    
    #check if there is any ambiguous symbol (identical symbol for different nodes)
    tot.nodes = NULL
    for(i in n.rank:1){
        tot.nodes = c(tot.nodes, setdiff(unique(anno.table[,i]),na.symbol))
    }
    if(length(tot.nodes)!=length(unique(tot.nodes))){
        stop("At least two nodes from different levels have the same name.")
    }
    
    #Build a "Tree" object
    for(i in n.rank:1){
        this.level.node = setdiff(unique(anno.table[,i]),na.symbol)
        this.level.len = length(this.level.node)
        if(this.level.len > 0){
            id = c(id, this.level.node)
            num = c(num, c((count+1): (count+this.level.len)))
            count = count + this.level.len
            level = c(level, rep(n.rank+1-i, this.level.len))
            for(j in 1:this.level.len){
                t = which(anno.table[,i] == this.level.node[j])
                k = i
                while(k > 1){
                    ts = unique(anno.table[t,k-1])
                    if(length(ts) > 1){
                        stop("This is not a tree.")
                    }else{
                        if(ts != na.symbol){
                            break
                        }else{
                            k = k-1
                        }
                    }
                }
                if(i > 1){
                    xparent = c(xparent, ts)
                }else{
                    xparent = c(xparent, "head")
                }
            }
        }
    }
    names(num) = id
    for(i in 1:count-1){
        parent[i] = as.double(num[xparent[i]])
    }
    parent = c(parent, count+1)
    
    
    leaf.inner.cutoff = min(parent)
    edge.list =  vector("list", count)
    for(i in 1:count){
        if(i < leaf.inner.cutoff){
            edge.list[[i]] = integer(0)
        }else{
            edge.list[[i]] = which(parent==i)
        }
    }
    
    names(num) = NULL
    dec.level = rep(NA, length(level))
    
    tree = new(Class = "Tree", node = data.frame(num,id, parent,
                                                 level,dec.level, stringsAsFactors=FALSE), edge = edge.list,
               test= list(reject=logical(0),pvalue=numeric(0), num= numeric(0)),
               info = list(alpha = numeric(0), nl = numeric(0), ql = numeric(0))
    )
    return(tree)
    
}



setMethod("summary","Tree", function(object){
    node = slot(object, "node")
    node.id = node$id
    test = slot(object, "test")
    cat("Number of hypotheses:", length(node.id), "\n")
    if(length(test$reject)==0 || length(test$pvalue)==0){
        cat("None hypothesis testing has been done.", "\n")
    }else{
        
        
        result.to.print = data.frame(rev(node.id[test$reject]), formatC(rev(test$pvalue[test$reject]),format = "e",digits =2))
        colnames(result.to.print) = c("node","p.value")
        cat("Number of tree discoveries:", sum(test$reject), "\n")
        cat("FDR/FAR controlled by:", object@test$num, "\n")
        n.to.print <- min(nrow(result.to.print), 10)
        print(result.to.print[1:n.to.print,], row.names=FALSE)
        if(n.to.print < nrow(result.to.print)) {
            cat('[only 10 top-level significant hypotheses shown]', '\n')
        }
    }
})



setGeneric(".sig", function(object) standardGeneric(".sig"))


setMethod(".sig","Tree", function(object){
    node = slot(object, "node")
    test = slot(object, "test")
    
    if(length(test$reject)==0){
        return("None hypothesis testing has been done")
    }
    return(test$reject)
})




setGeneric(".top", function(object) standardGeneric(".top"))


setMethod(".top","Tree", function(object){
    node = slot(object, "node")
    test = slot(object, "test")
    
    if(length(test$reject)==0){
        return("None hypothesis testing has been done")
    }
    nn = length(node$id)
    discoveries = which(test$reject)
    nd = length(discoveries)
    p.values = signif(test$pvalue[test$reject],3)
    top.discoveries = array(TRUE, dim = nd)
    
    for(i in 1:nd){
        t = discoveries[i]
        while(t < nn){
            t = node$parent[t]
            if(test$reject[t]){
                top.discoveries[i] = FALSE
                break
            }
        }
    }
    result = array(FALSE, nn)
    result[discoveries[top.discoveries]]=TRUE
    return(result)
})




setGeneric(".query", function(object, x, na.symbol) standardGeneric(".query"))


setMethod(".query","Tree", function(object, x, na.symbol){
    node = slot(object, "node")
    test = slot(object, "test")
    
    if(length(test$reject)==0){
        return("None hypothesis testing has been done")
    }
    nn = length(node$id)
    nq = length(x); curr = NULL;
    name.vec = array(NA, nq)
    name.mat = matrix(NA, nrow = nq, ncol = max(node$level))
    for(i in 1:nq){
        curr = x[i]
        Path.v = Path.m= NULL
        while(curr <= nn){
            if(is.null(Path.v)){
                Path.v = node$id[curr]
                Path.m = c(node$id[curr], rep(" ", node$level[curr]-1))
            }else{
                Path.v = paste(node$id[curr], Path.v, sep=',')
                Path.m = c(node$id[curr], Path.m)
            }
            curr.next = node$parent[curr]
            if(curr < nn){
                level.curr = node$level[curr]; level.next =node$level[curr.next]
                if(level.next > level.curr + 1){
                    Path.v = paste(paste0(rep(na.symbol,level.next-level.curr-1),collapse=','), Path.v, sep = ',')
                    Path.m = c(rep(na.symbol, level.next-level.curr-1), Path.m)
                }
            }
            curr = curr.next
        }
        
        name.vec[i] =  Path.v
        name.mat[i,] = Path.m
        
    }
    return(list(name.vec = name.vec, name.mat = name.mat))
})

.UnweightedStepdownControl<-function(p.value, n.discovery, q.spend, tau = 0.3){ ### Proposed FDR method
    len.p.value = length(p.value)
    alpha.i = (n.discovery  + (1:len.p.value))/(len.p.value+1-(1:len.p.value))*q.spend
    alpha.i = pmin(alpha.i/(1+alpha.i),tau)   # thresholding cutoff values by tau
    index.x = which(sort(p.value)/alpha.i > 1)
    if(length(index.x) > 0){
        cut.off = alpha.i[min(index.x)]
        output = (p.value <= cut.off)
    }else{
        cut.off = tau
        output = rep(TRUE,len.p.value)
    }
    return(list(cut.off = cut.off, output = output))
}



.WeightedStepdownControl<-function(p.value, weight, n.discovery, q.spend, tau = 0.3){ ### Proposed FDR method
    len.p.value = length(p.value)
    alpha.i = array(NA, dim = len.p.value)
    for(i in 1:len.p.value){
        alpha.i[i] = (n.discovery + sum(weight[1:i]))/(sum(weight[i:len.p.value]))*q.spend
    }
    alpha.i = pmin(alpha.i/(1+alpha.i),tau)   # thresholding cutoff values by tau
    index.x = which(sort(p.value)/alpha.i > 1)
    if(length(index.x) > 0){
        cut.off = alpha.i[min(index.x)]
        output = (p.value <= cut.off)
    }else{
        cut.off = tau
        output = rep(TRUE,len.p.value)
    }
    return(list(cut.off = cut.off, output = output))
}

bouth<-function(anno.table, pvalue.leaves, na.symbol ='unknown', far = 0.10, tau = 0.3, is.weighted = TRUE){
    # this version supports one-stage & two-stage tests
    
    if(is.null(anno.table) || is.null(pvalue.leaves)) stop("Incorrect number of parameters")
    if(!is.null(colnames(anno.table))){
        level.name = rev(colnames(anno.table))
    }else{
        level.name = rev(paste("level", c(1:ncol(anno.table)), sep='-'))
    }
    if(is.null(names(pvalue.leaves))){
        names(pvalue.leaves) = anno.table[,ncol(anno.table)]
    }
    # build a tree object based on the input table
    tree = .build.tree(anno.table = anno.table, na.symbol = na.symbol)
    if(max(pvalue.leaves)==1){
        warning("There are p-values that equal exactly 1.")
        pvalue.leaves[pvalue.leaves==1] = 0.5*(1+max(c(0,pvalue.leaves[pvalue.leaves<1])))
    }
    
    # one-stage approach
    if(length(far)==1){
        if(is.weighted){
            results = .bottom.up.weighted(tree, pvalue.leaves, far, tau)
        }else{
            results = .bottom.up.unweighted(tree, pvalue.leaves, far, tau)
        }
    }
    # two-stage approach
    if(length(far)==2){
        results = .two.stage.test(anno.table, pvalue.leaves, na.symbol, far.leaves = far[1], far.inner = far[2],
                                  tau)
        
    }
    
    
    options(digits=3)
    
    ## create a table for all nodes tested in a bottom-up procedure
    n.node = length(results@node$level); max.level = max(results@node$level)
    show.node = .query(results, c(1:n.node), na.symbol); colnames(show.node$name.mat) = rev(level.name)
    show.pvalue = results@test$pvalue
    show.detected = .sig(results)
    show.driver =  .top(results)
    results.by.node = data.frame(show.node$name.mat, show.pvalue, show.detected, show.driver,stringsAsFactors = FALSE)
    colnames(results.by.node) = c(rev(level.name),"pvalue", "is.detected", "is.driver")
    results.by.node = results.by.node[order(show.node$name.vec),]
    rownames(results.by.node) = NULL
    
    ## create a table to save characteristics for the bottom-up procedure
    detections.per.level = sapply(c(1:max.level),function(x)return(sum(results@node$level[results@test$reject]==x)))
    results.by.level = data.frame(c(1:max.level),level.name,
                                  results@info$ql, results@info$alpha,
                                  results@info$nl, detections.per.level,
                                  stringsAsFactors=FALSE)
    colnames(results.by.level) = c("level","level.name", "q_l", "pvalue.cutoff", "all.nodes", "detected.nodes")
    
    
    return(list(tree = results, results.by.node = results.by.node, results.by.level = results.by.level))
}