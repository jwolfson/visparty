
#Gets node ID for all the left splits
ptree_left <- function(newtree, start_id){
  nodes(newtree, start_id)[[1]]$left$nodeID
}

#Gets node ID for all the right splits
ptree_right <- function(newtree, start_id){
  nodes(newtree, start_id)[[1]]$right$nodeID
}

# Get prediction for the relevant node
ptree_y <- function(newtree, node_id) {
  # Picks out the prediction from the ctree structure
  p <- nodes(newtree, node_id)[[1]]$prediction
  
  if (length(p)==2) {
    return(p[2])
  }
  return (p)
}

# criteria for this node as a string
# If this is a terminal node, then error returned as being a terminal node
# If this is an ordered variable, use the primary split and then run through the tree bringing out the left and right splits

ptree_criteria <- function(newtree, node_id, left) {
  if (nodes(newtree, node_id)[[1]]$terminal) # Check if this is a terminal node
  {
    return("(error: terminal node)");
  } 
  if (nodes(newtree, node_id)[[1]]$psplit$ordered)
  {
    sp <- nodes(newtree, node_id)[[1]]$psplit$splitpoint
    vn <- nodes(newtree, node_id)[[1]]$psplit$variableName
    # Left being true then the left sting of variables with split points are returned 
    if (left) {
      op <- '<='   
    } else {
      op <- '>'
    }
    return(paste(vn, op, sp))
  } else {
    psplit <- nodes(newtree, node_id)[[1]]$psplit
    if (left){
      l <- as.logical(psplit$splitpoint)
    } else {
      l <- as.logical(!psplit$splitpoint)
    }
    
    r <- paste(attr(psplit$splitpoint, 'levels')[l], sep='', collapse="','")
    return(paste(psplit$variableName, " in ('", r,"')", sep=''))
  }
}

list_node <- function(newtree, node_id = 1, start_criteria = character(0)) {
  if (nodes(newtree, node_id)[[1]]$terminal) {
    prediction <- ptree_y(newtree, node_id)
    ypred <- paste( start_criteria, ',y =',prediction,';')
    return (ypred)
  }
  
  left_node_id <- ptree_left(newtree, node_id)
  right_node_id <- ptree_right(newtree, node_id)
  
  if (is.null(left_node_id) != is.null(right_node_id)) {
    print('left node ID != right node id')
  }
  ypred <- character(0)
  if (!is.null(left_node_id)) {
    new_criteria <- paste(start_criteria, ptree_criteria(newtree, node_id, T), sep=',')
    if (1 == node_id)
      new_criteria <- ptree_criteria(newtree, node_id, T)
    ypred <- list_node(newtree, left_node_id, new_criteria)
  }
  if (!is.null(right_node_id)) {
    new_criteria <- paste(start_criteria, ptree_criteria(newtree, node_id, F), sep=',')
    if (1 == node_id)
      new_criteria <- ptree_criteria(newtree, node_id, F)
    ypred <- paste(ypred, list_node(newtree, right_node_id, new_criteria))
  }
  return(ypred)
}

trim <- function (x) gsub("^\\s+|\\s+$", "", x) ## Function to get rid of spaces in strings

makeTransparent = function(..., alpha=0.5) {
  ## Helper function to make colors transparent
  if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
  alpha = floor(255*alpha)  
  newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
  .makeTransparent = function(col, alpha) {
    rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
  }
  newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
  return(newColor)
}

minmax.mat <- function(str,varnms) {
  ## Helper function to create a matrix of ranges for each variable in a path to a node
  comps <- strsplit(str,",")
  MMM <- matrix(data=rep(c(-Inf,Inf),length(varnms)),nrow=length(varnms),ncol=2,byrow=TRUE) ### min-max matrix
  rownames(MMM) <- varnms
  for(i in 1:(length(comps[[1]])-1)) {
    nodestr <- strsplit(trim(comps[[1]][i])," ")
    node.varnm <- trim(nodestr[[1]][1])
    node.dir <- trim(nodestr[[1]][2])
    node.split <- trim(nodestr[[1]][3])
    var.row <- which(varnms==node.varnm)
    if(node.dir == "<=") {
      MMM[var.row,2] <- as.numeric(node.split)
    }
    else {
      MMM[var.row,1] <- as.numeric(node.split)
    }
  }
  y <- comps[[1]][length(comps[[1]])]
  return(list(M=MMM,y=y))
}

plot.minmax <- function(My,X,Y) {
  ## Main function which plots the bars for each variable along with a histogram of the outcome
  mymat <- My$M
  my.y <- My$y
  
  if(!is.factor(Y)) {
    my.y.val <- as.numeric(strsplit(trim(my.y)," ")[[1]][3])
    my.y.pct <- ecdf(Y)(my.y.val)
  }
  
  var.nms <- rownames(mymat)
  act.vars <- apply(mymat,1,function(x) { !all(abs(x)==Inf) })
  
  max.y <- sum(act.vars)+1  
  rbw <- rainbow(n=nrow(mymat))
  
  ## Find the y's which "belong" in this node
  node.index <- 1:length(Y)
  for(i in 1:nrow(mymat)) {
    node.index <- intersect(node.index,which(X[,i]>mymat[i,1]&X[,i]<=mymat[i,2]))
  }
  
  ## Create the underlying histogram, but don't plot it yet
  if(is.factor(Y)) {
   wdth <- 1/length(levels(Y))
    H <- hist(as.integer(Y[node.index])/length(levels(Y)),breaks=seq(0,1,length.out=length(levels(Y))+1),plot=FALSE)
   ## Scale the histogram so it fits vertically on the plot.
   scale.factor <- max.y/max(H$density)
   ## Set up an empty plot of the correct size
   plot(NA,xlim=c(0,1),ylim=c(0,max.y),ylab="",xlab="Percentile",main=paste("Mode =",names(which.max(table(Y[node.index]))),", n =",length(node.index)),bty="n",yaxt="n")
   ## Plot the background histogram
    barplot(scale.factor*H$density,width=wdth,col=rgb(0,0,0,0.05),border=rgb(0,0,0,0.1),add=TRUE,space=0,yaxt="n")
   ## Add the category labels
   text(seq(wdth/2,1-wdth/2,by=wdth),rep(0,length(levels(Y))),levels(Y),pos=3,adj=0.5,cex=1,col=gray(0.5))
  }
  else{
    H <- hist(ecdf(Y)(Y[node.index]),breaks=seq(0,1,by=0.1),plot=FALSE)
    ## Scale the histogram so it fits vertically on the plot.
    scale.factor <- max.y/max(H$density)
    ## Set up an empty plot of the correct size
    plot(NA,xlim=c(0,1),ylim=c(0,max.y),ylab="",xlab="Percentile",main=paste("Mean =",signif(my.y.val,2),", n =",length(node.index)),bty="n",yaxt="n")
    ## Plot the background histogram
    barplot(scale.factor*H$density,width=0.1,col=rgb(0,0,0,0.05),border=rgb(0,0,0,0.1),add=TRUE,space=0,yaxt="n")
    ## Draw in a line for the mean
    mu.Y <- mean(Y[node.index])
    segments(ecdf(Y)(mu.Y),0,ecdf(Y)(mu.Y),max.y,col=rgb(0,0,0,0.5),lwd=2)    
  }
  
  ## Now plot the horizontal bars corresponding to each variable.
  j <- 1
  for(i in which(act.vars)) {
    F.x <- ecdf(X[,var.nms[i]])
    lo <- ifelse(mymat[i,1]==-Inf,0,F.x(mymat[i,1]))
    hi <- ifelse(mymat[i,2]==Inf,1,F.x(mymat[i,2]))
    polygon(c(lo,lo,hi,hi),c(j-0.5,j+0.5,j+0.5,j-0.5),col=makeTransparent(rbw[i],0.5),border=NA)
    text(mean(c(lo,hi)),j,rownames(mymat)[i]) ## Label the variables
    j <- j+1
  }
  
}

visTree <- function(cond.tree,rng=NULL) {
  ## Wrapper function to produce plots from a conditional inference tree
  ## 'range' parameter can restrict plotting to a particular set of nodes
  splittree<-list_node(cond.tree)
  structure<-strsplit(splittree, split=";")
  if(is.factor(Y)) {
    n.terminals <- length(structure[[1]])
    prob.mat <- matrix(data=unlist(lapply(structure,function(S) {
      
      unlist(lapply(strsplit(S,","),function(split.S) {
        seg <- split.S[length(split.S)]
        as.numeric(trim(strsplit(seg,"=")[[1]][2]))
      })) 
      })), nrow=n.terminals)
    
    y.list <- lapply(1:nrow(prob.mat),function(j) {
      paste0("y=",paste0(prob.mat[j,],collapse="|"))
    })

    x.list <- sapply(structure[[1]],function(s.row) {
      seg <- strsplit(s.row,",")[[1]]
      paste0(seg[-length(seg)],collapse=",")
    })
    
    structure <- lapply(1:length(y.list),function(i) {
      paste0(x.list[[i]],", ",y.list[[i]])
    })
  }
  if(length(unlist(structure))==1) { stop("Tree has only a single node; nothing to visualize.") }

  n.terminals <- ifelse(is.null(rng),length(unlist(structure)),length(rng))
  if(is.null(rng)) { 
    index <- 1:n.terminals } else { 
    index <- min(rng):min(max(rng),length(unlist(structure))) } ## Should probably do some range checking
  par(mfrow=c(2,ceiling(length(index)/2)),mar=c(2,1,3,1))
  X <- data.frame(cond.tree@data@get("input"))
  Y <- (cond.tree@data@get("response"))[,1]
  
  sapply(unlist(structure)[index],function(S) { plot.minmax(minmax.mat(S,colnames(X)),X,Y)})
}

testing <- function() {
  X1 <- rnorm(100)
  X2 <- rnorm(100)
  Y <- cut(rnorm(100,mean=X1+X2,sd=1),4)
  
  cond.tree <- party::ctree(Y~X1+X2)
  visTree(cond.tree)
}
