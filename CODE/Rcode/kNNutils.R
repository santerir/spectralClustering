library("parallel")


## -- returns indices of neighbours and distances to them -- ##

getKNearestNeighbors <- function(data,k) {

    if (k>=NROW(data)/2+5) {
      print("ERROR: this would return over half the datapoints as neighbours. I cant do that for you");
      return(-1);
    }

    b <- ceiling(log2(NROW(data)/(k+5))-1);

    sp_tree <- constructSpTree(data,b);

    cores <- detectCores();
    limit <- log2(length(sp_tree)+1)-1
    res <- mclapply(apply(data, 1, FUN=list),function(x) getKnn(unlist(x),sp_tree,k,limit),mc.cores=cores);
    
    result <- matrix(0,NROW(data),k*2);

    for(i in 1:length(res)) {
      result[i,] <- as.vector(c(res[[i]]$ind,res[[i]]$dist));
    }
    return(result);
}


getKnn <- function(x, sp_tree, k, limit) {
  searchBound <- recursiveSearchBound(sp_tree,x,1,limit,k);

  candidateSet <- recursiveSearch(sp_tree,x,k,1,limit,searchBound);

  return(candidateSet);
}


recursiveSearch <- function(sp_tree,Vector,k,level,limit,s_b)  {
  if (level>=(2^limit)){
    a <- L2Metric(sp_tree[[level]]$dat,Vector);    
    ord <- order(a);
  
    return(data.frame(cbind(dist=a[ord[1:k]],sp_tree[[level]][ord[1:k],])));
  }
  
  r = projections(sp_tree[[level]][[1]],Vector)
  d = (L2Metric(r,sp_tree[[level]][[2]]));

  a <- 0
  b <- 0

  if (r>=sp_tree[[level]][[2]] || 2*s_b>d) {
    a <- recursiveSearch(sp_tree,Vector,k,2*level,limit,s_b);
  }
  if (r<sp_tree[[level]][[2]] || 2*s_b>d) {
    b <- recursiveSearch(sp_tree,Vector,k,2*level+1,limit,s_b);
  }

  if (!is.data.frame(a)) {
    return(b);
  }

  else if (!is.data.frame(b)) {
    return(a);
  }
  
  else {
    ret <- rbind(a,b);
    ret <- ret[order(ret$dist),];
    ret <- ret[1:k,];
    return(ret);
  }
}



recursiveSearchBound <- function(sp_tree,Vector,level,limit,k)  {
  if (level>=(2^limit)){
    a <- L2Metric(sp_tree[[level]]$dat,Vector);
    a <- a[order(a)];
    if(min(a)==0) {
      return(a[k+1]);
    }
    else {
      return(-1);
    }
  }
  r = projections(sp_tree[[level]][[1]],Vector)
  a <- b <- -1
  if (r>=sp_tree[[level]][[2]]) {
    a <- recursiveSearchBound(sp_tree,Vector,2*level,limit,k);
  }
  if (r<sp_tree[[level]][[2]]) {
    b <- recursiveSearchBound(sp_tree,Vector,2*level+1,limit,k);
  }
  
  return(max(a,b))
}

## -----------------------
## -- SUPPORT FUNCTIONS --
## -----------------------


constructSpTree <- function(data, b)  {
  df <- data.frame(ind=1:NROW(data),dat=rep(0,NROW(data)));
  df$dat <- data;

  sp_tree <- list(df)

  limit <- 2^(b-1)-1;

  for(i in 1:limit) {
    
    p <- partition(sp_tree[[i]]$dat);

    p <- list(p);

      
    pr <- projections(p[[1]],sp_tree[[i]]$dat);
      
    SU <- summary(pr);
    sp_tree[[2*i]] <- sp_tree[[i]][which(pr>=SU[3]),];
    sp_tree[[(2*i)+1]] <- sp_tree[[i]][which(pr<SU[3]),];
      
    p[[2]] <- SU[3];
    sp_tree[[i]] <- p;

  }
  return(sp_tree)
}



partition <- function(data) {

  A <- L2Metric(data,data[1,]);
  a <- data[which.max(A),];
  A <- L2Metric(data,a);
  b <- data[which.max(A),];

  a <- a - b;

  return(a);   
}




#Calculates scalar projections
projections <- function(pr_v, data) {
  k <- sqrt(sum(pr_v*pr_v));
  if (is.matrix(data)) {
    return(apply(t(t(data)*pr_v),1,sum)/k);
  }
  return(sum(data*pr_v)/k);
}



L2Metric <- function(a, b) {        
  if (NROW(a)==NROW(b)&&!is.matrix(a)&&!is.matrix(b)) {
    
    return(sqrt(sum((a-b)^2)));
  }
  else {
    
    return(sqrt(apply(t(t(a)-b)^2,1,sum)));
  }
}
