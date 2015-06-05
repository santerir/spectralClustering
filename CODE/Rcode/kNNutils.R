library("parallel")


## -- returns indices of neighbours and distances to them -- ##

getKNearestNeighbors <- function(data,k,param=0) {
    if(param==0) {
      parame <- rep(1,NCOL(data));
    }
    param <<- param;


    sp_tree <- constructSpTree(data,max(k+10,20));
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
    a <- LInfMetric(sp_tree[[level]]$dat,Vector);    
    ord <- order(a);
    copy <- sp_tree[[level]][ord,];
  
    if (is.data.frame(copy)) {
      return(cbind(dist=a[ord[1:k]],copy[1:k,]));
    }
    else {
      print("huzzaah! a wild error appears!");
      browser();
    }
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
    a <- LInfMetric(sp_tree[[level]]$dat,Vector);
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

constructSpTree <- function(data, b_size)  {
  df <- data.frame(ind=1:NROW(data),dat=rep(0,NROW(data)));
  df$dat <- data;

  sp_tree <- list(df)

  i <- 1;
  
  while(1==1) {
    
    p <- partition(sp_tree[[i]]$dat);

    p <- list(p);
    if (is.matrix(data)) {
      
      pr <- projections(p[[1]],sp_tree[[i]]$dat);
      
      SU <- summary(pr);
      sp_tree[[2*i]] <- sp_tree[[i]][which(pr>=SU[3]),];
      sp_tree[[(2*i)+1]] <- sp_tree[[i]][which(pr<SU[3]),];
      
      p[[2]] <- SU[3];
      sp_tree[[i]] <- p;

      
      if ((log2(2*(i+1)))==floor(log2(2*(i+1))) && (NROW(sp_tree[[2*i]])<=b_size)){
        return(sp_tree);;                                                                      
      }
      i<-i+1;
    }
  }
}



LInfMetric <- function(vector1,vector2) {

  resultV <- c();
  lim <- 1;
  if (is.matrix(vector1)) {
    for(i in (param-1)) {
      resultV <- cbind(resultV,L2Metric(vector1[,lim:(lim+i)],vector2[lim:(lim+i)]));
      lim <- lim + i + 1;
    }
  }
  else  {
    for(i in param) {
      resultsV <- cbind(resultV,L2Metric(vector1[lim:(lim+i)],vector2[lim:(lim+i)]));
      lim <- lim + i + 1;
    }
  }
  if (is.matrix(resultV)) {
    return(apply(resultV,1,max));
  }
  else {
    return(max(resultV));
  }
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
