library("parallel")

getKnnRadii <- function(data,param=0,k) {
  ## -- if Data is 1-dimensional
  if(NCOL(data) == 1) {
     radii <- (getKnnRadii1D(data,k));
    return(unlist(radii));
	}

  param <<- param;


	sp_tree <- constructSpTree(data,20);

  cores <- detectCores();
  limit <- log2(length(sp_tree)+1)-1;


  res <- mclapply(apply(data, 1, FUN=list),function(x) singleRadii(unlist(x),sp_tree,k,limit),mc.cores=cores);
  if (NROW(param)==1) {
    return(unlist(res));
  } 
  return(t(sapply(res, unlist)));
}



singleRadii <- function(x,sp_tree,k,limit) {
  candidateSet <- getKnn(x,sp_tree,k,limit);
  res <- c();
  lim <- 1;
  if (is.matrix(candidateSet)){
    for (i in (param-1)) {
      a <- L2Metric(candidateSet[,lim:(lim+i)],x[lim:(lim+i)]);
      lim <- lim+i+1;
      res <- c(res,max(a));
    }
  }
  else {
    for (i in (param-1)) {
      a <- L2Metric(candidateSet[lim:(lim+i)],x[lim:(lim+i)]);
      lim <- lim+i+1;
      res <- c(res,max(a));
    }
  }
  return(res);
}


getKNearestNeighbors <- function(data,param=0,k) {
    param <<- param;


    sp_tree <- constructSpTree(data,max(k+1,20));
    cores <- detectCores();
    limit <- log2(length(sp_tree)+1)-1
    res <- mclapply(apply(data, 1, FUN=list),function(x) as.vector(t(getKnn(unlist(x),sp_tree,k,limit))),mc.cores=cores);

    return(t(sapply(res, unlist)));
}

getKNearestNeighbor <- function(data,param=0,k) {
    param <<- param;


    sp_tree <- constructSpTree(data,max(k+1,20));
    cores <- detectCores();
    limit <- log2(length(sp_tree)+1)-1
    res <- mclapply(apply(data, 1, FUN=list),function(x) getKnn(unlist(x),sp_tree,k,limit)[k,],mc.cores=cores);

    return(t(sapply(res, unlist)));
}


getKnn <- function(x, sp_tree, k, limit) {
  searchBound <- recursiveSearchBound(sp_tree,x,1,limit,k);

  candidateSet <- recursiveSearch(sp_tree,x,k,1,limit,searchBound)[,2:(NROW(x)+1)];

  return(candidateSet);
}


countPointsInRadius <- function(data,param,radii) {
  param <<- param;
  if (NCOL(data)==1) {
    return(unlist(countPointsInRadius1D(data,radii)));
  }

  spTree <- constructSpTree(data,50);
  limit <- log2(length(spTree)+1)-1

  cores <- detectCores();

  res <- mclapply(apply(cbind(data,radii), 1, FUN=list), function(x) recursiveCount(spTree,
    unlist(x)[1:NCOL(data)],
    unlist(x)[(NCOL(data)+1):NROW(unlist(x))],
    1,
    limit),
    mc.cores=cores);

return(unlist(res));

}

recursiveCount <- function(sp_tree,Vector,radius,level,limit)  {
  if (level>=(2^limit)){
    res <- c()
    lim <- 1;
    iterator <- 1;
    for(i in (param-1)) {
      res <- cbind(res,as.numeric(L2Metric(sp_tree[[level]][,lim:(lim+i)],Vector[lim:(lim+i)])<=radius[iterator]));
      lim <- lim + i + 1; 
      iterator <- iterator + 1;
    }
    return(sum(apply(res,1,FUN=prod)));
  }

  r = projections(sp_tree[[level]][[1]],Vector)
  d = (L2Metric(r,sp_tree[[level]][[2]]));
  
  a <- 0;
  b <- 0;

  if (r>=sp_tree[[level]][[2]] || (1.2*(sqrt(sum(radius^2)))) > d) {
    a <- recursiveCount(sp_tree,Vector,radius,2*level,limit);
  }
  if (r<sp_tree[[level]][[2]] || (1.2*(sqrt(sum(radius^2)))) > d) {
    b <- recursiveCount(sp_tree,Vector,radius,2*level+1,limit);
  }
  
  return(a+b);
}

countPointsInRadius1D <- function(data,radii) {
  data <- cbind(data,radii);
  ord <- order(data[,1]);
  orderedData <- cbind(1:NROW(data),data[ord,]);
  cores <- detectCores();
  nPoints <- mclapply(apply(orderedData, 1, FUN=list), function(x,dataOrdered) {
        x <- unlist(x);
        radius <- x[3];
        vector <- x[2];
        index <- x[1];
        iterator <- index;
        n <- 0;
        while(iterator > 0) {
          if(L2Metric(orderedData[iterator,2],vector)<=radius) {
            n <- n+1;
            iterator <- iterator - 1;
          }
          else {
            break;
          }
        }
         iterator <- index+1;
         while(iterator < NROW(orderedData)+1) {
          if(L2Metric(orderedData[iterator,2],vector)<=radius) {
            n <- n+1;
            iterator <- iterator + 1;
          }
          else {
            break;
          }
        }
        return(n);


    },mc.cores=cores);
  nPoints[order(ord)];
  return(nPoints);
}



getKnnRadii1D <- function(data,k) {
	cores <- detectCores()
	dataOrdered <- data[order(data)];
	radii <- mclapply(data, function(x){
				y <- which(dataOrdered==x);
				a <- 100000;
				b <- 100000;
				if(y>=k) {
					a <- L2Metric(x,dataOrdered[y-(k-1)]);
				}
				if (y<=NROW(dataOrdered)-k) {
					b <- L2Metric(x,dataOrdered[y+(k-1)])
				}

				return(min(a,b));	
				
				},mc.cores=cores)
}


recursiveSearch <- function(sp_tree,Vector,k,level,limit,s_b)  {
  if (level>=(2^limit)){
    a <- LInfMetric(sp_tree[[level]],Vector);    
    ord <- order(a);
    copy <- sp_tree[[level]][ord,];
  
    if (is.matrix(copy)) {
      return(cbind(a[ord[1:k]],copy[1:k,]));
    }
    else {
      return(cbind(a[ord[1:k]],copy[1:k]));
    }
  }
  
  r = projections(sp_tree[[level]][[1]],Vector)
  d = (L2Metric(r,sp_tree[[level]][[2]]));
  a <- c(-1);
  b <- c(-1);

  
  if (r>=sp_tree[[level]][[2]] || 2*s_b>d) {
    a <- recursiveSearch(sp_tree,Vector,k,2*level,limit,s_b);
  }
  if (r<sp_tree[[level]][[2]] || 2*s_b>d) {
    b <- recursiveSearch(sp_tree,Vector,k,2*level+1,limit,s_b);
  }
  
  if (!is.matrix(a)) {
    return(b);
  }
  if (!is.matrix(b)) {
    return(a);
  }
  else {
    ret <- rbind(a,b);
    ret <- ret[order(ret[,1]),];
    ret <- ret[1:k,];
    return(ret);
  }
}



recursiveSearchBound <- function(sp_tree,Vector,level,limit,k)  {
  if (level>=(2^limit)){
    a <- LInfMetric(sp_tree[[level]],Vector);
    a <- a[order(a)];
    if(min(a)==0) {
      return(a[k]);
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


constructSpTree <- function(data, b_size)  {      
  sp_tree <- list(data)
  i <- 1;
  
  while(1==1) {
    
    p <- partition(sp_tree[[i]]);

    p <- list(p);
    if (is.matrix(data)) {
      
      pr <- projections(p[[1]],sp_tree[[i]]);
      
      SU <- summary(pr);
      sp_tree[[2*i]] <- rbind(sp_tree[[i]][which(pr>=SU[3]),]);
      sp_tree[[(2*i)+1]] <- rbind(sp_tree[[i]][which(pr<SU[3]),]);
      
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
