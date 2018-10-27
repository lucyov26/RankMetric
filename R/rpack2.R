#' @useDynLib RankMetric, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats as.dist
#' @importFrom LIStest lis
NULL





#' @rdname RankMetric
#' @title Metrics for Rankings
#' @description Calculates the distances between a set of rankings or permutations,
#' using one of six metrics.
#' @author Lucy Small, \email{lucy.small@@ucdconnect.ie}
#'
#' @param x matrix where each row is a ranking or a permutation
#' @param metric a distance metric, one of
#' \itemize{
#' \item \code{"kendall"}
#' \item \code{"ulam"}
#' \item \code{"spearrho"},  spearman's rho
#' \item \code{"spearfoot"}, spearman's footrule
#' \item \code{"hamming"}
#' \item \code{"cayley"}
#' }
#' @param perm \code{TRUE} for a matrix of permutations, \code{FALSE} for a matrix of rankings
#' @param ranktype indicates the type of ranking
#'    \itemize{
#'    \item \code{"full"} full ranking
#'    \item \code{"partial"},  partial ranking with no ties
#'    \item \code{"sametimes"}, ties allowed, the same number must be ranked 1st, 2nd etc
#'    \item \code{"anyties"}, ties allowed any number can be ranked 1st, 2nd etc
#'}
#' @param k the number of items ranked in a partial ranking
#' @return Returns a matrix of the distances between all the rankings/ permutations
#' in the inputted matrix. If \code{ranktype} is \code{"anyties"} the function will work for any
#' type of partial ranking or ranking with ties. However using \code{"full"} for complete rankings
#' or \code{"partial"} for partial ranking with no ties for that kind of data will run faster.
#'
#' If the ranking of some of the items is not known or given they can be left as \code{NA}.

#' @references \url{http://www.springer.com/gp/book/9780387962887}

#' @examples
#' x = t(matrix(replicate(10,sample(1:5,5)),ncol=10))
#' distance(x,metric = "spearfoot",perm = FALSE,"full")
#' @export
distance = function(x, metric, perm = TRUE, ranktype  = "full", k=2){



  d1 = dim(x)[1]
  d2 = dim(x)[2]

  if(ranktype == "full"){
    if( metric == "kendall"){
      f = kend
    }
    if(metric == "ulam"){
      f = ulam
    }

    if(metric == "spearrho"){
      f = spear
    }

    if(metric == "spearfoot"){
      f = spearfoot
    }

    if(metric == "hamming"){
      f = ham
    }

    if(metric == "cayley"){
      f = cayley
    }
  }

  if(ranktype == "partial"){
    if( metric == "kendall"){
      f = kendP
    }
    if(metric == "ulam"){
      f = ulamP
    }

    if(metric == "spearrho"){
      f = spearP
    }

    if(metric == "spearfoot"){
      f = spearfootP
    }

    if(metric == "hamming"){
      f = hamP
    }

    if(metric == "cayley"){
      f = cayleyP
    }
  }

  dis = matrix(NA,d1,d1)

  if(ranktype == "full"){

    if(perm == F){
      if(metric == "kendall"||metric == "spearfoot"||metric == "spearrho"||metric == "ulam"){
        x = t(apply(x,1,inv))
      }
    }


         for(i in 1:d1){
           for(j in 1:i){
             dis[i,j] = f(x[i,],x[j,])
           }
         }
  }

  if(ranktype == "partial"){

    if(perm==F){
      x = t(apply(x,1,inv))
    }

    for(i in 1:d1){
      for(j in 1:i){
        dis[i,j] = f(x[i,],x[j,],k)
      }
    }

    # fills NA values
    x[is.na(x)] =  rep((k+1):d2,each=d1)

  }



  if( ranktype == "sameties"|| ranktype == "anyties"){
    if(perm == F){
      print("For this distance measure, the data must be in permutation form, with each entry denoting the group membership of each item")
      stop
    }

    if(ranktype == "sameties"){
      if( metric == "kendall"){
        f = kendE
      }
      if(metric == "ulam"){
        f = ulamG
      }

      if(metric == "spearrho"){
        f = spearE
      }

      if(metric == "spearfoot"){
        f = spearfootE
      }

      if(metric == "hamming"){
        f = hamE
      }

      if(metric == "cayley"){
        f = cayleyE
        print("only works for groups of size 2 or 3")
        return()
      }
    }

    if(ranktype == "anyties"){
      if( metric == "kendall"){
        f = kendG
      }
      if(metric == "ulam"){
        f = ulamG
      }

      if(metric == "spearrho"){
        f = spearG
      }

      if(metric == "spearfoot"){
        f = spearfootG
      }

      if(metric == "hamming"){
        f = hamG
      }

      if(metric == "cayley"){
        print("no cayley's distance function for this kind of ranking")
      }
    }

    #fill NA values

    ind = which(is.na(x), arr.ind=TRUE)
    m=apply(x,1,max,na.rm = T)+1
    x[ind] <- m[ind[,1]]


    for(i in 1:d1){
      for(j in 1:i){
        dis[i,j] = f(x[i,],x[j,])
      }
    }
  }




  dis = as.dist(dis)
  return(dis)
  }





#some preliminary functions



#========================================

#' @rdname inv
#' @title Inverse Permutation
#' @description Computes the inverse of a permutation.
#' @author Lucy Small, \email{lucy.small@@ucdconnect.ie}
#'
#' @param x an integer vector
#' @return Returns the inverse permutation of a vector.

#' @examples
#' a = c(3,1,2,5,4)
#' inv(a)
#' @export
inv = function(x){
  names(x) = 1:length(x)
  as.numeric(names(sort(x)))
}

#' @rdname nig3
#' @title nig3
#' @description computes..
#' @author Lucy Small, \email{lucy.small@@ucdconnect.ie}
#'
#' @param a,b integer vectors
#' @keywords internal

nig3 = function(a,b){

  la= length(a)
  ma = max(a)
  mb = max(b)
  nig3=matrix(NA,ma,mb)

  for(i in 1:ma){
    for(j in 1:mb){
      nig3[i,j] = nij(a,b,i,j,la)
    }
  }
  return(nig3)
}


#Code below computes various measures of distance for
#=====================================================================
#Spearman's rho


#' @rdname spear
#' @title Spearman's Rho
#' @description Computes Spearman's rho between two full permutations.
#' @author Lucy Small, \email{lucy.small@@ucdconnect.ie}
#'
#' @param a,b integer vectors
#' @return Returns Spearman's rho between the two permutations.

#' @examples
#' a = c(3,1,2,5,4)
#' b = 1:5
#' spear(a,b)
#' @export
spear = function(a,b){
  sqrt(sum((a-b)^2))
}

#=====================================================================
#Spearman's Footrule
#' @rdname spearfoot
#' @title Spearman's Footrule
#' @description Computes Spearman's Footrule between two full permutations.
#' @author Lucy Small, \email{lucy.small@@ucdconnect.ie}
#'
#' @param a,b integer vectors
#' @return Returns Spearman's footrule between the two permutations.

#' @examples
#' a = c(3,1,2,5,4)
#' b = 1:5
#' spearfoot(a,b)
 #' @export
spearfoot = function(a,b){
  sum(abs(a-b))
}


#' @rdname kend
#' @title Kendall's Tau
#' @description Computes Kendall's tau between two full permutations.
#' @author Lucy Small, \email{lucy.small@@ucdconnect.ie}
#'
#' @param a,b integer vectors
#' @return Returns Kendall's tau between two full permutations.

#' @examples
#' a = c(3,1,2,5,4)
#' b = 1:5
#' kend(a,b)
#' @export

kend = function(a,b){
  a=inv(a)
  b=inv(b)
  n = length(a)
  return(kendR(a,b,n))
}

#=====================================================================
#Computes Hamming distance of perms a,b
#' @rdname ham
#' @title Hamming Distance
#' @description Returns the Hamming distance between two full permutations.
#' @author Lucy Small, \email{lucy.small@@ucdconnect.ie}
#'
#' @param a,b integer vectors
#' @return The Hamming distance between the two permutations.

#' @examples
#' a = c(3,1,2,5,4)
#' b = 1:5
#' ham(a,b)
#' @export

ham = function(a,b){
  c = sum(a != b)
  return(c)
}

#=====================================================================
# computes Ulam's Distance for perms a,b
#' @rdname ulam
#' @title Ulam's Distance
#' @description Computes Ulam's distance between two full permutations.
#' @author Lucy Small, \email{lucy.small@@ucdconnect.ie}
#'
#' @param a,b integer vectors
#' @return Returns Ulam's distance between the two permutations.

#' @examples
#' a = c(3,1,2,5,4)
#' b = 1:5
#' ulam(a,b)
#' @export
ulam = function(a,b){
  u = length(a)-lis(pmult(b,inv(a)))
  return(u)
}



# returns the partial ranking metrics from pg 18 (25)

#=====================================================================
#Computes Hamming distance of perms a,b for partial ranking
#' @rdname hamP
#' @title Hamming Distance for Partial permutations
#' @description Computes the Hamming distance between two partial permutations.
#' @author Lucy Small, \email{lucy.small@@ucdconnect.ie}
#' @param a,b integer vectors
#' @param k integer
#' @return Returns the Hamming distance between the two partial permutations,
#' where only the first \code{k} items have been ranked. No ties are permitted.

#' @examples
#' a = c(3,1,2,5,4)
#' b = 1:5
#' k=3
#' hamP(a,b,k)
#' @export
hamP = function(a,b,k){
  a1 = a[1:k]
  b1 = b[1:k]
  c = sum(a1 != b1)
  h = length(a1[!(a1 %in% b1)])

  return(c+h)
}

#=====================================================================
#Computes Kendall's distance of perms a,b for partial permutations
#' @rdname kendP
#' @title Kendall's Distance for Partial permutations
#' @description Computes Kendall's distance between two partial permutations.

#' @author Lucy Small, \email{lucy.small@@ucdconnect.ie}
#' @param a,b integer vectors
#' @param k integer
#' @return Returns Kendall's distance between two permutations x and y, where we care about the first k ranked items only

#' @examples
#' a = c(3,1,2,5,4)
#' b = 1:5
#' k=3
#' kendP(a,b,k)
#' @export
kendP = function(a,b,k){
  n = length(a)
  v = 1:n

  A = v[a<=k &b<=k]
  B = v[a<=k &b>k]
  D = v[a>k &b<=k]

  if(k == 1){p1=0}else{p1 = kend(a[A],b[A])}

  h = length(B)
  p4 = h*(n+k-(h-1)/2)

  return(p1 +p4 - sum(a[B]) - sum(b[D]))
}

#=====================================================================
# Spearman's footrule
#' @rdname spearfootP
#' @title Spearman's Footrule for Partial permutations
#' @description Computes Spearman's footrule between two partial permutations.
#' @author Lucy Small, \email{lucy.small@@ucdconnect.ie}
#' @param a,b integer vectors
#' @param k integer
#' @return Returns Spearman's footrule between  two partial permutations,
#' where only the first \code{k} items have been ranked. No ties are permitted.

#' @examples
#' a = c(3,1,2,5,4)
#' b = 1:5
#' k=3
#' spearfootP(a,b,k)
#' @export

spearfootP = function(a,b,k){

  n = length(a)
  v = 1:n

  A = v[a<=k &b<=k]
  B = v[a<=k &b>k]
  D = v[a>k &b<=k]
  h=length(B)

  return(h*(2*n+1-h)+sum(abs(a[A]-b[A]))-sum(a[B])-sum(b[D]))

}

#=====================================================================
# Cayley's distance

ncy3 = function(x,k){

  n = length(x)
  t = matrix(0,n,n)

  m1 = 1:n
  m2 = x
  r = 1

  for(i in 1:n){
    if(m1[i] %in% t == FALSE){

      c = 3
      t[r,1] = m1[i]
      t[r,2] = m2[i]
      col = i
      a = match(m2[col],m1)

      while(m2[a] %in% t[r,] == FALSE){

        t[r,c] = m2[a]
        c=c+1
        col = a
        a = match(m2[col],m1)
      }
      r = r+1

    }
  }

  v = c((k+1):n)

  t=t[(rowSums(matrix(t %in% v,nrow(t)))==0) ,]

  nz =   which( t[,1]==0, arr.ind=TRUE)[1]
  return(nz-1)
  #return(t)
}



#====================================================
# code for Cayley's distance from the package 'Rankcluster'
# it's quicker than the code I've written above
#' @rdname cayley
#' @title Cayley's Distance
#' @description \code{cayley} returns Cayley's distance between two full permutations.
#' @author Lucy Small, \email{lucy.small@@ucdconnect.ie}
#'
#' @param x,y integer vectors
#' @return Returns Cayley's distance between \code{x} and \code{y}

#' @examples
#' x = c(3,1,2,5,4)
#' y = 1:5
#' cayley(x,y)


#' @export
cayley = function(x,y)
{

  d=0
  for (i in 1:(length(x) - 1)) {
    if (y[i]!=x[i]) {
      d=d+1
      y[which(y==x[i])] = y[i]
      y[i]=x[i]
    }
  }
  return(d)
}


#' @rdname cayleyP
#' @title Cayley's Distance for Partial permutations
#' @description Computes Cayley's distance between two partial permutations.

#' @author Lucy Small, \email{lucy.small@@ucdconnect.ie}
#' @param a,b integer vectors
#' @param k integer
#' @return Returns Cayley's distance between  two partial permutations,
#' where only the first \code{k} items have been ranked. No ties are permitted.

#' @examples
#' a = c(3,1,2,5,4)
#' b = 1:5
#' k=3
#' cayleyP(a,b,k)
#' @export

cayleyP = function(a,b,k){
  if(k == length(a)){
  return(cayley(a,b))}else if(all(inv(a)[1:k]==inv(b)[1:k])){
  return(0)}else{
  c = k - ncy3(pmult(inv(b),a),k)
    return(c)}
}


#===============================================
# Spearman's rho
#' @rdname spearP
#' @title Spearman's Rho  for Partial permutations
#' @description Computes Spearman's rho  between two partial permutations.
#' @author Lucy Small, \email{lucy.small@@ucdconnect.ie}
#' @param a,b integer vectors
#' @param k integer
#' @return Returns Spearman's rho  between two permutations x and y, where we care about the first k ranked items only

#' @examples
#' a = c(3,1,2,5,4)
#' b = 1:5
#' k=3
#' spearP(a,b,k)
#' @export


spearP = function(a,b,k){
  n = length(a)
  v = 1:n

  A = v[a<=k &b<=k]
  B = v[a<=k &b>k]
  D = v[a>k &b<=k]

  h = length(B)

  p = sort(a[B])
  s = sort(b[D])
  l = length(p)

  return((sum((a[A]-b[A])^2)+h*h*(n-k-h)+
            max(sum((n+1-1:l-p)^2)+sum((k+1:l-s)^2),sum((k-1:l-p)^2)+
               sum((n+1-1:l-s)^2)))^0.5)
}


#==============================================
# Ulam's distance
#' @rdname ulamP
#' @title Ulam's Distance for Partial permutations
#' @description Computes Ulam's distance  between two partial permutations.
#' @author Lucy Small, \email{lucy.small@@ucdconnect.ie}
#' @param a,b integer vectors
#' @param k integer
#' @return Returns Ulam's distance between  two partial permutations,
#' where only the first \code{k} items have been ranked. No ties are permitted.

#' @examples
#' a = c(3,1,2,5,4)
#' b = 1:5
#' k=3
#' ulamP(a,b,k)
#' @export

ulamP=function(a,b,k){
  if(k == length(a)){
    return(ulam(a,b))
  }else{
    n = length(a)
    alph = vector()
    beta = rep(NA,n)

    pm = pmult(a,inv(b))
    ipm = inv(pm)


    n = length(a)
    v = 1:n

    B = v[a<=k &b>k]
    h=length(B)
    # part (1)

    for(i in 1:n){
      if(i <= k && pm[i] <=k){
        alph[pm[i]] =  i
        beta[i] = pm[i]
      }
    }

    # part (2)

    j = sort((1:k)[1:k %in% pm[(k+1):n]])
    jp = sort(((k+1):n)[pm[1:k] >k] )



    for(i in 1:h){
      beta[jp[i]] = n+1-i
      alph[j[i]] = n+1-i
    }

    # part (3)

    alph = c(alph,c(1:n)[!(1:n %in% alph)])

    be = 1:n*(is.na(beta)==T)
    be = be[be!=0]
    beta[be] = c(1:n)[!(1:n %in% beta)]

    return(n-min(lis(alph),lis(beta)))
  }
}


# metrics from pg 36(43) of metric methods book
# items with equal ranking are now allowed, but no. in each group must
# be the same



#======================================================
# Hamming distance
# Hamming distance
#' @rdname hamE
#' @title Hamming Distance with Ties
#' @description Computes the Hamming distance between two permutations, where
#'  items with equal ranking are now permitted. The number of items placed in the
#'  ith category must be the same.


#' @author Lucy Small, \email{lucy.small@@ucdconnect.ie}
#' @param x,y integer vectors
#' @return Returns the Hamming distance between the two permutations.
#' @examples
#' a = c(3,1,2,2,3)
#' b = c(1,2,2,3,3)
#' hamE(a,b)
#' @export

hamE = function(x,y){
  length(x)-length(which((x-y)==0))
}


#======================================================
# Kendall's tau
#' @rdname kendE
#' @title Kendall's Tau for Tankings with Ties
#' @description Computes Kendall's tau between two permutations, where
#'  items with equal ranking are now permitted. The number of items placed in the
#'  ith category must be the same.


#' @author Lucy Small, \email{lucy.small@@ucdconnect.ie}
#' @param x,y integer vectors
#' @return Returns Kendall's tau between the two permutations.

#' @examples
#' a = c(3,1,2,2,3)
#' b = c(1,2,2,3,3)
#' kendE(a,b)
#' @export


kendE = function(x,y){
  s=0
  r=max(x)
  l = length(x)
  n = nig3(x,y)

  for(i in 1:(r-1)){
    for(ip in (i+1):r){
      for(jp in 1:r){
        for(j in jp:r){
          s=s+(n[i,j]*n[ip,jp])
        }
      }
    }
  }
  return(s)
}

#======================================================
# Spearman's footrule
#' @rdname spearfootE
#' @title Spearman's Footrule with Ties
#' @description Computes Spearman's footrule between two permutations, where
#'  items with equal ranking are now permitted. The number of items placed in the
#'  ith category must be the same.


#' @author Lucy Small, \email{lucy.small@@ucdconnect.ie}
#' @param x,y integer vectors
#' @return Returns Spearman's footrule between the two permutations.
#' @examples
#' a = c(3,1,2,2,3)
#' b = c(1,2,2,3,3)
#' spearfootE(a,b)
#' @export


spearfootE = function(x,y){
  r=max(x)
  N = as.vector(table(x))
  t1=0
  t2=0
  n=data.frame(r,r)

  n = nig3(x,y)


  for(i in 1:r){
    for(j in 1:r){
      if(n[i,j]!=0){

        if(i>1){a=sum(N[1:(i-1)])
                d2 = sum(n[1:(i-1),j])
        }else{a=0
              d2=0}

        if(j>1){
          b=sum(N[1:(j-1)])
          c =  sum(n[i, 1: (j-1)])

        }else{
          b=0
          c=0
        }

        if(i<r){d=sum(n[(i+1):r,j])}else{d=0}
        if(j<r){c2= sum(n[i, (j+1):r])}else{c2=0}

        t1 = t1 + n[i,j]*abs(a-b+ c- d)
        t2 = t2 + n[i,j]*abs(a- b +  c2-d2)

    }
  }
}

  return(max(t1,t2))

}


#======================================================
# Spearman's rho

#' @rdname spearE
#' @title Spearman's Rho for permutations with Ties
#' @description Computes Spearman's rho between two permutations, where
#'  items with equal ranking are now permitted. The number of items placed in the
#'  ith category must be the same.


#' @author Lucy Small, \email{lucy.small@@ucdconnect.ie}
#' @param x,y integer vectors
#' @return Returns Spearman's rho between the two permutations.
#' @examples
#' a = c(3,1,2,2,3)
#' b = c(1,2,2,3,3)
#' spearE(a,b)
#' @export

spearE = function(x,y){
  r=max(x)
  N = as.vector(table(x))
  t1=0
  t2=0
  n=data.frame(r,r)

  n = nig3(x,y)



  for(i in 1:r){
    for(j in 1:r){
      if(n[i,j]!=0){

        if(i>1){a=sum(N[1:(i-1)])
                d2 = sum(n[1:(i-1),j])
        }else{a=0
              d2=0}
        if(j>1){
          b=sum(N[1:(j-1)])
          c =  sum(n[i, 1: (j-1)])

        }else{
          b=0
          c=0
        }

        if(i<r){d=sum(n[(i+1):r,j])}else{d=0}
        if(j<r){c2= sum(n[i, (j+1):r])}else{c2=0}

        t1 = t1 + n[i,j]*(a-b+ c- d)^2
        t2 = t2 + n[i,j]*(a- b +  c2-d2)^2
      }
    }
  }

  return(max(t1,t2)^0.5)

}

#======================================================
# Cayley's distance, only works for size of groups = 2 or 3
#' @rdname cayleyE
#' @title Cayley's Distance with Ties
#' @description Computes Cayley's distance between two permutations, where
#'  items with equal ranking are now permitted. The number of items placed in the
#'  ith category must be the same. The number of groups must be 2 or 3, as Cayley's distance is undefined in other cases.
#' @author Lucy Small, \email{lucy.small@@ucdconnect.ie}
#' @param x,y integer vectors
#' @return Returns Cayley's distance between the two permutations.
#' @export
#' @examples
#' a = c(3,1,2,2,3)
#' b = c(1,2,2,3,3)
#' cayleyE(a,b)
cayleyE = function(x,y){
  r=max(x)
  n=data.frame(r,r)

  l = length(x)
  s=0
  s2=0

  n = nig3(x,y)

  if(r == 2){return(n[1,2])}

  if(r == 3){

    for(i in 1:3){
      p1 = s2+n[i,i]
      for(j in i:3){
        s= s + min(n[i,j],n[j,i])
      }
    }


    return(l -(s2+s+abs(n[1,2]-n[2,1])))
  }

}

#===================================================
#ulam - didn;t do this as the code for the more general case is
# v similar



#

# metrics from pg 48 (55) of metric methods book
# items with equal ranking are now allowed, no. in each group
# can be different
#======================================================
# Kendall's tau
#' @rdname kendG
#' @title Kendall's Tau for any Number of Ties
#' @description Computes Kendall's tau between two permutations, where any number of items
#' with equal ranking are now permitted in each ranking.
#' @author Lucy Small, \email{lucy.small@@ucdconnect.ie}
#' @param x,y integer vectors
#' @return Returns Kendall's tau between the two permutations.
#' @examples
#' a = c(3,1,2,2,3)
#' b = c(1,2,3,4,4)
#' kendG(a,b)
#' @export

kendG = function(x,y){
  s=0
  s2=0
  r=max(x)
  rp = max(y)
  n = nig3(x,y)

  for(i in 1:(r-1)){
    for(ip in (i+1):r){
      for(jp in 1:rp){
        for(j in jp:rp){
          s=s+(n[i,j]*n[ip,jp])
        }
      }
    }
  }

  for(i in 1:r){
    for(ip in i:r){
      for(jp in 1:(rp-1)){
        for(j in (jp+1):rp){
          s2=s2+(n[i,j]*n[ip,jp])
        }
      }
    }
  }
  return(max(s,s2))
}

#======================================================
# Spearman's footrule
#' @rdname spearfootG
#' @title Spearman's Footrule for any Number of Ties
#' @description Computes Spearman's footrule between two permutations, where any number of items with equal permutations
#' are now permitted in each ranking.
#' @author Lucy Small, \email{lucy.small@@ucdconnect.ie}
#' @param x,y integer vectors
#' @return Returns Spearman's footrule between the two permutations.
#' @examples
#' a = c(3,1,2,2,3)
#' b = c(1,2,3,4,4)
#' spearfootG(a,b)
#' @export

spearfootG = function(x,y){
  r=max(x)
  rp = max(y)
  N = as.vector(table(x))
  Np = as.vector(table(y))
  t1=0
  t2=0

  n = nig3(x,y)

  for(i in 1:r){
    for(j in 1:rp){
      if(n[i,j] == 0){

      }else{

        if(i>1){a=sum(N[1:(i-1)])
                d2 = sum(n[1:(i-1),j])
        }else{a=0
              d2=0}
        if(j>1){
          b=sum(Np[1:(j-1)])
          c =  sum(n[i, 1: (j-1)])

        }else{
          b=0
          c=0
        }

        if(i<r){d=sum(n[(i+1):r,j])}else{d=0}
        if(j<rp){c2= sum(n[i, (j+1):rp])}else{c2=0}

        t1 = t1 + n[i,j]*abs(a-b+ c- d)
        t2 = t2 + n[i,j]*abs(a- b +  c2-d2)
      }
    }
  }

  return(max(t1,t2))

}



#======================================================
# Spearman's rho
#' @rdname spearG
#' @title Spearman's Rho for any Number of Ties
#' @description Computes Spearman's rho between two permutations x and y, where any number of items with equal rankings
#' are now permitted in each ranking. The number of items
#'  ranked r for the two permutations can vary.
#' @author Lucy Small, \email{lucy.small@@ucdconnect.ie}
#' @param x,y integer vectors
#' @return Returns Spearman's rho between two permutations x and y
#' @examples
#' a = c(3,1,2,2,3)
#' b = c(1,2,3,4,4)
#' spearG(a,b)
#' @export

spearG = function(x,y){

  r=max(x)
  rp = max(y)
  N = as.vector(table(x))
  Np = as.vector(table(y))

  t1=0
  t2=0
  n = nig3(x,y)

  for(i in 1:r){
    for(j in 1:rp){

      if(i>1){a=sum(N[1:(i-1)])
              d2 = sum(n[1:(i-1),j])
      }else{a=0
            d2=0}
      if(j>1){
        b=sum(Np[1:(j-1)])
        c =  sum(n[i, 1: (j-1)])

      }else{
        b=0
        c=0
      }
      n
      rp
      if(i<r){d=sum(n[(i+1):r,j])}else{d=0}
      if(j<rp){c2= sum(n[i, (j+1):rp])}else{c2=0}

      t1 = t1 + n[i,j]*(a-b+ c- d)^2
      t2 = t2 + n[i,j]*(a- b +  c2-d2)^2

    }

  }

  return(max(t1,t2)^0.5)

}


#======================================================
# Hamming distance
#' @rdname hamG
#' @title Hamming Distance for any Number of Ties
#' @description Computes Hamming distance between two permutations, where any number of items with equal rankings
#'  are now permitted in each ranking.
#' @author Lucy Small, \email{lucy.small@@ucdconnect.ie}
#' @param x,y integer vectors
#' @return Returns Hamming distance between the two permutations.
#' @examples
#' a = c(3,1,2,2,3)
#' b = c(1,2,3,4,4)
#' hamG(a,b)
#' @export

hamG = function(x,y){


  R=max(x,y)
  r=max(x)
  rp = max(y)

  N = as.vector(table(x))
  Np = as.vector(table(y))
  Nv= list()
  Nvp=list()

  n = nig3(x,y)

  Nv[[1]] = 1:N[1]
  Nvp[[1]] = 1:Np[1]
  t1 = N[1]
  t2 = Np[1]

  for(i in 2:r){

    Nv[[i]] = (t1+1):(t1+N[i])

    #Nvp[[i]] = (t2+1):(t2+Np[i])
    t1 =t1 +N[i]
    #t2 = t2+Np[i]
  }
  for(i in 2:rp){

    #Nv[[i]] = (t1+1):(t1+N[i])

    #error is here:
    Nvp[[i]] = (t2+1):(t2+Np[i])
    #t1 =t1 +N[i]
    t2 = t2+Np[i]
  }

  v = matrix(0,r,rp)

  for(i in 1:r){
    for(j in 1:rp){
      v[i,j] = length(intersect(Nv[[i]],Nvp[[j]]))
    }
  }
  a=0
  b=0

  for(i in 1:r){
    for(j in 1:rp){
      a = a+max(0,n[i,j]+v[i,j]-Np[j])
      b = b + max(0,n[i,j]+v[i,j]-N[i])

    }
  }

  return(length(x)-min(a,b))

}

#==========
# Ulam's distance
#' @rdname ulamG
#' @title Ulam's distance for any Number of Ties
#' @description Computes Ulam's distance between two permutations, where any number of items with equal rankings
#'  are now permitted in each ranking.
#' @author Lucy Small, \email{lucy.small@@ucdconnect.ie}
#' @param x,y integer vectors
#' @return Returns Ulam's distance between the two permutations.
#' @examples
#' a = c(3,1,2,2,3)
#' b = c(1,2,3,4,4)
#' ulamG(a,b)
#' @export

ulamG = function(x,y){

  l = length(x)
  r=max(x)
  rp = max(y)

  n = nig3(x,y)


  p = perm(rp,r)
  p2 = matrix(0,dim(p)[1],rp)

  for(i in 1:(dim(p)[1])){
    for(j in 1:rp){
      p2[i,j] = n[p[i,j],j]
    }
  }


  p=perm(r,rp)
  p3 = matrix(0,dim(p)[1],r)

  for(i in 1:(dim(p)[1])){
    for(j in 1:r){
      p3[i,j] = n[j,p[i,j]]
    }
  }

  ret1 = max(rowSums(p2))
  ret2 = max(rowSums(p3))

  return(l-min(ret1,ret2))


}

#' @rdname ulamE
#' @title  Ulam's distance with Ties
#' @description Computes  Ulam's distance between two permutations, where
#'  items with equal ranking are now permitted. The number of items placed in the
#'  ith category must be the same.


#' @author Lucy Small, \email{lucy.small@@ucdconnect.ie}
#' @param x,y integer vectors
#' @return Returns  Ulam's distance between the two permutations.
#' @examples
#' a = c(3,1,2,2,3)
#' b = c(1,2,2,3,3)
#' ulamE(a,b)
#' @export
ulamE= function(x,y){

  l = length(x)
  r=max(x)
  rp = max(y)

  n = nig3(x,y)


  p = perm(rp,r)
  p2 = matrix(0,dim(p)[1],rp)

  for(i in 1:(dim(p)[1])){
    for(j in 1:rp){
      p2[i,j] = n[p[i,j],j]
    }
  }


  p=perm(r,rp)
  p3 = matrix(0,dim(p)[1],r)

  for(i in 1:(dim(p)[1])){
    for(j in 1:r){
      p3[i,j] = n[j,p[i,j]]
    }
  }

  ret1 = max(rowSums(p2))
  ret2 = max(rowSums(p3))

  return(l-min(ret1,ret2))


}



#' @title Voting data from 2016 FIFA Best Player of the Year
#'
#' @docType data
#' @keywords datasets
#' @name fifa16
#' @usage data(fifa16)
#' @description This is the voting data with voter covariates from the  FIFA Best Player of the Year 2016 award. There were 
#' 23 candidates on the shortlist and 450 
#' voters. Each voter provides their top-3 choices. The voters were the national captains,
#' manager, and one media representative from each country. 
#' 
#' FIFA, the world football governing body,
#' divides member countries into six continental confederations, which each organise
#' continental national and club competitions. The confederation of each voter is
#' given.
#' @format A data frame with 450 voters (rows) and 30 variables (columns). The first four columns give the voter \code{name}, \code{role} (captain, manager or media), \code{country} of origin and \code{confederation}(AFC, CAF, CONCAF,CONMEBOL or UEFA).
#' Columns \code{vote1}, \code{vote2}, \code{vote3} give the names of the candidates chosen by each voter as their top-3 ranking. 
#' The remaining columns (8:30) give the full partial rankings in permutation form. The votes are arbitrarily filled in after the top-3.
#'
#' @source \url{http://resources.fifa.com/mm/Document/the-best/PlayeroftheYear-Men/02/86/27/05/faward_MenPlayer2016_Neutral.pdf}
NULL

#' @title Voting data from the 2010 UK Labour leadership election.
#' @description There are five candidates 266 rankings, some of which are partial rankings
#' @docType data
#' @keywords datasets
#' @name labour
#' @usage data(labour)
#' @format A data frame with 234 rows and 11 variables
#' 
#' @references \url{https://web.archive.org/web/20110101171158/http://www2.labour.org.uk/leadership-mps-and-meps}
#'
#' @source \url{https://docs.google.com/spreadsheets/d/1e-gx4Km2ywG85kJCf_byJdMZvdP4QkPHGjPKy_meO30/edit?hl=en&hl=en#gid=0}
NULL