#================================
#' @rdname perm
#' @title non decreasing permutations
#' @description generates all non decreasing permutations of 1:b of length k
#' @param k,b integers
#' @return returns all non decreasing permutations of 1:b of length k
#' @keywords internal
#'
perm = function(k,b){

  if( k >= b){
    v = c(1:k)
    L = b
    M = matrix(0,k,L^k)

    c2=1
    del=0

    for(i in 1:k){
      c = 1
      for(j in 1:(L^(i-1))){
        for(q in 1:L){
          for(m in 1:(L^(k-i))){
            M[i,c] = v[q]

            if(i>1){
              if( v[q] < M[i-1,c]){
                del[c2] = c
                c2=c2+1
              }
            }
            c = c + 1
          }
        }
      }
    }

    M = t(M)
    M2 = M[-del,]
  }
  else{

    k1=k
    k=b
    v = c(1:k)
    L = b
    M = matrix(0,k1,L^k1)

    c2=1
    del=0

    for(i in 1:k1){
      c = 1
      for(j in 1:(L^(i-1))){
        for(q in 1:L){
          for(m in 1:(L^(k1-i))){
            M[i,c] = v[q]

            if(i>1){
              if( v[q] < M[i-1,c]){
                del[c2] = c
                c2=c2+1
              }
            }
            c = c + 1
          }
        }
      }
    }

    M = t(M)
    M2 = M[-del,]
    M2 = M2[,1:k1]
    M2 = unique(M2)

  }

  return(M2)
}


#========================================
#' @rdname pmult
#' @title Metrics for Rankings
#' @description Multiplies two Permutations.
#' @author Lucy Small, \email{lucy.small@@ucdconnect.ie}
#'
#' @param x,y permutations
#' @return Returns the composition of two permutations.

#' @examples
#' a = 1:5
#' b = c(3,1,2,5,4)
#' pmult(a,b)
#' @keywords internal
pmult<-function(x,y)
{
  x[y]
}

