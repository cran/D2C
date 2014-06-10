#' @import RBGL gRbase 


#' @docType methods
setGeneric("compute", function(object,...) {standardGeneric("compute")})


#' An S4 class to represent a directed acyclic graph and its functional dependency.
#' @examples
#' require(D2C)
#' require(gRbase)
#' require(RBGL)
#' H_Rn <- function(n)
#' {
#'   a = runif(n+1,min = -1,max = 1)
#'   f <- function(x)
#'   {
#'     X  = x^(0:n)
#'     return ( sum(X * a))
#'   }
#'   return(Vectorize(f))
#' }
#' H = function() return(H_Rn(1))
#' #H is a function that return a linear function example : 
#' 
#' a = H() ; b= H() ; 
#' print(a(-10:10)) ;  
#' print(b(-10:10)) ;
#' #
#' 
#' sdn=0.5
#' sigma=function(x) return(rnorm(n = 1,sd = sdn))
#' DAG = new("D2C.network",network=random_dag(1:50,maxpar = 5,wgt=0.8),H=H,sdn=sdn,sigma=sigma)
#' X = compute(DAG,N=150)
#' 
#' 
#' #DAG is a network that contains the graph and the functional dependency
#' plot(DAG@@network)
#' #DAG contains node and edges attribute.
#' #nodes attributes : bias the bias of a node , sigma the noise 
#' #edges attributes :  H the dependency function  
setClass("D2C.network",  slots = list(network = "graph"))



#' create a D2C.network
#'
#' \code{new("D2C.network",network=random_dag(1:50,maxpar = 5,wgt=0.8),H=H,sdn=sdn,sigma=sigma)} returns D2C.network object
#' @param network a DAG network 
#' @param sdn the standart deviation of the noise 
#' @param sigma a noise function
#' @param H the dependency function 
#' @param .Object a D2C.network
#' @return N*nNodes matrix

#' @examples
#' require(D2C)
#' require(gRbase)
#' require(RBGL)
#' H=function() function(x) return(x)
#' sdn=0.5
#' sigma=function(x) return(rnorm(n = 1,sd = sdn))
#' DAG = new("D2C.network",network=random_dag(1:50,maxpar = 5,wgt=0.8),H=H,sdn=sdn,sigma=sigma)
#' X = compute(DAG,N=200)
#' #H is a function which return a function. 
#' @references Gianluca Bontempi, Maxime Flauder (2014) From dependency to causality: a machine learning approach. Under submission
#' @export
setMethod("initialize", signature="D2C.network",function(.Object, network,sdn=0.5,sigma=function(x) return(rnorm(n = 1,sd = sdn)),H=function() function(x) return(x) )
{
  DAG = network
  if(!is.DAG(DAG))
  {
    stop("it is not a DAG")
  }
  else
  {
    nodeDataDefaults(DAG,"bias") <-0
    nodeDataDefaults(DAG,"sigma") <-sigma
    edgeDataDefaults(DAG,"H") <- function(x) return(x)
    for( edge in edgeList(DAG))
    {
      edgeData(DAG, from=edge[1], to=edge[2], attr="weight") <- rnorm(1)
      edgeData(DAG, from=edge[1], to=edge[2],attr="H") <- H()
      
    }     
  }   	 
  .Object@network <- DAG
  
  return(.Object)
})





#' compute N samples according to the network distribution
#' @param N  numeric. the number of samples generated according to the network
#' @param object a D2C.network object
#' @return a N*nNodes matrix
#' @examples
#' require(RBGL)
#' require(gRbase)
#' N=150
#' DAG = new("D2C.network",network=random_dag(1:50,maxpar = 5,wgt=0.8))
#' X = compute(DAG,N=N)
#' plot(DAG@@network)
#' @export
setMethod("compute", signature="D2C.network", function(object,N=50)
{
  if(!is.numeric(N))
  {
    stop("N is not numeric")
  }
  DAG = object@network
  nNodes <- numNodes(DAG)
  topologicalOrder <-tsort(DAG)
  
  D <- matrix(NA,nrow=N,ncol=nNodes)
  colnames(D) <- topologicalOrder  
  
  
  for (i in topologicalOrder){   
    bias = nodeData(DAG,n=i,attr="bias")[[1]] 
    sigma = nodeData(DAG,n=i,attr="sigma")[[1]] 
    inEdg <-  inEdges(node=i,object=DAG)[[1]]
    if (length(inEdg)==0){
      D[,i]<-bias + replicate(N,sigma())
    }
    else
    {
      D[,i]<-bias
      for(j in  inEdg)
      {
        inputWeight = edgeData(self=DAG,from=j,to=i,attr="weight")[[1]]
        H = edgeData(self=DAG,from=j,to=i,attr="H")[[1]]
        D[,i]<- D[,i] + H(D[,j]) *  inputWeight
        D[,i] <- scale(D[,i]) + replicate(N,sigma())
      }
    }    
  }
  return(D)
})








