#' @import RBGL gRbase 


#' @docType methods
setGeneric("compute", function(object,...) {standardGeneric("compute")})



#' An S4 class to store DAG.network
#' @param network : object of class "graph"

setClass("DAG.network",  slots = list(network = "graph"))



#' creation of a DAG.network
#' @param .Object : DAG.network object
#' @param sdn : standard deviation of aditive noise. 
#' @param sigma : function returning the additive noise
#' @param H : function describing the type of the dependency. 
#' @param network : object of class "graph"
setMethod("initialize", signature="DAG.network",
          function(.Object, network,
                   sdn=0.5,
                   sigma=function(x) return(rnorm(n = 1,sd = sdn)),
                   H=function(x) return(x) ){
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
                edgeData(DAG, from=edge[1], to=edge[2], attr="weight") <- runif(1,0.5,1)*sample(c(-1,1),1)
                edgeData(DAG, from=edge[1], to=edge[2],attr="H") <- H()
                
              }     
            }   	 
            .Object@network <- DAG
            
            return(.Object)
          })





#' compute N samples according to the network distribution
#' @param N  numeric. the number of samples generated according to the network
#' @param object a DAG.network object
#' @return a N*nNodes matrix
#' @export
setMethod("compute", signature="DAG.network", function(object,N=50)
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
        ## it computes the linear combination of the inputs
      {
        inputWeight = edgeData(self=DAG,from=j,to=i,attr="weight")[[1]]
        H = edgeData(self=DAG,from=j,to=i,attr="H")[[1]]
        D[,i]<- D[,i] + H(D[,j]) *  inputWeight
        D[,i] <- scale(D[,i]) + replicate(N,sigma())
      }
    }    
  }
  col.numeric<-as(colnames(D),"numeric")
  D<-D[,topologicalOrder[order(col.numeric)]]
  
  return(D)
})



#' An S4 class to store a list of DAGs and associated observations
#' @param list.DAGs : list of stored DAGs
#' @param list.observationsDAGs : list of observed datasets, each sampled from the corresponding member of list.DAGs 
#' @param NDAG  : number of DAGs. 
#' @param functionType : type of the dependency. It is of class "character" and is one of  ("linear", "quadratic","sigmoid")
#' @param seed : random seed
setClass("simulatedDAG",
         slots = list(list.DAGs="list",list.observationsDAGs="list",
                      NDAG="numeric", functionType="character",seed="numeric"))




#' creation of a "simulatedDAG" containing a list of DAGs and associated observations 
#' @param .Object : simulatedDAG object
#' @param NDAG : number of DAGs to be created and simulated
#' @param noNodes  : number of Nodes of the DAGs. If it is a two-valued vector , the value of Nodes is randomly sampled in the interval 
#' @param N  : number of sampled observations for each DAG. If it is a two-valued vector [a,b], the value of N is randomly sampled in the interval [a,b]
#' @param sdn : standard deviation of aditive noise. If it is a two-valued vector, the value of N is randomly sampled in the interval
#' @param seed : random seed
#' @param verbose : if TRUE it prints out the state of progress
#' @param functionType : type of the dependency. It is of class "character" and is one of  ("linear", "quadratic","sigmoid")
#' @references Gianluca Bontempi, Maxime Flauder (2014) From dependency to causality: a machine learning approach. Under submission
#' @examples
#' require(RBGL)
#' require(gRbase)
#' descr=new("D2C.descriptor")
#'descr.example<-new("D2C.descriptor",bivariate=FALSE,ns=3,acc=TRUE)
#'trainDAG<-new("simulatedDAG",NDAG=10, N=c(50,100),noNodes=c(15,40),
#'              functionType = "linear", seed=0,sdn=c(0.45,0.75)) 
#' @export
#' 
#' 
setMethod("initialize",
          "simulatedDAG",
          function(.Object, NDAG=10,
                   noNodes=sample(10:20,size=1),functionType="linear",
                   verbose=TRUE,N=sample(100:500,size=1),
                   seed=1234,sdn=0.5) {
            
            #generate a training set
            # NDAG the number of network to use
            #functionType example : "R1" "R2" "sigmoid1"
            set.seed(seed)
            
            
            
            .Object@NDAG=NDAG
            .Object@functionType=functionType
            .Object@seed=seed
            X=NULL
            Y=NULL
            list.DAGs=NULL
            list.observationsDAGs=NULL
            for(i in 1:NDAG){
              N.i<-N
              if (length(N)>1)
                N.i<-sample(N[1]:N[2],1)
              
              noNodes.i<-noNodes
              if (length(noNodes)>1)
                noNodes.i<-sample(noNodes[1]:noNodes[2],1)
              
              sdn.i<-sdn
              if (length(sdn)>1)
                sdn.i<-runif(1,sdn[1],sdn[2])
              
              
              
              
              wgt = runif(n = 1,min = 0.65,max = 0.85)
              V=1:noNodes.i
              
              maxpar = sample(2:round(noNodes.i/2),size=1)
              
              
              if(functionType=="linear"){
                H = function() return(H_Rn(1))
                
              }else if(functionType=="quadratic"){
                H = function() return(H_Rn(2))
                
              }else if(functionType=="sigmoid"){
                H = function() return(H_sigmoid(1))
                
              }
              netwDAG<-random_dag(V,maxpar = maxpar,wgt)
              cnt<-2
              while (sum(unlist(lapply(edges(netwDAG),length)))<(length(V)) & cnt<10){
                netwDAG<-random_dag(V,maxpar = maxpar,wgt)
                cnt<-cnt+1
                
              }
              
              DAG = new("DAG.network",
                        network=netwDAG,H=H,sdn=sdn.i)
              
              
              observationsDAG = compute(DAG,N=N.i)
              
             
              
              list.DAGs<-c(list.DAGs,list(netwDAG))
              list.observationsDAGs<-c(list.observationsDAGs,list(observationsDAG))
              if (verbose)
                cat("simulatedDAG: DAG number:",i,"generated: #nodes=", length(V), 
                    "# edges=",sum(unlist(lapply(edges(netwDAG),length))), "# samples=", N.i, "\n")
              
            }
            
            
            .Object@list.DAGs=list.DAGs
            .Object@list.observationsDAGs=list.observationsDAGs
            
            
            .Object
          }
)










