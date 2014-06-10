#' @import RBGL gRbase randomForest




#' Dataset example
#'@title stored D2C object
#'@description D2C object for testing D2C functionalities
#' @name example
#' @docType data
#' @keywords data
#' #' @examples
#' require(RBGL)
#' require(gRbase)
#' data(example)
#' print(example@@mod)
#' ## Random Forest
#' print(dim(example@@X))
#' ## dimension of the training set
NULL


#' An S4 class to store the descriptor parameters
setClass("D2C.descriptor",
         slots = list(lin="logical", acc="logical",
                      struct="logical",pq="numeric",
                      bivariate="logical",ns="numeric"))

#' creation of a D2C.descriptor 
#' @param .Object : the D2C.descriptor object
#' @param lin :	TRUE OR FALSE: if TRUE it uses a linear model to assess a dependency, otherwise a local learning algorithm
#' @param acc : TRUE OR FALSE: if TRUE it uses the accuracy of the regression as a descriptor
#'  @param struct	: TRUE or FALSE to use the ranking in the markov blanket as a descriptor
#'  @param pq :a vector of quantiles used to compute the descriptors
#'  @param bivariate :TRUE OR FALSE: if TRUE it includes also the descriptors of the bivariate dependence
#'  @param ns : size of the Markov Blanket returned by the mIMR algorithm
#' @references Gianluca Bontempi, Maxime Flauder (2014) From dependency to causality: a machine learning approach. Under submission
#' @examples
#' require(RBGL)
#' require(gRbase)
#'descr.example<-new("D2C.descriptor",bivariate=FALSE,ns=3,acc=TRUE)
#'trainDAG<-new("simulatedDAG",NDAG=2, N=50,noNodes=10,
#'              functionType = "linear", seed=0,sdn=0.5) 
#'
#' @export
setMethod("initialize",
          "D2C.descriptor",
          function(.Object, lin=TRUE, acc=TRUE,
                   struct=TRUE,pq=c(0.1, 0.25, 0.5, 0.75, 0.9),
                   bivariate=FALSE,ns=4) 
          {
            
            .Object@lin <- lin
            .Object@acc <- acc
            .Object@struct <- struct
            .Object@bivariate <- bivariate
            .Object@pq <- pq
            .Object@ns <- ns
            .Object
          })

setOldClass("randomForest")

#' An S4 class to store the RF model trained on the basis of the descriptors of NDAG DAGs
setClass("D2C",
         slots = list(mod="randomForest", X="matrix",Y="numeric",                    
                      descr="D2C.descriptor",features="numeric",rank="numeric",
                      allEdges="list"
         ))

#' creation of a D2C object which preprocesses the list of DAGs and observations contained in sDAG and fits a  Random Forest classifier
#' @param .Object : the D2C object
#' @param sDAG : simulateDAG object
#' @param descr  : D2C.descriptor object containing the parameters of the descriptor
#' @param max.features  : maximum number of features used by the Random Forest classifier \link[randomForest]{randomForest}. The features are selected by the importance returned by the function \link[randomForest]{importance}.
#' @param ratioEdges  : percentage of existing edges which are added to the training set
#' @param ratioMissingNode  : percentage of existing nodes which are not considered. This is used to emulate latent variables.
#' @param verbose  : if TRUE it prints the state of progress
#' @references Gianluca Bontempi, Maxime Flauder (2014) From dependency to causality: a machine learning approach. Under submission
#' @examples
#' require(RBGL)
#' require(gRbase)
#' descr=new("D2C.descriptor")
#'descr.example<-new("D2C.descriptor",bivariate=FALSE,ns=3,acc=TRUE)
#'trainDAG<-new("simulatedDAG",NDAG=2, N=50,noNodes=10,
#'              functionType = "linear", seed=0,sdn=0.5) 
#' example<-new("D2C",sDAG=trainDAG, descr=descr.example)
#' @export
setMethod("initialize",
          "D2C",
          function(.Object, sDAG, 
                   descr=new("D2C.descriptor"),
                   verbose=TRUE, 
                   ratioMissingNode=0,
                   ratioEdges=1,max.features=20) {
            
            #generate a training set
            # NDAG the number of network to use
            #functionType example : "R1" "R2" "sigmoid1"
            
            .Object@descr=descr
            
            
            X=NULL
            Y=NULL
            allEdges=NULL
            for(i in 1:sDAG@NDAG){
              
              DAG = sDAG@list.DAGs[[i]]
              observationsDAG =sDAG@list.observationsDAGs[[i]]
              
              Nodes = nodes(DAG)  
              keepNode = sort(sample(Nodes,size = ceiling(length(Nodes)*(1-ratioMissingNode)),replace = F))
              
             
              
              DAG2 =subGraph(keepNode, DAG)
              
              
              ##choose wich edge to train / predict and find the right label
              nEdge = length(edgeList(DAG))
              
              edges = matrix(unlist(sample(edgeList(DAG2),size = round(nEdge*ratioEdges),replace = F)),ncol=2,byrow = TRUE)  
              edges = rbind(edges,t(replicate(n =round(nEdge*ratioEdges) ,sample(keepNode,size=2,replace = FALSE))))
              
              nEdges =  NROW(edges)
              labelEdge = numeric(nEdges)
              for(j in 1:nEdges)
              {
                I =edges[j,1] ; 
                J =edges[j,2] ;
                labelEdge[j] = as.numeric(I %in% inEdges(node = J,DAG2)[[1]])
              } 
              
              
              #compute the descriptor for the edges 
              nNodes = length(labelEdge)
              X.out = NULL
              for(j in 1:nNodes)
              {
                I =as(edges[j,1],"numeric") ; 
                J =as(edges[j,2],"numeric") ; 
                
                d<-descriptor(observationsDAG,I,J,lin=descr@lin,acc=descr@acc,
                              struct=descr@struct,bivariate=descr@bivariate,
                              pq=descr@pq,ns=descr@ns)
                
                X.out = rbind(X.out,d)
              }
              
              X=rbind(X,X.out)
              Y=c(Y,labelEdge)
              allEdges<-c(allEdges,list(edges))
              
              
              
              if (verbose)
                cat("D2C:  DAG", i, " processed: # positive examples=",length(which(Y==1)),
                    "# negative examples=",length(which(Y==0)),"\n")
              
              
            }
            
            
            
            
            features<-1:NCOL(X)
            wna<-which(apply(X,2,sd)<0.01)
            if (length(wna)>0)
              features<-setdiff(features,wna)
            
            X<-scale(X[,features])
            .Object@features=features
            .Object@X=X
            .Object@Y=Y
            .Object@allEdges=allEdges
            RF <- randomForest(x =X ,y = factor(Y),importance=TRUE)
            IM<-importance(RF)[,"MeanDecreaseAccuracy"]
            rank<-sort(IM,decr=TRUE,ind=TRUE)$ix[1:min(max.features,NCOL(X))]
            RF <- randomForest(x =X[,rank] ,y = factor(Y))
            .Object@rank=rank
            .Object@mod=RF
            
            .Object
          }
)



#' predict if there is a connection between node i and node j 
#' @param object : a D2C object 
#' @param i :  index of putative cause (\eqn{1 \le i \le n})
#' @param j  : index of putative effect (\eqn{1 \le j \le n})
#' @param data : dataset of observations from the DAG    
#' @return list with  response and prob of the prediction
#' @examples  
#' require(RBGL)
#' require(gRbase)
#' data(example)
#'## load the D2C object
#' testDAG<-new("simulatedDAG",NDAG=1, N=50,noNodes=5,
#'            functionType = "linear", seed=1,sdn=c(0.25,0.5))
#' ## creates a simulatedDAG object for testing
#' plot(testDAG@@list.DAGs[[1]])
#' ## plot the topology of the simulatedDAG 
#' predict(example,1,2, testDAG@@list.observationsDAGs[[1]])
#' ## predict if the edge 1->2 exists
#' predict(example,4,3, testDAG@@list.observationsDAGs[[1]])
#' ## predict if the edge 4->3 exists
#' predict(example,4,1, testDAG@@list.observationsDAGs[[1]])
#' ## predict if the edge 4->1 exists
#' @references Gianluca Bontempi, Maxime Flauder (2014) From dependency to causality: a machine learning approach. Under submission
#' @export
setMethod("predict", signature="D2C",
          function(object,i,j,data)
          {
            out = list() 
            
            
            X_descriptor = descriptor(data,i,j,lin = object@descr@lin, acc = object@descr@acc,ns=object@descr@ns,
                                      struct = object@descr@struct, pq = object@descr@pq, bivariate =object@descr@bivariate)
            
            
           
            X_descriptor=X_descriptor[object@features]
            
            
            X_descriptor=scale(array(X_descriptor,c(1,length(X_descriptor))),attr(object@X,"scaled:center"),attr(object@X,"scaled:scale"))
            
            out[["response"]] = predict(object@mod, X_descriptor[,object@rank], type="response")
            out[["prob"]] = predict(object@mod, X_descriptor[,object@rank], type="prob")
            
            
            return(out)  
          })



