#' @import RBGL gRbase randomForest 



setClass("D2C.Predict",  slots = list(data = "matrix" , lin="logical"  ))




#' create a data structure that store a sample of a network  distribution 
#'
#' @param .Object a D2C.predict object 
#' @param data : matrix of samples
#' @param lin : logical
#' @return D2C.Predict object
setMethod("initialize", signature="D2C.Predict",function(.Object, data ,lin=stop("lin must be specified"))
{
  .Object@lin = lin ; 
  .Object@data = data ;
  return(.Object) ; 
})




#' predict if there is a connection between node i and node j 
#' @param object : a D2C.Predict object 
#' @param i : node name
#' @param j  : node name
#' @param networkType "small" , "medium" , "large'
#' @param functionType  "linear" , "quadratic", "sigmoid"   
#' @param RF a randomForest object
#' @return list with  response and prob of the prediction
#' @references Gianluca Bontempi, Maxime Flauder (2014) From dependency to causality: a machine learning approach. Under submission
#' @export
setMethod("predict", signature="D2C.Predict",function(object,i,j,networkType=stop("networkType must be defined")
                                                      ,functionType=stop("stop function type must be defined"),RF =NULL)
{
  out = list() 
  if( is.null(RF))
  {
    if (! functionType %in% c("linear","quadratic","sigmoid"))
    {
      stop("bad functionType name")
    }
    
    if (! networkType %in% c("small","medium","large"))
    {
      stop("bad networkType name")
    }
    RF = RF.package[[paste0(networkType,functionType)]]
  }
 
  
  
  X_descriptor = descriptor(object@data,i,j,lin=object@lin)
  out[["response"]] = predict(RF, X_descriptor, type="response")
  out[["prob"]] = predict(RF, X_descriptor, type="prob")
  return(out)  
})
