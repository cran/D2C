% Generated by roxygen2 (4.0.2): do not edit by hand
\docType{methods}
\name{initialize,simulatedDAG-method}
\alias{initialize,simulatedDAG-method}
\title{creation of a "simulatedDAG" containing a list of DAGs and associated observations}
\usage{
\S4method{initialize}{simulatedDAG}(.Object, NDAG = 1,
  noNodes = sample(10:20, size = 1), functionType = "linear",
  quantize = FALSE, verbose = TRUE, N = sample(100:500, size = 1),
  seed = 1234, sdn = 0.5, goParallel = FALSE)
}
\arguments{
\item{.Object}{: simulatedDAG object}

\item{NDAG}{: number of DAGs to be created and simulated}

\item{noNodes}{: number of Nodes of the DAGs. If it is a two-valued vector , the value of Nodes is randomly sampled in the interval}

\item{functionType}{: type of the dependency. It is of class "character" and is one of  ("linear", "quadratic","sigmoid")}

\item{quantize}{: if TRUE it discretize the observations into two bins. If it is a two-valued vector [a,b], the value of quantize is randomly sampled in the interval [a,b]}

\item{verbose}{: if TRUE it prints out the state of progress}

\item{N}{: number of sampled observations for each DAG. If it is a two-valued vector [a,b], the value of N is randomly sampled in the interval [a,b]}

\item{seed}{: random seed}

\item{sdn}{: standard deviation of aditive noise. If it is a two-valued vector, the value of N is randomly sampled in the interval}

\item{goParallel}{: if TRUE it uses parallelism}
}
\description{
creation of a "simulatedDAG" containing a list of DAGs and associated observations
}
\examples{
require(RBGL)
require(gRbase)
require(foreach)
descr=new("D2C.descriptor")
descr.example<-new("D2C.descriptor",bivariate=FALSE,ns=3,acc=TRUE)
trainDAG<-new("simulatedDAG",NDAG=10, N=c(50,100),noNodes=c(15,40),
             functionType = "linear", seed=0,sdn=c(0.45,0.75))
}
\references{
Gianluca Bontempi, Maxime Flauder (2014) From dependency to causality: a machine learning approach. Under submission
}

