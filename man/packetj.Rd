% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/packetj.R
\docType{methods}
\name{packetj}
\alias{packetj}

\title{
Converts an imwd object to a wst2D object at a given level. 
}
\description{
Reorders the informtion at a given level of an \code{imwd} object to be in the form required for a \code{wst2D} object.
}
\usage{
packetj(imwd, level, o)
}

\arguments{
  \item{imwd}{
An \code{imwd} object from a non-decimated wavelet transform, i.e. \code{type} is \code{station}.
}
  \item{level}{
The level we wish to convert. 
}
  \item{o}{
Computes weaving permutation for conversion from imwd to wst2D.
}
}
\details{
This function combines the information from an \code{imwd} object at a given level and reorders it to be in the approriate form for a \code{wst2D} object. This function is required by \code{\link{convertimwd}}.   
}
\value{
Returns a matrix to be put into a \code{wst2D} object.
}

\note{Not intended to be used}

\author{
Matt Nunes
}

\examples{

#Obtain an imwd object
#
testimage <- HaarMontage(256, "diagonal")
testimageIMWD<- imwd(testimage, type="station")
#
#Specify our weaving permutation
#
arrvec <- getarrvec(9, sort=FALSE)
#
#Convert level 6 coefficients into packet ordered object
#
o <- arrvec[,2]
packmat <- packetj(testimageIMWD, 6, o)
#
}

\keyword{classes}

