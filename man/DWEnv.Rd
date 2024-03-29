% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DWEnv.R
\docType{data}
\name{DWEnv}
\alias{DWEnv}
\title{
Environment containing precomputed objects.
}
\description{
Package environment for commonly used objects in the LS2W methodology.
}
%\usage{ls(DWEnv)}
%\format{
%  The format is:
%<environment: 0x2014238> 
%}
\details{
The \code{DWEnv} environment is used to store precomputed objects for the LS2W methodology.  Certain computations can be fairly long or intensive, and thus it is useful to store precomputed objects for reuse. 
These objects are namely autocorrelation wavelets (output from \code{PsiJ} and \code{PhiJ}), and bias correction matrices produced by \code{D2Amat}. 
}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
%\examples{
%data(DWEnv)
%## maybe str(DWEnv) ; plot(DWEnv) ...
%}
\keyword{data}
