\name{CTWOS}
\alias{CTWOS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Scalable calibration transfer without standards
}
\description{
Scalable calibration transfer without standards via dynamic time warping for near-infrared spectroscopy
}
\usage{
CTWOS (MX,Y,SX,ytest,width,order,deriv,span,A0)
}
\arguments{
  \item{MX}{master spectra,n x p matrix(Sample matrix)}
  \item{Y}{measured property of master spectra}
  \item{SX}{slave spectra,m x p matrix(Sample matrix)}
  \item{ytest}{measured property of slave spectra}
  \item{width}{the optional variables specifying the number of points in filter}
  \item{order}{the order of the polynomial}
  \item{deriv}{the derivative}
  \item{span}{the width of the moving window, integer. We recommend using a dilation to form a penalty for use in VPdtw. A dilation is a moving local maximum over a window of specific fixed width specified by span.}
  \item{A0}{The maximal number of latent variables for cross-validation}
}
\details{
With the good transfer ability and the standard-free advantage, CTWOS has the potential to be a widely used method for solving the calibration transfer problem in near-infrared spectroscopy.
}
\value{
Return a list of result
  \item{RMSEP}{Root-Mean-Square Error of Prediction}
  \item{Q2}{A parameter for evaluating model prediction performance, the closer Q2 is to 1, the better the model}
  \item{ypred}{The value obtained by model prediction}
}
\references{
D. Clifford, G. Stone, I. Montoliu, S. Rezzi, F.-P. Martin, P. Guy, S. Bruce, S. Kochhar, Alignment Using Variable Penalty Dynamic Time Warping, Analytical Chemistry, 81 (2009) 1000-1007.
}
\examples{
data(wheat)
MX = wheat$A.Test
SX = wheat$vliadA[c(1:55,57:57,60:60,64:121,122:135,136:150),]
ytest=wheat$y.Valid[c(1:55,57:57,60:60,64:121,122:135,136:150)]
Y=wheat$Y.CAL
width=21
order=2
deriv=2
span=5
A0=20
result <- CTWOS (MX,Y,SX,ytest,width,order,deriv,span,A0)
RMSEP=result$RMSEP
}
\keyword{ CTWOS }
\keyword{ calibration transfer }
