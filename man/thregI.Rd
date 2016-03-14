\name{thregI}
\alias{thregI}
\title{fit threshold regression for interval-censored data}
\description{
Extend right-censored data of Xiao et al. (2015) into interval-censored data. Fit a threshold regression model for interval-censored data and calculate the prediction of the initial health status value and the drift value of the health process.
 }
\usage{
thregI(formula, data)
}

\arguments{
  \item{formula}{a formula object, which has the response on the left of a ~ operator, and the explaination variables on the right. The response must be a survival object as returned by the Surv(). The | operator on the left-hand side is used for the linear regression function of \eqn{ln{y_0}} and right-hand side is for specifing explaination variables of \eqn{\mu}. A constant \eqn{ln{y_0}} or \eqn{\mu} is allowed, let 0 or 1 as a placeholder on the left or right of the | operator. }

  \item{data}{a data frame must contains two time points (\eqn{L} and \eqn{R}) in which to through the Surv() of R package survival by type = interval2 and interpret the variables named in the formula.}

}
\author{
Man-Hua Chen

Maintainer: Man-Hua Chen < mchen@tku.edu.tw >
}
\references{
Xiao T, Whitmore GA, He X and Lee M-LT. (2015) The R Package threg to Implement Threshold Regression Models. \emph{Journal of Statistical Software} 66, 1141-1156.

Xiao T. (2013) \emph{threg: Threshold Regression}. R~package version 1.0.3, URL \url{http://CRAN.R-project.org/package=threg}.
}


\examples{
#load the data "bcos"
data("bcos", package="thregI")

#transform the "treatment" variable into factor variable f.treatment
bcos$f.treatment=factor(bcos$treatment)

#fit the threshold regression model on the factor variable f.treatment
fit<-thregI(Surv(left, right, type='interval2')~f.treatment|f.treatment, data=bcos)
fit
}

\keyword{interval-censored data}
\keyword{survival analysis}
\keyword{threshold regression}
