% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rpack2.R
\name{kend}
\alias{kend}
\title{Kendall's Tau}
\usage{
kend(a, b)
}
\arguments{
\item{a, b}{integer vectors}
}
\value{
Returns Kendall's tau between two full rankings.
}
\description{
Computes Kendall's tau between two full rankings.
}
\examples{
a = c(3,1,2,5,4)
b = c(1,2,3,4,5)
kend(a,b)
}
\author{
Lucy Small, \email{lucy.small@ucdconnect.ie}
}
