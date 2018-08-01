\name{adkernel}
\alias{adkernel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Kernel estimation using an adaptive kernel.
}
\description{
Common area estimation using an adaptive kernel.
}
\usage{
adkernel(x, kernel = dnorm, bw = bw.nrd0, adjust = 1, n = 500, from, to, cut = 3, na.rm = TRUE, alpha = 0.5)
}

\arguments{
  \item{x}{
A vector of values.
}
  \item{kernel}{
The kernel function.
}
  \item{bw}{
Band-width method.
}
  \item{adjust}{
Parameter to adjust tha variability of the bw.
}
  \item{n}{
Number of cuts over the values.
}
  \item{from}{
Minimun value to cut.
}
  \item{to}{
Maximum value to cut
}
  \item{cut}{
Size of the cuts.
}
  \item{na.rm}{
Treatment of the NA's.
}
  \item{alpha}{
An parameter of adjust teh adaptability.
}
}

\value{
Return a "density" class objetc.
}



