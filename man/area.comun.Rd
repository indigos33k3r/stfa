\name{area.comun}
\alias{area.comun}

\title{ Common area estimation for the SFA-parameters
}
\description{
It lets use the common area method to obtain a robust estimation of the SFA-parameters. This function include the Normal-HalfNormal and Normal-Exponencial models. In addition several options of kernel estimation has been included.
}
\usage{
area.comun(x, y, model = "Normal-HalfNormal", method = "Kernel")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A vector or matrix with the inputs.}
  \item{y}{
A vector of outputs }
  \item{model}{
  The model that we want to suposed. It can be the "Normal-HalfNormal" or "Normal-Exponencial".
}
  \item{method}{
The method of kernel estimation that we want to use. It can be  "Kernel" or "Adaptative-Kernel".
}
}

\value{
It returns a list with the results.
}


\examples{

set.seed(3)
b0<-100
b1<-5
landa<-1
n<-100
sigmav <- sqrt(sigma2v)
sigma2u <- landa^2 * sigma2v
sigmau <- sqrt(sigma2u)
x <- runif(n, 1, 30)
v <- rnorm(n, 0, sigmav)
u <- abs(rnorm(n, 0, sigmau))
e.real <- v - u
y <- b0 + b1 * x + e.real
asim<-skewness(e.real)

l1<-area.comun(x,y,model="Normal-HalfNormal",method='Adaptative-Kernel')$landa_ac
l2<-area.comun(x,y,model="Normal-HalfNormal",method='Kernel')$landa_ac
cat("Para un valor del parámetro =",landa,"los resultados son:","\n","adkernel=",l1,'\n','kernel=',l2)

}
