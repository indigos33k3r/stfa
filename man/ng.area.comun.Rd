\name{ng.area.comun}
\alias{ng.area.comun}

\title{
Common area estimation using a procedure based on negentropy to remove outliers.
}
\description{
This function provide a tool to obtain a estimation of SFA-parameters using the area common method. Furthermore it uses a functionality to remove outliers. It lets us use the Normal-HalfNormal and Normal-Exponencial models which are well-known in this field of research.
}
\usage{
ng.area.comun(x, y)
}
\arguments{
  \item{x}{
  The vector or matrix of inputs in the production process.
}
  \item{y}{
 The vector of the ouputs. It is our reponse variable. 
}
}

\value{

\item{landa_ng }{Value of the estimation using the remove-outliers procedure}
\item{landa_ac}{Value of the estimation using the common area method without remove outliers}

}



