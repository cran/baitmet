\name{computeRI}
\alias{computeRI}
\title{Retetion Index computation}
\description{Computes the Retetion Index by elastic curve or reference standards interpolation}
\usage{
computeRI(Experiment, ms.library=mslib, IS.alignid=NULL)
}
\arguments{
  \item{Experiment}{A 'MetaboSet' S4 object containing the experiment data previously created by newExp, and after being deconvolved using decBaitMet.}
  \item{ms.library}{The same MS library used in the deconvolution step.} 
  \item{IS.alignid}{If NULL, the elastic curve modification will be used to compute the RI (no internal standards needed). Otherwise, IS.alignid contains the AlignID number of the compounds to be used as a reference to compute the retention indexes.} 
}
\details{
See BaitMet vignette for more details. To open the vignette, execute the following code in R:
vignette("BaitMetManual", package="baitmet")
}
\value{
The function returns an updated S4 'MetaboSet' class, where the compounds RI has been computed.
}
\references{
[1] Targeting the untargeted: BaitMet, an R package for GC-MS library-driven compound profiling in metabolomics. Xavier Domingo-Almenara, Alexandre Perera-Lluna, Gabriel Vivo-Truyols, Gabriela Venturini, Maria Vinaixa, Jesus Brezmes. (2016) Submitted. 
}
\author{
Xavier Domingo-Almenara. xavier.domingo@urv.cat
}

\seealso{
\code{\link{idList}}
}



