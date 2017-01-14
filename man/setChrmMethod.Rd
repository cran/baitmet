\name{setChrmMethod}
\alias{setChrmMethod}
\title{Sets the chromatographic method to be used}
\description{Sets the relation between the retention times and indexes for a chromatographic method to be used}
\usage{
setChrmMethod(method=c("alk.var5","alk.mdn35",
"fame.var5","fame.mdn35"), rt, ri, name="ChromMethod #1")
}
\arguments{
  \item{method}{The chromatographic method type to be used, currently, only ALK-VAR5 (alk.var5) type is available. However, this method is compatible with both Alkalanes and FAMEs, but only for chromatographic methods using VAR5-type columns.}
  \item{rt}{The retention times (RT) of the reference standards.}
  \item{ri}{The retention indexes (RI) of the reference standards.} 
  \item{name}{The name of the chromatographic method, for further information.} 
}
\details{
This function allows initializing a chromatographic method to be used with BaitMet. Please, see BaitMet vignette for more details. To open the vignette, execute the following code in R:
vignette("BaitMetManual", package="baitmet")
}
\references{
[1] Targeting the untargeted: BaitMet, an R package for GC-MS library-driven compound profiling in metabolomics. Xavier Domingo-Almenara, Alexandre Perera-Lluna, Gabriel Vivo-Truyols, Gabriela Venturini, Maria Vinaixa, Jesus Brezmes. (2016) Submitted. 
}
\author{
Xavier Domingo-Almenara. xavier.domingo@urv.cat
}

\seealso{
\code{\link{decBaitMet}}
}
\examples{

# To set a new chromatographic method, simply 
# indicate the retention time and indexes of the
# compounds to be used as a reference for a given
# chromatographic method.

ChrM.Plasma <- setChrmMethod(method="alk.var5", 
rt=c(8.39, 10.759, 13.271, 15.604, 17.751, 19.685, 
21.471, 23.126, 24.645), ri=c(1225.27, 1326.95, 
1526.31, 1729.57, 1928.52, 2131.7, 2342.06, 
2548.56, 2739.86))

}


