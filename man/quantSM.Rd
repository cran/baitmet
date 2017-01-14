\name{quantSM}
\alias{quantSM}
\title{Quantification of selective masses}
\description{Quantification of the compounds selective masses for quantitative validation.}
\usage{
quantSM(Experiment, ms.library, AlignID=NULL, pre.process=TRUE, fit.gaussian = FALSE)
}
\arguments{
  \item{Experiment}{A 'MetaboSet' S4 object containing the experiment data previously created by newExp, and after being deconvolved using decBaitMet.}
  \item{ms.library}{The same MS library used in the deconvolution step.} 
  \item{AlignID}{If NULL, all the compounds of the experiment are re-quantified by their masses. Otherwise, AlignID (vector) contains the AlignID number of the compounds to be used as a reference to compute the retention indexes.} 
  \item{pre.process}{If TRUE (highly recommended), a pre-processing is conducted prior to quantification of masses.} 
  \item{fit.gaussian}{If TRUE, a gaussian is fitted to the selected masss.} 

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
\code{\link{mzList}}, \code{\link{plotMZ}}
}



