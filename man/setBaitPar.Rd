\name{setBaitPar}
\alias{setBaitPar}
\title{Set BaitMet Software Parameters}
\description{Sets BaitMet Software Parameters.}
\usage{
setBaitPar(ri.error=0.05, min.peak.width, 
min.peak.height=500, noise.threshold=500, 
avoid.processing.mz=c(1:69,73:75,147:149), 
compression.coef=2, analysis.time=0)
}
\arguments{
  \item{ri.error}{The Retention Index error in which compounds are going to be searched. A value of 0.05 corresponds to a 5 percent.}
  \item{min.peak.width}{Minimum compound peak widht (in seconds).}
  \item{min.peak.height}{Minimum compound peak height}
  \item{noise.threshold}{Data above this threshold will be considered as noise}
  \item{avoid.processing.mz}{The masses that do not want to be considered for processing. Typically, in GC-MS those masses are 73,74,75,147,148 and 149, since they are they are ubiquitous mass fragments typically generated from compounds carrying a trimethylsilyl moiety.}
  \item{compression.coef}{Data will be compressed when using OSD algorithm according to this value. This parameter is recomended to bet to 2.}
  \item{analysis.time}{The chromatographic retention time window to process. If 0, all the chromatogram is processed.}
 
 }
\details{
See BaitMet vignette for more details. To open the vignette, execute the following code in R:
vignette("BaitMetManual", package="baitmet")
} 
 
\references{
[1] Targeting the untargeted: BaitMet, an R package for GC-MS library-driven compound profiling in metabolomics. Xavier Domingo-Almenara, Alexandre Perera-Lluna, Gabriel Vivo-Truyols, Gabriela Venturini, Maria Vinaixa, Jesus Brezmes. (2016) Submitted. 
}
\author{
Xavier Domingo-Almenara. xavier.domingo@urv.cat
}

\seealso{
\code{\link{newExp}}, \code{\link{decBaitMet}}, \code{\link{subSetLib}}, \code{\link{setChrmMethod}}
}



