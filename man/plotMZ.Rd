\name{plotMZ}
\alias{plotMZ}
\title{Plotting compounds selected masses}
\description{Plots the chromatophic profiles of the masses (selective) for each compound.}
\usage{
plotMZ(Experiment, AlignId, mz.ind=1, per.class=T, aligned=FALSE, xlim=NULL)
}
\arguments{
 \item{Experiment}{A 'MetaboSet' S4 object containing the experiment after the masses have been quantified by quantSM() function. }
   \item{AlignId}{the Id identificator for the compound to be shown.}
   \item{mz.ind}{Integer. the -th selective mass, the 1 (st) selective mass, the 2 (nd) selective mass, etc..}
   \item{per.class}{logical. if TRUE the profiles are shown one color per class, if FALSE one color per sample.}
   \item{aligned}{logical. if TRUE the profiles (masses) are shown aligned for a better visual comparison, if FALSE they are shown as they are.}
   \item{xlim}{x axsis (retention time) limits (see \code{\link{plot.default}}).}
}
\author{
Xavier Domingo-Almenara. xavier.domingo@urv.cat
}

\seealso{
\code{\link{mzList}}, \code{\link{quantSM}}
}
