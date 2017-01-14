\name{decBaitMet}
\alias{decBaitMet}
\title{Library-driven deconvolution of compounds by BaitMet}
\description{Library-driven deconvolution of GC-MS data by BaitMet}
\usage{
decBaitMet(Experiment, BaitParameters, 
ms.library=mslib, chrom.method, 
samples.to.process=NULL)
}
\arguments{
  \item{Experiment}{A 'MetaboSet' S4 object containing the experiment data previously created by newExp.}
  \item{BaitParameters}{The BaitMet deconvolution parameters object previously created by \code{\link{setBaitPar}}}
  \item{ms.library}{The mass-spectra library to be used for retention time and spectral comparison. By default, the MassBank [2] - Mass Bank of North America (MoNa) database is employed. However, the Golm Metabolome Database [3] is highly recommended.}
  \item{chrom.method}{The chromatographic method previously created by \code{\link{setChrmMethod}}.} 
  \item{samples.to.process}{Vector indicating which samples are to be processed.} 
}
\details{
See BaitMet vignette for more details. To open the vignette, execute the following code in R:
vignette("BaitMetManual", package="baitmet")
}
\value{
The function returns an updated S4 'MetaboSet' class, where the compounds in the library have been searched in the GC-MS samples and deconvolved.
}
\references{
[1] Targeting the untargeted: BaitMet, an R package for GC-MS library-driven compound profiling in metabolomics. Xavier Domingo-Almenara, Alexandre Perera-Lluna, Gabriel Vivo-Truyols, Gabriela Venturini, Maria Vinaixa, Jesus Brezmes. (2016) Submitted. 

[2] MassBank: A public repository for sharing mass spectral data for life sciences, H. Horai, M. Arita, S. Kanaya, Y. Nihei, T. Ikeda, K. Suwa. Y. Ojima, K. Tanaka, S. Tanaka, K. Aoshima, Y. Oda, Y. Kakazu, M. Kusano, T. Tohge, F. Matsuda, Y. Sawada, M. Yokota Hirai, H. Nakanishi, K. Ikeda, N. Akimoto, T. Maoka, H. Takahashi, T. Ara, N. Sakurai, H. Suzuki, D. Shibata, S. Neumann, T. Iida, K. Tanaka, K. Funatsu, F. Matsuura, T. Soga, R. Taguchi, K. Saito and T. Nishioka, J. Mass Spectrom., 45 (2010) 703-714.

[3] Hummel J, Strehmel N, Selbig J, Walther D, Kopka J. Decision tree supported sub-structure prediction of metabolites from GC-MS profiles. Metabolomics, 6 (2010) 322-333.

}
\author{
Xavier Domingo-Almenara. xavier.domingo@urv.cat
}

\seealso{
\code{\link{newExp}}, \code{\link{setChrmMethod}}, \code{\link{setBaitPar}}
}
\examples{

# Library-driven deconvolution from an experiment created by \code{\link{newExp}}.

# ex <- newExp(instrumental="path")

# Now, we have to set up a chromatographic method:

Chrm.PLASMA <- setChrmMethod(method="alk.var5", 
rt=c(8.39, 10.759, 13.271, 15.604, 17.751, 19.685, 
21.471, 23.126, 24.645), ri=c(1225.27, 1326.95, 1526.31, 
1729.57, 1928.52, 2131.7, 2342.06, 2548.56, 2739.86), 
name="Test Chrm Method")

# The following will set BaitMet for analyzing the chromatograms
# without taking into account the masses 1:69,73:75,147:149, 
# (since is the mass range of the Golm Metabolome Database)
# with a minimum peak width of 2 seconds.

ext.par <- setBaitPar(ri.error=0.05, min.peak.width=2,  
min.peak.height=1000, noise.threshold=100, 
avoid.processing.mz=c(1:69,73:75,147:149))

# An now deconvolve the compounds in the samples:
# ex <- decBaitMet(ex, ext.par, chrom.method=Chrm.PLASMA)

}


