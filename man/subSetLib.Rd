\name{subSetLib}
\alias{subSetLib}
\title{Creation of a subset MS library}
\description{Creates a new MS library containing a subset of a larger MS library}
\usage{
subSetLib(database, indexes)
}
\arguments{
  \item{database}{A mass-spectral/retention index library.}
  \item{indexes}{The indexes to be included in this sub-set version (see Details/Examples)}
}
\details{
This function is useful to reduce the library to a number of compounds of interest, for example, those having a KEGG number (see examples below)
	
See BaitMet vignette for more details. To open the vignette, execute the following code in R:
vignette("BaitMetManual", package="baitmet")
}
\value{
The function returns the new subset library, see the examples below.
}
\author{
Xavier Domingo-Almenara. xavier.domingo@urv.cat
}

\examples{
# this function is useful if we want to reduce a 
# bigger MS library to a certain compounds of interest.
# For example, we can reduce the MS library included
# in BaitMet ("mslib" object containing the MassBank 
# library), to a smaller subset containing only those
# compounds with KEGG number:

kegg.ind <- which(lapply(mslib@database, function(x) x$KEGG)!="")
mslib.kegg <- subSetLib(mslib, kegg.ind)

# This is a naive example, as in fact, all the compounds
# in the mslib object have a KEGG number. However, this
# not occurs in the Golm Metabolome Database

}


