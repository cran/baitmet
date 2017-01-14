\name{mzList}
\alias{mzList}
\title{MZ list}
\description{The list of compounds selective masses and their relative quantification for each sample in a given experiment}
\usage{
mzList(Experiment, by.area=TRUE)
}
\arguments{
  \item{Experiment}{A 'MetaboSet' S4 object containing the experiment data. The experiment has to be previously deconvolved, aligned and (optionally) identified.}
    \item{by.area}{if TRUE (default), the function outputs the quantification by the area of the mz of each compound. If FALSE, the function outputs the peak intensity.}
}
\details{
Returns an alignment table containing the list of the selected MZ and their relative quantification for each sample in a given experiment.}

\seealso{
\code{\link{quantSM}}, \code{\link{plotMZ}}

}