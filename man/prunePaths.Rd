\name{prunePaths}
\alias{prunePaths}
\title{Summarize the paths obtained by clipper according to their similarity.}
\description{This function allows the user to chose only one
  representant of those paths that have more than 1-thr similarity. The
  best scoring path is choosen.}
\usage{prunePaths(pathSummary, thr=NULL, clust=NULL, sep=";")}
\arguments{
  \item{pathSummary}{a matrix resulting from clipper function.}
  \item{thr}{a dissimilarity threshold. NULL means no pruning.}
  \item{clust}{filename where path-cluster is saved. NULL means no
    cluster saved.}
  \item{sep}{the separator to split genes for similarity
    computation. Default = ;}
}

\value{a matrix}

\seealso{\code{\link{clipper}}}

\examples{
toyEx <- matrix(c(1,1,5,3,5,2,5,3,8.2,3,2,1,0.3,0.1,2,1,"1;2;3;4;5","1;2;3","1;2;3;4;5","1;2;3","1;2;3;4;5","1;2;3","1;2;3;4;5","1;2;3"),2,12)
row.names(toyEx) <- c("1;5","1;3")
toyEx
prunePaths(toyEx, thr=0.1)
}