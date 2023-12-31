% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SCOM.R
\name{Overlap.hub}
\alias{Overlap.hub}
\title{Overlap.hub}
\usage{
Overlap.hub(hub, overlap.num = 1, type = c("equal", "least"))
}
\arguments{
\item{hub}{List object, the list of hubs.}

\item{overlap.num}{The minimum number of hubs existing in multiple lists of hubs}

\item{type}{The overlapped hubs in overlap.num hub lists ("equal") or at least overlap.num hub lists ("least").}
}
\value{
A vector: The overlapped hubs.
}
\description{
Identifying the overlap between multiple lists of hubs.
}
\examples{
hub <- list(c("ncRNA1", "ncRNA2", "ncRNA3", "ncRNA4"), c("ncRNA1", "ncRNA2", "ncRNA4", "ncRNA5"), c("ncRNA1", "ncRNA2", "ncRNA5", "ncRNA6"), c("ncRNA1", "ncRNA3", "ncRNA4", "ncRNA6")) 
hub_2 <- Overlap.hub(hub, overlap.num = 2, type = "least") 

}
\author{
Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng-Zhang-2})
}
