% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getLRinfosummary.R
\name{getLRinfosummary}
\alias{getLRinfosummary}
\title{Summarize LR info (expression, interaction score, and ligand score)}
\usage{
getLRinfosummary(exprinfo, LRIscorelist, Lscoredf, lr_network, cl_to, condpair)
}
\arguments{
\item{exprinfo}{output of "getExprInfo".}

\item{LRIscorelist}{output of "getLRIScore".}

\item{Lscoredf}{output of "getLScore".}

\item{lr_network}{a ligand_receptor network matrix with columns "L", "Lgene", "R", and "Rgene".}

\item{cl_to}{receiver cell-type/cluster.}

\item{condpair}{sample conditions under comparison, in the form of "cond1_vs_cond2".}
}
\value{
A list of data.frames of LR-info for each sender cell-type/cluster to the specified receiver cell-type/cluster
}
\description{
Summarize LR info (expression, interaction score, and ligand score)
}
