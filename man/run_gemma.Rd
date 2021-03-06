% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_gemma.R
\name{run_gemma}
\alias{run_gemma}
\title{\code{Run serial gemma jobs on slurm-based HPC}}
\usage{
run_gemma(inputdf, cmdno = 1, workdir = "largedata/",
  shid = "slurm-script/run_gemma_array.sh", email = NULL,
  runinfo = c(FALSE, "batch", "1", "2G", "1:00:00"))
}
\arguments{
\item{inputdf}{An input data.frame for this function, including all the necessary pheno, geno file pathes,
and other parameters. [data.frame, [df, cols=bfile, trait, kin, output]]}

\item{cmdno}{Number of commands per CPU, i.e. number of rows per inputdf. [num, default=1].}

\item{workdir}{The dir that will be cd to. [chr, default="largedata/"].}

\item{shid}{The shell ID for this job. [chr, default is "slurm-script/run_gensel_array.sh"].}

\item{email}{Your email address that farm will email to once the jobs were done/failed.
[chr, default is NULL]}

\item{runinfo}{Parameters specify the array job partition information.
A vector of c(FALSE, "batch", "1", "2G", "8:00:00"): 1) run or not, default=FALSE
2) -p partition name, default=batch, 3) --ntasks=3, 4) --mem=30G, and 5) --time=8:00:00.}
}
\value{
return a batch of shell scripts.
}
\description{
Gemma Zhou lab: http://www.xzlab.org/software.html
}
\examples{


}
