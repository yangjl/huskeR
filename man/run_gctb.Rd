% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_gctb.R
\name{run_gctb}
\alias{run_gctb}
\title{\code{Run GCTB in a wrapper developed by Jinliang Yang Lab at UNL.}}
\usage{
run_gctb(df, bayes, shid = "run_gtcb", email = NULL,
  runinfo = c(TRUE, "jclarke", "1", "8G", "8:00:00"))
}
\arguments{
\item{df}{Input files with all the parameters for the Bayesian alphabet. [data.frame, col: bfile, ... (see examples for more details)]}

\item{shid}{The prefix of the shell code name in slurm-script folder.}

\item{email}{Your email address that farm will email to once the job was done or failed.}

\item{runinfo}{Parameters control the array job partition.
A vector of c(TRUE, "batch", "1", "8G", 8:00:00"): 1) run within R or not, 2) -p partition name,
3) --ntasks, 4) --mem and 5) --times.}
}
\value{
return a batch of shell scripts to run.
}
\description{
you need to download and install GTCB before running this wrapper.
}
\details{
see more detail about GTCB:
\url{http://cnsgenomics.com/software/gctb/#Overview}
}
\examples{
### build kallisto idex
 kallisto index -i ZmB73_5a_WGS_exons.kallisto.idx ZmB73_5a_WGS_exons.fasta.gz
 kallisto quant -i ~/dbcenter/AGP/AGPv2/ZmB73_5a_WGS_exons.kallisto.idx --plaintext -o .
--single -l 200 -s 20 SRR957415.sra.fastq SRR957416.sra.fastq > test.txt

}
