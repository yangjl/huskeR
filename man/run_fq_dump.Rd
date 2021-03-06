% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_fq_dump.R
\name{run_fq_dump}
\alias{run_fq_dump}
\title{\code{Run fastq-dump.}}
\usage{
run_fq_dump(inputdf, runnum = 100, dumppath, rmsra = FALSE,
  gzip = FALSE, email = NULL, slurmsh = "slurm-script/run_dump_",
  runinfo = c(FALSE, "bigmemh", 5, "5G", "16:00:00"))
}
\arguments{
\item{inputdf}{sra files in df as input. no need path info. [df, df["sra"]]}

\item{runnum}{Number of running jobs in parallel. [int, =100]}

\item{dumppath}{The absolute path of the SRA files to dump. [chr, "largedata/dump/"]}

\item{rmsra}{Remove the original SRA file after dumpping. [logical, =TRUE]}

\item{gzip}{GZIP the fastq files. [logical, =FALSE]}

\item{email}{Your email address that farm will email to once the job was done/failed. [chr, "y@gmail.com"]}

\item{slurmsh}{File name of the output shell command. [chr, ="slurm-script/run_dump_"]}

\item{runinfo}{Parameters control the array job partition.
A vector of c(TRUE, "bigmemh", "8196", "1"): 1) run or not, 2) -p partition name, 3) --mem, adn 4) --ntasks.}
}
\value{
return a single shell script to run.
}
\description{
fastq-dump to dump SRA file.
}
\details{
see more detail about SRA with Aspera downloading:
\url{http://www.ncbi.nlm.nih.gov/books/NBK158899/#SRA_download.downloading_sra_data_using}
}
\examples{
## run a single node job:
run_fq_dump(filepath="/group/jrigrp4/BS_teo20/WGBS",
            slurmsh="slurm-script/dump_WGBS.sh", rmsra=TRUE, email=NULL)

##  run array job:
run_fq_dump2(filepath="test", rmsra=TRUE, gzip=TRUE, email=NULL, run=c(TRUE, "bigmemh", "8196", "1"))

}
