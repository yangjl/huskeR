% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_GBSv2.R
\name{run_GBSv2}
\alias{run_GBSv2}
\alias{GBSSeqToTag}
\alias{TagExportToFastq}
\alias{run_bowtie2}
\alias{SAMToGBSdb}
\alias{DiscoverySNP}
\alias{SNPQuality}
\alias{ProductionSNPCaller}
\title{\code{Run TASSEL GBSv2 on slurm-based HPC}}
\usage{
run_GBSv2(outdir = "largedata/gbs", shfile, fqdir, keyfile, bt2_idx, sam,
  ref, snpqc_file, h5, db = NULL, mem, cpu, kmerlen = 64,
  enzyme = "ApeKI", mnlcov = 0.05, mnmaf = 0.01,
  deleteolddata = "true", production = FALSE, seq2tag = TRUE,
  tag2fq = TRUE, bt2 = TRUE, snpcall = TRUE)

GBSSeqToTag(shfile, mem, fqdir, db, keyfile, kmerlen, enzyme)

TagExportToFastq(shfile, mem, db, tagexport, mindepth = 1)

run_bowtie2(shfile, cpu, bt2_idx, tagexport, sam)

SAMToGBSdb(shfile, mem, sam, db)

DiscoverySNP(shfile, mem, db, ref, mnlcov = 0.05, mnmaf = 0.001,
  deleteolddata)

SNPQuality(shfile, mem, db, snpqc_file)

ProductionSNPCaller(shfile, mem, db, enzyme, fqdir, keyfile, kmerlen, h5)
}
\arguments{
\item{shfile}{A shell script in "slurm-script/". [chr, "slurm-script/run_it.sh"]}

\item{fqdir}{Absolute dir for the fastq files, can be nested within the dir. [chr, "my/fq/filedir"]}

\item{keyfile}{The key file. The GBSv2 key file has 4 required headers: Flowcell, Lane, Barcode and FullSampleName. [chr, "data/file.txt]}

\item{bt2_idx}{The idx file for bowtie2. [chr, "/lustre/work/jyanglab/jyang21/dbcenter/AGP/AGPv4/AGPv4_bt2_index/AGPv4_bt2_index"]}

\item{sam}{The output sam file. [chr, "tagexport/tagsForAlignFullvs_F1s.sam"]}

\item{ref}{Fastq file of the ref. [chr, "/lustre/work/jyanglab/jyang21/dbcenter/AGP/AGPv4/AGPv4_bt2_index/AGPv4_bt2.fa"]}

\item{snpqc_file}{SNP qc file. [chr, "largedata/gbs/outputStats.txt"]}

\item{h5}{Output hdf5 format. [chr, "gbs/discovery/output.h5"]}

\item{db}{Relative file name of the db. [chr, =NULL, or "data/gbs.db"]}

\item{mem}{Memory usage. [num, 50]}

\item{cpu}{Number of CPUs to use. [num, 50]}

\item{kmerlen}{Kmer length [num, =64]}

\item{enzyme}{Enzyme used for cutting. [chr, ="ApeKI"]}

\item{mnlcov}{mnLCov. [num, =0.05]}

\item{mnmaf}{Min MAF. [num, =0.001]}

\item{deleteolddata}{Delete Old Data. [chr, ="true"]}

\item{production}{SNP production. [logical, =F]}

\item{seq2tag}{Generated tags from fastq files. [logical, =T]}

\item{tag2fq}{extract tags from db for alignment to reference. [logical, =T]}

\item{bt2}{Run bt2 alignment. [logical, =T]}

\item{snpcall}{Call SNP. [logical, =T]}

\item{tagexport}{Fastq.gz file to export tags. [chr, "tagexport/tagsForAlign.fa.gz"]}

\item{mindepth}{The minimum count of reads for a tag to be output. [num, =1]}
}
\value{
return a batch of shell scripts.
}
\description{
a single slurm job.
help doc: https://bitbucket.org/tasseladmin/tassel-5-source/wiki/Tassel5GBSv2Pipeline
}
\examples{


}
