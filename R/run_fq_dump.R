#' \code{Run fastq-dump.}
#'
#' fastq-dump to dump SRA file.
#'
#' see more detail about SRA with Aspera downloading:
#' \url{http://www.ncbi.nlm.nih.gov/books/NBK158899/#SRA_download.downloading_sra_data_using}
#'
#' @param inputdf sra files in df as input. no need path info. [df, df["sra"]]
#' @param dumppath The absolute path of the SRA files to dump. [chr, "largedata/dump/"]
#' @param slurmsh File name of the output shell command. [chr, ="slurm-script/run_dump_"]
#' @param rmsra Remove the original SRA file after dumpping. [logical, =TRUE]
#' @param email Your email address that farm will email to once the job was done/failed. [chr, "y@gmail.com"]
#' @param gzip GZIP the fastq files. [logical, =FALSE]
#' @param runinfo  Parameters control the array job partition.
#' A vector of c(TRUE, "bigmemh", "8196", "1"): 1) run or not, 2) -p partition name, 3) --mem, adn 4) --ntasks.
#'
#' @return return a single shell script to run.
#'
#' @examples
#' ## run a single node job:
#' run_fq_dump(filepath="/group/jrigrp4/BS_teo20/WGBS",
#'             slurmsh="slurm-script/dump_WGBS.sh", rmsra=TRUE, email=NULL)
#'
#' ##  run array job:
#' run_fq_dump2(filepath="test", rmsra=TRUE, gzip=TRUE, email=NULL, run=c(TRUE, "bigmemh", "8196", "1"))
#'
#' @export
run_fq_dump <- function(inputdf, dumppath, rmsra=FALSE, gzip=FALSE, email=NULL,
                        slurmsh="slurm-script/run_dump_", runinfo=c(FALSE, "bigmemh", 5, "5G", "16:00:00")){

  files <- inputdf$sra
  dir.create("slurm-script", showWarnings = FALSE)
  for(i in 1:length(files)){

    shid <- paste0(slurmsh, i, ".sh")
    cmd1 <- paste0("cd ", dumppath)
    cmd2 <- paste0("fastq-dump --split-spot --split-3 -A ", files[i])
    cmd <- c(cmd1, cmd2)
    if(rmsra){
      cmd3 <- paste0("rm ", files[i])
      cmd <- c(cmd, cmd3)
    }
    if(gzip){
      cmd4 <- paste0("gzip ", paste0(files[i], "_1.fastq"))
      cmd5 <- paste0("gzip ", paste0(files[i], "_2.fastq"))
      cmd <- c(cmd, cmd4, cmd5)
    }
    cat(cmd, file=shid, sep="\n", append=FALSE)
  }

  shcode <- paste0("module load SRAtoolkit/2.8; sh ",
                   slurmsh, "$SLURM_ARRAY_TASK_ID.sh", sep="\n")
  set_array_job(shid=paste0(slurmsh, "array.sh"),
                shcode=shcode, arrayjobs=paste("1", length(files), sep="-"),
                wd=NULL, jobid="aspera", email=email, runinfo=runinfo)
}
