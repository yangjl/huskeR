#' \code{Run GCTB in a wrapper developed by Jinliang Yang Lab at UNL.}
#'
#' you need to download and install GTCB before running this wrapper.
#'
#' see more detail about GTCB:
#' \url{http://cnsgenomics.com/software/gctb/#Overview}
#'
#' @param df Input files with all the parameters for the Bayesian alphabet. [data.frame, col: bfile, ... (see examples for more details)]
#' @param shid The prefix of the shell code name in slurm-script folder.
#' @param email Your email address that farm will email to once the job was done or failed.
#' @param runinfo  Parameters control the array job partition.
#' A vector of c(TRUE, "batch", "1", "8G", 8:00:00"): 1) run within R or not, 2) -p partition name,
#' 3) --ntasks, 4) --mem and 5) --times.
#'
#' @return return a batch of shell scripts to run.
#'
#' @examples
#' ### build kallisto idex
#'  kallisto index -i ZmB73_5a_WGS_exons.kallisto.idx ZmB73_5a_WGS_exons.fasta.gz
#'  kallisto quant -i ~/dbcenter/AGP/AGPv2/ZmB73_5a_WGS_exons.kallisto.idx --plaintext -o .
#' --single -l 200 -s 20 SRR957415.sra.fastq SRR957416.sra.fastq > test.txt
#'
#' @export
run_gctb <- function(df, bayes, shid="run_gtcb",
                     email=NULL, runinfo=c(TRUE, "jclarke", "1", "8G", "8:00:00")){

    # create the slurm-script folder if it is not there.
  dir.create("slurm-script", showWarnings = FALSE)

    # gctb --bfile /common/HapMap3/set282/allchr__maf001_geno1 --pheno largedata/gctb_pheno/Zn66.phen
    # --bayes NS --pi 0.05 --hsq 0.5 --S 0 --wind 0.1 --chain-length 410000 --burn-in 10000
    # --out largedata/bayesNS/test > largedata/bayesNS/test.log
  for(i in 1:nrow(df)){

    shid_i<- paste0("slurm-script/", shid, "_", i, ".sh")
    #dir.create(df$out[i], showWarnings = FALSE)

    if(bayes == "S"){
        cmd <- paste("gctb --bfile", df$bfile[i], "--pheno", df$pheno[i],
                     "--bayes", bayes,
                     "--pi", df$pi[i],
                     "--hsq", df$hsq[i],
                     "--S", df$S[i],
                     "--chain-length", df$chainlength[i],
                     "--burn-in", df$burnin[i],
                     "--out", df$out[i],
                     ">", df$log)

    }else if(bayes == "NS"){
      cmd <- paste("gctb --bfile", df$bfile[i], "--pheno", df$pheno[i],
                   "--bayes", bayes,
                   "--pi", df$pi[i],
                   "--hsq", df$hsq[i],
                   "--S", df$S[i],
                   "--wind", df$wind[i],
                   "--chain-length", df$chainlength[i],
                   "--burn-in", df$burnin[i],
                   "--out", df$out[i],
                   ">", df$log[i])

    }else{
        message(sprintf("[Warning!] stay tuned, [%s] have not been wrapped yet!", bayes))
    }
    cat(cmd, file=shid_i, sep="\n", append=FALSE)
  }

  shcode <- paste0("sh slurm-script/", shid, "_$SLURM_ARRAY_TASK_ID.sh", sep="\n")
  myshid <- paste0("slurm-script/", shid, "_array.sh")

  set_array_job(shid=myshid,
                shcode=shcode, arrayjobs=paste("1", nrow(df), sep="-"),
                wd=NULL, jobid=shid, email=email, runinfo = runinfo)
}
