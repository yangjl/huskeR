#' \code{Run kallisto quant.}
#'
#' kallisto quant
#'
#' see more detail about fastq-dump with Aspera downloading:
#' \url{https://pachterlab.github.io/kallisto/manual}
#' 
#' @param df Input files, [data.table, col: outdir, fq (SE) or fq1 and fq2 (PE)]
#' @param single Single (or PE) end fastq [TRUE(FALSE)]. 
#' @param l Estimated average fragment length [number, default=200].
#' @param s Estimated standard deviation of fragment length [number, default=20].
#' @param run  Parameters control the array job partition.
#' A vector of c(TRUE, "bigmemh", "1", "8:00:00"): 1) run or not, 2) -p partition name, 3) --ntasks and 4) --times.
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
run_kallisto <- function(df, single=TRUE, l=200, s=20,
                         idx="~/dbcenter/AGP/AGPv2/ZmB73_5a_WGS_exons.kallisto.idx",
                         email=NULL,
                         jobid="run_kquant", runinfo=c(TRUE, "bigmemh", "1", "8:00:00")){
  
  #files <- list.files(path=filepath, pattern="sra$")
  dir.create("slurm-script", showWarnings = FALSE)
  for(i in 1:nrow(df)){
    
    shid <- paste0("slurm-script/run_kquant_", i, ".sh")
    dir.create(df$outdir[i], showWarnings = FALSE)
    if(single){
      cmd <- paste("kallisto quant -i", idx, "--plaintext -o", df$outdir[i],
                     "--single -l", l, "-s", s, df$fq[i])
    }else{
      cmd <- paste("kallisto quant -i", idx, "--plaintext -o", df$outdir[i],
                   df$fq1[i], df$fq2[i])
    }
    
    cat(cmd, file=shid, sep="\n", append=FALSE)
  }
  
  shcode <- paste("sh slurm-script/run_kquant_$SLURM_ARRAY_TASK_ID.sh", sep="\n")
  myshid <- "slurm-script/run_kquant_array.sh"
  
  set_array_job(shid=myshid,
                shcode=shcode, arrayjobs=paste("1", nrow(df), sep="-"),
                wd=NULL, jobid=jobid, email=email, runinfo = runinfo)
}
