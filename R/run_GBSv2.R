#' \code{Run TASSEL GBSv2 on slurm-based HPC}
#'
#' a single slurm job.
#' help doc: https://bitbucket.org/tasseladmin/tassel-5-source/wiki/Tassel5GBSv2Pipeline
#'
#'
#' @param inputdf An input data.frame for this function, including all the necessary pheno, geno file pathes,
#' and other parameters. [data.frame, see exaples]
#' @param cv Cross-validation experiment [TRUE/FALSE, default=FALSE].
#'           If TRUE, inputdf must contain `trainpheno` and `testpheno`; if FALSE, only need `pheno`.
#' @param bayesType Analysis types for Bayes Methods. [chr, default=BayesC; BayesA, BayesB, BayesC, BayesCPi, RBR]
#' @param inpdir Dir for the GS .inp file. [chr, default="largedata/"].
#' @param cmdno Number of commands per CPU, i.e. number of rows per inputdf. [num, default=1].
#' @param remove Whether to remove the intermediate GS files or not. [TRUE/FALSE, default is TRUE].
#' @param shid The shell ID for this job. [chr, default is "slurm-script/run_gensel_array.sh"].
#' @param email Your email address that farm will email to once the jobs were done/failed.
#'           [chr, default is NULL]
#' @param runinfo Parameters specify the array job partition information.
#'           A vector of c(FALSE, "batch", "30G", 1, "8:00:00"): 1) run or not, default=FALSE
#'           2) -p partition name, default=batch, 3) --ntasks=3, 4) --mem=30G, and 5) --time=8:00:00.
#'
#' @return return a batch of shell scripts.
#'
#' @examples
#'
#' inputdf <- data.frame(pi=0.995,
#'    geno="/Users/yangjl/Documents/GWAS2_KRN/SNP/merged/geno_chr.newbin",
#'    trainpheno="/Users/yangjl/Documents/Heterosis_GWAS/pheno2011/reports/cd_GenSel_fullset.txt",
#'    testpheno="/Users/yangjl/Documents/Heterosis_GWAS/pheno2011/reports/cd_GenSel_fullset.txt",
#'    chainLength=1000, burnin=100, varGenotypic=1.4, varResidual=2,
#'    out="test_out")
#'
#' run_GenSel4(inputdf, cv=FALSE, inpdir="slurm-script/", cmdno=1,
#'             shid = "slurm-script/run_gensel_array.sh", remove=TRUE,
#'             email=NULL, runinfo = c(FALSE, "batch", 3, "30", "8:00:00") )
#'
#' @export
GBSSeqToTag <- function(){

}
#'
#'
run_GenSel4 <- function(
  inputdf, cv=FALSE, bayesType="BayesC", inpdir="largedata/", cmdno=1,
  shid = "slurm-script/run_gensel_array.sh",
  remove=TRUE,
  email=NULL, runinfo = c(FALSE, "batch", "30", 3, "8:00:00")
){

  #### create dir if not exist
  dir.create("slurm-script", showWarnings = FALSE)
  dir.create(inpdir, showWarnings = FALSE)

  for(i in 1:nrow(inputdf)){
    inpid <- paste0(inpdir, "/", inputdf$out[i], ".inp")
    ### output the inp file:

    if(cv){
      GS_cv_inp(
        inp= inpid, bayesType=bayesType, pi=inputdf$pi[i], geno=inputdf$geno[i],
        trainpheno=inputdf$trainpheno[i], testpheno=inputdf$testpheno[i],
        chainLength=inputdf$chainLength[i], burnin=inputdf$burnin[i],
        varGenotypic=inputdf$varGenotypic[i], varResidual=inputdf$varResidual[i]
      )
    }else{
      GS_regular_inp(
        inp= inpid, bayesType=bayesType, pi=inputdf$pi[i], geno=inputdf$geno[i],
        pheno=inputdf$trainpheno[i],
        chainLength=inputdf$chainLength[i], burnin=inputdf$burnin[i],
        varGenotypic=inputdf$varGenotypic[i], varResidual=inputdf$varResidual[i]
      )
    }

  }

  ### setup shell id
  set_GS(inputdf, cmdno, inpdir, remove)

  sh0 <- "module load compiler/gcc/6.1"
  shcode <- c(sh0, paste0("sh ", inpdir, "/run_gensel_$SLURM_ARRAY_TASK_ID.sh"))
  set_array_job(shid=shid, shcode=shcode, arrayjobs=paste("1", nrow(inputdf)/cmdno, sep="-"),
                wd=NULL, jobid="gensel", email=email, runinfo=runinfo)
}



#' \code{generate tags from fastq file}
#'
#' @param shfile A shell script.
#' @param mem Memory [num, =50]
#' @param fqdir Absolute dir for the fastq files, can be nested within the dir. [chr, "my/fq/filedir"]
#' @param db Relative file name of the db. [chr, "data/gbs.db"]
#' @param keyfile The key file. The GBSv2 key file has 4 required headers: Flowcell, Lane, Barcode and FullSampleName. [chr, "data/file.txt]
#'
#'
#' @return return a shell script.
#'
#' @rdname run_tassel
GBSSeqToTag <- function(shfile, mem=50, fqdir, db, keyfile){

  cat(paste("## Step1: GBSSeqToTag, Tassel5.2 GBSv2", Sys.time(), sep=" "),

      paste0("run_pipeline.pl -Xmx", mem, "g -fork1 -GBSSeqToTagDBPlugin -e ApeKI \\"),
      paste("-i", fqdir, " \\"),
      paste("-db", db, " \\"),
      paste("-k", keyfile, " \\"),

      "-c 10 \\",
      "-kmerLength 80 \\",
      "-minKmerL 20 \\",
      "-mnQS 20 \\",
      "-batchSize 8 \\",
      "-endPlugin -runfork1",

      file=shfile, sep="\n", append=FALSE
  )
}


#' \code{extract tags from db for alignment to reference}
#'
#' @param shfile A shell script.
#' @param tagexport Fastq.gz file to export tags. [chr, ""export/tagsForAlign_F1s.fa.gz"]
#' @param mindepth The minimum count of reads for a tag to be output. [num, =1]
#'
#' @return return a shell script.
#'
#' @rdname run_tassel
TagExportToFastq <- function(shfile, mem=50, db, tagexport, mindepth=1){

    cat(paste("## Step2: TagExportToFastq, Tassel5.2 GBSv2", Sys.time(), sep=" "),

        paste0("run_pipeline.pl -Xmx", mem, "g -fork1 -TagExportToFastqPlugin \\"),
        paste("-db", db, " \\"),
        paste("-o", tagexport, " \\"),
        paste("-c", mindepth, " \\"),
        "-endPlugin -runfork1",

        file=shfile, sep="\n", append=FALSE
    )
}


