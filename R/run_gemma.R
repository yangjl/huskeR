#' \code{Run GenSel on slurm-based HPC}
#'
#' GenSel4R
#'
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
run_gemma <- function(
  inputdf, cmdno=1,
  shid = "slurm-script/run_gensel_array.sh",
  remove=TRUE,
  email=NULL, runinfo = c(FALSE, "batch", "1", "2G", "1:00:00")
){

  #### create dir if not exist
  dir.create("slurm-script", showWarnings = FALSE)
  dir.create(inpdir, showWarnings = FALSE)

  ### setup shell id
  set_gemma(inputdf, cmdno, workdir)

  #sh0 <- "module load compiler/gcc/6.1"
  shcode <- paste0("sh slurm-script/run_gemma_$SLURM_ARRAY_TASK_ID.sh")
  set_array_job(shid=shid, shcode=shcode, arrayjobs=paste("1", ceiling(nrow(inputdf)/cmdno), sep="-"),
                wd=NULL, jobid="gensel", email=email, runinfo=runinfo)
}

##### set up serial gemma!
set_gemma <- function(inputdf, cmdno, workdir){
    ####
    for(j in 1:(ceiling(nrow(inputdf)/cmdno))){

        srow <- cmdno*(j-1) + 1
        erow <- cmdno*j
        if(erow > nrow(inputdf)){
            erwo <- nrow(inputdf)
        }

        shid <- paste0("slurm-script/run_gemma_", j, ".sh")

        sh1 <- c("cd ", workdir)
        sh2 <- paste0("gemma -bfile ", inputdf$bfile[srow:erow],
                      " -n ", inputdf$trait[srow:erow],
                      " -k ", inputdf$kin[srow:erow],
                      " -lmm 4 -o ", inputdf$output[srow:erow])

        cat(paste("### run gemma", Sys.time(), sep=" "),
            c(sh1, sh2),
            file=shid, sep="\n", append=FALSE)
    }
}



