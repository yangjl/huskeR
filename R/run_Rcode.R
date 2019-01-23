#' \code{Run array job of R codes}
#'
#'  Run array rcodes according to the parameters in the inputdf file.
#'  The array job number = nrow(inputdf).
#'
#' @param inputdf An input data.frame, with columns of file and out. [df, cols=file, out]
#'                Note, row parameter "j" will be passed to each array job!
#' @param outdir The dir of shell files. [chr, "largedata/"]
#' @param cmdno Number of commands to excute in each array. [num, =1]
#' @param rcodes The abosulte path of your R codes to run. [chr, ="lib/C_format.R"]
#' @param arrayshid The sbatch id. [chr, ="slurm-script/run_bcf_query_array.sh"]
#' @param email Your email address that farm will email to once the jobs were done/failed. [chr, =NULL]
#' @param runinfo [vector, runinfo = c(FALSE, "bigmemh", 5, "5G", "16:00:00")]
#'
#' @return return a batch of shell scripts.
#'
#' @examples
#' sh <- paste0('R --no-save --args ', j, ' < ', rcodes)
#' run_Rcode(inputdf=data.frame(file=1:11, out=10), outdir="slurm-script", cmdno=10,
#'            rcodes = "lib/C_format.R", arrayshid = "slurm-script/run_rcode_array.sh",
#'            email=NULL, runinfo = c(FALSE, "bigmemh", 1))
#'
#' @export
run_Rcode <- function(
    inputdf, outdir, cmdno=1,
    rcodes = "lib/C_format.R",
    arrayshid = "slurm-script/run_bcf_query_array.sh",
    email=NULL, runinfo = c(FALSE, "bigmemh", 5, "5G", "16:00:00", 1)
){

    #runinfo <- get_runinfo(runinfo)
    #### create dir if not exist
    dir.create("slurm-script", showWarnings = FALSE)
    dir.create(outdir, showWarnings = FALSE)

    tot <- ceiling(nrow(inputdf)/cmdno)
    for(j in 1:tot){

        shid <- paste0(outdir, "/run_rcode_", j, ".sh")

        ##chr:start-end
        #sh1 <- paste("cd", outdir)
        sh <- paste0('R --no-save --args ', j, ' < ', rcodes)

        cat(paste("### run Rcode", Sys.time(), sep=" "),
            sh,
            file=shid, sep="\n", append=FALSE)
    }

    shcode <- paste0("module load R; sh ", outdir, "/run_rcode_$SLURM_ARRAY_TASK_ID.sh")
    set_array_job(shid=arrayshid, shcode=shcode, arrayjobs=paste("1", tot, sep="-"),
                  wd=NULL, jobid="rcode", email=email, runinfo=runinfo)
}

