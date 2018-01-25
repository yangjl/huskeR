#' \code{Run VCAP array jobs on HPC}
#'
#' kallisto quant
#'
#' see more detail about fastq-dump with Aspera downloading:
#' \url{https://pachterlab.github.io/kallisto/manual}
#'
#' @param df Input files. [data.frame, cols: "bedfile1", "bedfile2", ..., "bedfilen", gz.lix, genome_kinship and res_kinship]
#'
#' @param runinfo  Parameters control the array job partition.
#' A vector of c(TRUE, "bigmemh", "1", "8:00:00"): 1) run or not, 2) -p partition name, 3) --ntasks and 4) --times.
#'
#' @return return a batch of shell scripts to run.
#'
#' @examples
#'
#'
#' @export
run_VCAP <- function(df,
                    email=NULL,
                    jobid="run_vcap", runinfo=c(TRUE, "jclarke", "1", "8G", "8:00:00")){

    #files <- list.files(path=filepath, pattern="sra$")
    dir.create("slurm-script", showWarnings = FALSE)
    for(i in 1:nrow(df)){

        shid <- paste0("slurm-script/run_vcap_", i, ".sh")

        bedidx <- grep("bedfile", names(df))

        cmd0 <- paste0("run_pipeline.pl -Xmx", runinfo[4], " -fork1 -importGuess", df$gz.lix[i])

        fid <- 2 #fork id
        cmd1 <- c()
        for(j in bedidx){
            # -fork2 -FilterSiteBuilderPlugin -bedFile path/to/myBedFile1.bed -endPlugin -input1 -KinshipPlugin -method Centered_IBS -endPlugin
            temp1 <- paste0("-fork", fid, " -FilterSiteBuilderPlugin -bedFile ", df[i, j],
                       " -endPlugin -input1 -KinshipPlugin -method Centered_IBS -endPlugin")
            # -fork3 -FilterSiteBuilderPlugin -bedFile path/to/myBedFile2.bed -endPlugin -input1 -KinshipPlugin -method Centered_IBS -endPlugin
            cmd1 <- c(cmd1, temp1)
            fid <- fid + 1
        }

        inputid <- 2
        cmd2 <- c()
        for(j in bedidx){
            # -fork4 -export myAptlyNamedCentered_IBS_kinship1 -exportType SqrMatrixBin -input2
            # -fork5 -export myAptlyNamedCentered_IBS_kinship2 -exportType SqrMatrixBin -input3
            temp2 <- paste0("-fork", fid, " -export kinship_", df[i, j], " -exportType SqrMatrixBin -input", inputid)
            cmd2 <- c(cmd2, temp2)
            inputid <- inputid + 1
            fid <- fid + 1
        }

        #-combine6 -input2 -input3 -SubtractDistanceMatrixPlugin -wholeMatrix myGenotypesWholeGenome_Centered_IBS_kinship.txt
        cmd3 <- paste0("-combine", fid)
        cmd4 <- paste0("-input", 2:(inputid -1), " ")
        cmd5 <- paste0("-SubtractDistanceMatrixPlugin -wholeMatrix ", df$genome_kinship[i])
        #-endPlugin -export myAptlyNamedCentered_IBS_kinshipRest -exportType SqrMatrixBin
        cmd6 <- paste0("-endPlugin -export ", df$res_kinship[i], " -exportType SqrMatrixBin")

        cat(cmd0, cmd1, cmd2, cmd3,
            cmd4, cmd5, cmd6,
            file=shid, sep=" ", append=FALSE)
    }

    shcode <- paste("module load java/1.8", "module load tassel/5.2", "sh slurm-script/run_vcap_$SLURM_ARRAY_TASK_ID.sh", sep="\n")
    myshid <- "slurm-script/run_vcap_array.sh"

    set_array_job(shid=myshid,
                  shcode=shcode, arrayjobs=paste("1", nrow(df), sep="-"),
                  wd=NULL, jobid=jobid, email=email, runinfo = runinfo)
}



