#' \code{Run Aspera Connect to download from SRA.}
#'
#' Downloading SRA using 'ascp' utility or Aspera Connect.
#' How to download and install 'ascp':
#' \url{http://downloads.asperasoft.com/en/downloads/8?list}
#'
#' see more detail about SRA with Aspera downloading:
#' \url{http://www.ncbi.nlm.nih.gov/books/NBK158899/#SRA_download.downloading_sra_data_using}
#'
#' @param sradf An input data.frame for SRA ids. Must contains column: SRR. [df["SRR"]]
#' @param maxspeed The max speed for downloading. [chr, =200m]
#' @param outdir The output directory. [chr, ="."]
#' @param cmdno Number of commands to excute in each array job. [num, =1]
#' @param email Your email address that farm will email to once the job was done/failed. [chr, =NULL]
#' @param runinfo [vector, = c(FALSE, "bigmemh", 5, "5G", "16:00:00")]
#'
#' @return return a batch of shell scripts.
#'
#' @examples
#' sra <- data.frame(SRR=c("ERR877647", "ERR877648"),
#' SRX=c( "ERX957210", "ERX957211"),
#' pid=c( "1_Base1_Bbreve-sc-2188486", "P1_ECvsrS_1-sc-2201977"))
#' run_aspera(sra, maxspeed="200m", outdir=".", arrayjobs="1-2", jobid="aspera", email=NULL)
#'
#' @export
run_aspera <- function(sradf,
                       maxspeed="100m",
                       outdir=".",
                       cmdno=1,
                       email=NULL,
                       runinfo){

    # create dir if not exist
    dir.create("slurm-script", showWarnings = FALSE)
    tot <- ceiling(nrow(sradf)/cmdno)

    for(j in 1:(tot-1)){
        srow <- (j-1)*cmdno + 1
        erow <- j*cmdno
        set_aspera(sradf, srow, erow, j, maxspeed, outdir)
    }
    set_aspera(sradf, srow=cmdno*(tot-1)+1, erow=nrow(sradf), j=tot, maxspeed, outdir)

    shcode <- paste("sh slurm-script/run_aspera_$SLURM_ARRAY_TASK_ID.sh", sep="\n")
    set_array_job(shid="slurm-script/run_aspera_array.sh",
                  shcode=shcode, arrayjobs=paste("1", tot, sep="-"),
                  wd=NULL, jobid="aspera", email=email, runinfo=runinfo)
}

#' \code{Set up aspera code in a job shell}
#'
#' @param srow Start row of the df. [num, 1]
#' @param erow End row of the df. [num, 10]
#'
#' @rdname run_aspera
set_aspera <- function(sradf, srow, erow, j, maxspeed, outdir){

    subdf <- sradf[srow:erow, ]

    ### setup the front line:
    shid <- paste0("slurm-script/run_aspera_", j, ".sh")
    cat(paste("## [huskeR: run_aspera]: set up aspera download", Sys.time(), sep=" "),
        file=shid, sep="\n", append=FALSE)

    for(i in 1:nrow(subdf)){

        # ascp root: vog.hin.mln.ibcn.ptf@ptfnona:
        #Remainder of path:
        #    /sra/sra-instant/reads/ByRun/sra/{SRR|ERR|DRR}/<first 6 characters of accession>/<accession>/<accession>.sra
        #Where
        #{SRR|ERR|DRR} should be either ‘SRR’, ‘ERR’, or ‘DRR’ and should match the prefix of the target .sra file
        #ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -QT -l 100m
        # anonftpftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/SRR161/SRR1610960/SRR1610960.sra
        sraid <- subdf$SRR[i]
        id1 <- substr(sraid, start=1, stop=3)
        id2 <- substr(sraid, start=1, stop=6)
        # genome ftp://ftp.ensemblgenomes.org/pub/plants/release-22/fasta/zea_mays/dna/README
        #http://www.ncbi.nlm.nih.gov/books/NBK158899/#SRA_download.downloading_sra_data_using

        cmd <- paste0("ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -QT -l ", maxspeed,
                      " anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/",
                      id1, "/", id2, "/", sraid, "/", sraid, ".sra ", outdir)
        cat(c(cmd, ""), file=shid, sep="\n", append=TRUE)
    }
}
