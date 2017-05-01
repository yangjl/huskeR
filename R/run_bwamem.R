#' \code{Run BWA-mem job on HPC}
#'
#' Heng Li's bwa manual: \url{http://bio-bwa.sourceforge.net/bwa.shtml}
#'
#' HCC help doc about running BWA Commands
#' \url{https://hcc-docs.unl.edu/display/HCCDOC/Running+BWA+Commands}
#'
#' BWA index:
#'   bwa index input_reference.fasta index_prefix
#'
#' see more detail about fermikit by Li, Heng:
#' \url{https://github.com/lh3/fermikit}
#'
#' @param fq An input data.frame for fastq files. Must contains fq1, fq2 and out.
#' @param kitpath The absolute or relative path of the fermi.kit directory that can invoke the pipeline.
#' @param genome The full path of genome with bwa indexed reference fasta file.
#' @param s Approximate genome size, default=3g.
#' @param t Number of threads, default=16.
#' @param l Primary read length, default=100.
#' @param arrayjobs A character specify the number of array you try to run, i.e. 1-100.
#' @param jobid The job name show up in your sq NAME column.
#' @param email Your email address that farm will email to once the job was done/failed.
#'
#' @return return a batch of shell scripts.
#'
#' @examples
#' fq <- data.frame(fq1=c("f_1.fq", "t_1.fq"), fq2=c("f_1.fq", "t_2.fq"), out=c("t1", "t2"))
#' run_fermikit(fq, kitpath="/home/jolyang/bin/fermikit/",
#' genome="/home/jolyang/dbcenter/AGP/AGPv2", s='3g', t=16, l=100, arrayjobs="1-2",
#' jobid="fermi", email=NULL)
#'
#' @export
run_bwamem <- function(fq,
                         kitpath="/home/jolyang/bin/fermikit",
                         s='3g', t=16, l=100,
                         arrayjobs="1-2",
                         jobid="fermi",
                         email=NULL){

  # create dir if not exist
  dir.create("slurm-script", showWarnings = FALSE)
  for(i in 1:nrow(fq)){

    shid <- paste0("slurm-script/run_bwamem_", i, ".sh")


    #-t INT     number of threads [1]
    #-M         mark shorter split hits as secondary (for Picard/GATK compatibility)
    #-T INT     minimum score to output [30]
    # genome ftp://ftp.ensemblgenomes.org/pub/plants/release-22/fasta/zea_mays/dna/README

    rg <- paste0("\'@RG\\tID:", fq$group[i], "\\tSM:", fq$sample[i],
                 "\\tPL:", fq$PL[i], "\\tLB:", fq$LB[i], "\\tPU:", fq$PU[i], "\'")

    cmd0 <- "module load bwa samtools"
    cmd1 <- paste("bwa mem -t", t, #mark shorter split hits as secondary
                  "-T", minscore, #minimum score to output [30]
                  "-R", rg, #read group header line such as '@RG\tID:foo\tSM:bar' [null]
                   indexfile, fq1, fq2,
                   "| samtools view -bSh - > $outfile.bam")
    #cmd2 <- paste0("samtools index sorted.$output.bam")
    #cmd3 <- paste0("samtools sort -m 10G -@ 2 $bam sorted.$output)
    cat(cmd1, file=shid, sep="\n", append=FALSE)
  }

  shcode <- paste("sh slurm-script/run_bwamem_$SLURM_ARRAY_TASK_ID.sh", sep="\n")

  set_array_job(shid="slurm-script/run_bwamem_array.sh",
                shcode=shcode, arrayjobs=arrayjobs,
                wd=NULL, jobid=jobid, email=email)
}


#http://gatkforums.broadinstitute.org/gatk/discussion/2799/howto-map-and-mark-duplicates
cat(paste("### Generate a SAM file containing aligned reads"),
    paste("bwa mem -M -R", rg, "-T", minscore, "-t", run[3], ref.fa, fq$fq1[i], fq$fq2[i], ">", aligned_sam),
    paste(""),

    ### http://broadinstitute.github.io/picard/
    paste0("java -Xmx", floor(as.numeric(run[4])/1024),
           "g -jar ", picardpwd, " SortSam \\"),
    paste0("    INPUT=", aligned_sam, " \\"),
    paste0("    OUTPUT=", sorted_bam, " \\"),
    "    SORT_ORDER=coordinate",
    paste("#rm", aligned_sam),
    paste(""),
    file=shid, sep="\n", append=TRUE)
message("###>>> set up BWA mem and then sort to bam using picard-tools!")
