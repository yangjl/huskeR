#' \code{Run GATK array jobs on HPC}
#'
#' GATK Best Practices: recommended workflows for variant discovery analysis.
#'
#' see more detail about GATK:
#' \url{https://www.broadinstitute.org/gatk/guide/bp_step.php?p=1}
#'
#' idxing:
#' bwa index Zea_mays.AGPv2.14.dna.toplevel.fa
#'
#' module load java/1.8
#' module load bwa/0.7.9a
#'
#' local programs:
#' bwa Version: 0.7.5a-r405
#' picard-tools-2.1.1
#' GenomeAnalysisTK-3.5/
#'
#' @param inputdf An input data.frame for fastq files. Must contains fq1, fq2, out (and/or bam).
#' If inputdf contained bam, bwa alignment will be escaped.
#' Additional columns: group (group id), sample (sample id), PL (platform, i.e. illumina),
#' LB (library id), PU (unit, i.e. unit1). These strings (or info) will pass to BWA mem through -R.
#'
#' @param runbwa Set up BWA-mem, default=TRUE.
#' @param markDup Mark Duplicates, default=TRUE.
#' @param addRG Add or replace Read Groups using Picard AddOrReplaceReadGroups, default=FALSE.
#' @param rungatk Setup GATK, default=FALSE.
#'
#' @param ref.fa The full path of genome with bwa indexed reference fasta file.
#' @param gatkpwd The absolute path of GenomeAnalysisTK.jar.
#' @param picardpwd The absolute path of picard.jar.
#' @param minscore Minimum score to output, default=5, [bwa 30]. It will pass to bwa mem -T INT.
#'
#' @param realignInDels Realign Indels, default=FALSE. IF TRUE, a golden indel.vcf file should be provided.
#' @param indels.vcf The full path of indels.vcf.
#' @param recalBases Recalibrate Bases, default=FALSE. IF TRUE, a golden snps.vcf file should be provided.
#' @param dbsnp.vcf The full path of dbsnp.vcf.
#'
#' @param shbase Base for the shell id, i.e. "slurm-script/run_gatk_". [chr]
#' @param jobid Job ID, default="runarray". [chr]
#' @param email Your email address that farm will email to once the jobs were done/failed.
#' @param runinfo Parameters specify the array job partition information.
#' A vector of c(FALSE, "bigmemh", "1"): 1) run or not, default=FALSE
#' 2) -p partition name, default=bigmemh and 3) --cpus, default=1. 4) mem, default=1.5, in Gb.
#' It will pass to \code{set_array_job}.
#'
#' @return return a batch of shell scripts.
#'
#' @examples
#' inputdf <- data.frame(fq1="fq_1.fq", fq2="f1_2.fq", out="mysample",
#'                  group="g1", sample="s1", PL="illumina", LB="lib1", PU="unit1")
#'
#' run_GATK(inputdf, runbwa=TRUE, markDup=TRUE, addRG=FALSE,rungatk=FALSE,
#'          ref.fa="~/dbcenter/Ecoli/reference/Ecoli_k12_MG1655.fasta",
#'          gatkpwd="$HOME/bin/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar",
#'          picardpwd="$HOME/bin/picard-tools-2.1.1/picard.jar",
#'          minscore=5,
#'          realignInDels=FALSE, indels.vcf="indels.vcf",
#'          recalBases=FALSE, dbsnp.vcf="dbsnp.vcf",
#'          shbase=NULL, jobid="runarray",
#'          email=NULL, runinfo = c(FALSE, "batch", 1, "1.5", "10:00:00"))
#'
#' @export
run_GATK <- function(inputdf, runbwa=TRUE, markDup=TRUE, addRG=FALSE, rungatk=FALSE,

                     ref.fa="~/dbcenter/Ecoli/reference/Ecoli_k12_MG1655.fasta",
                     gatkpwd="$HOME/bin/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar",
                     picardpwd="$HOME/bin/picard-tools-2.1.1/picard.jar",
                     minscore=5,
                     realignInDels=FALSE, indels.vcf="indels.vcf",
                     recalBases=FALSE, dbsnp.vcf="dbsnp.vcf",

                     shbase=NULL, jobid="runarray",
                     email=NULL, runinfo = c(FALSE, "batch", 1, "1.5", "10:00:00")){

  ##### prepare parameters:
  fq <- inputdf
  ### determine memory based on partition
  #run <- get_runinfo(runinfo)

  #### create dir if not exist
  dir.create("slurm-script", showWarnings = FALSE)
  for(i in 1:nrow(fq)){

    if(is.null(shbase)){
        shid <- paste0("slurm-script/run_gatk_", i, ".sh")
    }else{
        shid <- paste0(shbase, i, ".sh")
    }

    ### header of the shell code
    cat("### GATK pipeline created by huskeR",
        paste("###", format(Sys.time(), "%a %b %d %X %Y")),
        paste(""),
        file=shid, sep="\n", append=FALSE)

    if(runbwa){
        ### alignment and sorting using picard-tools
        inputbam <- set_bwa(fq, runinfo, minscore, picardpwd, i, ref.fa, shid)
    }else if(sum(names(fq) %in% "bam") > 0){
        inputbam <- fq$bam[i]
    }else{
        stop("### no bam files, either set runbwa=TRUE, or add a `bam` column on inputdf.")
    }

    #### mark duplicates
    if(markDup) inputbam <- set_markDup(fq, picardpwd, inputbam, i, runinfo, shid)

    #### Add read groups, sort
    if(addRG) inputbam <- set_addRG(fq, picardpwd, inputbam, i, runinfo, shid)

    ### Perform local realignment around indels
    if(realignInDels) inputbam <- set_realignInDels(fq, inputbam, i, indels.vcf, ref.fa, gatkpwd, runinfo, shid)

    ### Recalibrate Bases
    if(recalBases) inputbam <- set_recalBases(fq, inputbam, i, indels.vcf, dbsnp.vcf, ref.fa, gatkpwd, runinfo, shid)

    ### Variant Discovery using HaplotypeCaller
    if(rungatk) vcaller(fq, inputbam, i, ref.fa, gatkpwd, runinfo, shid)
  }

  if(is.null(shbase)){
      shid2 <- paste0("slurm-script/run_gatk_$SLURM_ARRAY_TASK_ID.sh")
  }else{
      shid2 <- paste0(shbase, "$SLURM_ARRAY_TASK_ID.sh")
  }
  shcode <- paste0("module load java", "\n",
                   "module load bwa", "\n",
                   "module load samtools", "\n",
                   "sh ", shid2)

  if(is.null(shbase)){
      ash <- paste0("slurm-script/run_gatk_array.sh")
  }else{
      ash <- paste0(shbase, "array.sh")
  }

  set_array_job(shid= ash,
                shcode=shcode, arrayjobs=paste("1", nrow(inputdf), sep="-"),
                wd=NULL, jobid=jobid, email=email, runinfo=runinfo)
  #  sbatch -p bigmemh --mem 32784 --ntasks=4  slurm-script/run_gatk_array.sh
}


##########
set_bwa <- function(fq, run, minscore, picardpwd, i, ref.fa, shid){
    #Generate a SAM file containing aligned reads
    #http://gatkforums.broadinstitute.org/gatk/discussion/2799/howto-map-and-mark-duplicates
    rg <- paste0("\'@RG\\tID:", fq$group[i], "\\tSM:", fq$sample[i],
                 "\\tPL:", fq$PL[i], "\\tLB:", fq$LB[i], "\\tPU:", fq$PU[i], "\'")

    aligned_sam <- paste0(fq$out[i], ".aln.sam")
    sorted_bam <- paste0(fq$out[i], ".sorted.bam")
    #Generate a SAM file containing aligned reads
    #http://gatkforums.broadinstitute.org/gatk/discussion/2799/howto-map-and-mark-duplicates
    cat(paste("### Generate a SAM file containing aligned reads"),
        paste("bwa mem -M -R", rg, "-T", minscore, "-t", run[3], ref.fa, fq$fq1[i], fq$fq2[i], ">", aligned_sam),
        paste(""),

        ### http://broadinstitute.github.io/picard/
        paste0("java -Xmx", floor(as.numeric(run[4])/1.024),
               "g -jar ", picardpwd, " SortSam \\"),
        paste0("    INPUT=", aligned_sam, " \\"),
        paste0("    OUTPUT=", sorted_bam, " \\"),
        "    SORT_ORDER=coordinate",
        paste("#rm", aligned_sam),
        paste(""),
        file=shid, sep="\n", append=TRUE)
    message("###>>> set up BWA mem and then sort to bam using picard-tools!")
    return(sorted_bam)
}


##########
set_markDup <- function(fq, picardpwd, inputbam, i, run, shid){
  ### http://broadinstitute.github.io/picard/
  sorted_bam <- inputbam
  dedup_bam <- paste0(fq$out[i], ".sorted.dedup.bam")
  metrics <- paste0(fq$out[i], "_metrics.txt")
  log <- paste0(fq$out[i], ".sorted.bam.log")

  cat("### Mark duplicates",
      paste0("java -Xmx", floor(as.numeric(run[4])/1.024), "g ",
             "-jar ", picardpwd, " MarkDuplicates \\"),
      paste0("INPUT=", sorted_bam, " \\"),
      paste0("OUTPUT=", dedup_bam, " \\"),
      paste0("VALIDATION_STRINGENCY=SILENT \\"),
      paste0("METRICS_FILE=", metrics),
      #paste(""),
      paste0(""),
      paste0("java -Xmx", floor(as.numeric(run[4])/1.024), "g ",
             "-jar $HOME/bin/picard-tools-2.1.1/picard.jar BuildBamIndex \\"),
      paste0("INPUT=", dedup_bam),
      "",
      paste0("samtools flagstat ", sorted_bam, " > ", log ),
      "",
      #paste("rm", sorted_bam),
      file=shid, sep="\n", append=TRUE)
  message("###>>> set up Mark Duplicates using picard-tools!")
  return(dedup_bam)
}

###
set_addRG <- function(fq, picardpwd, inputbam, i, run, shid){
  ### http://broadinstitute.github.io/picard/
  dedupRG_bam <- gsub("bam$", "RG.bam", inputbam)
  #rg <- paste0("\'@RG\\tID:", fq$group[i], "\\tSM:", fq$sample[i],
  #             "\\tPL:", fq$PL[i], "\\tLB:", fq$LB[i], "\\tPU:", fq$PU[i], "\'")

  cat("### add read groups",
      paste0("java -Xmx", floor(as.numeric(run[4])/1.024), "g ",
             "-jar ", picardpwd, " AddOrReplaceReadGroups \\"),
      paste0("INPUT=", inputbam, " \\"),
      paste0("OUTPUT=", dedupRG_bam, " \\"),
      paste0("SORT_ORDER=coordinate \\"),
      paste0("CREATE_INDEX=true \\"),
      paste0("RGID=", fq$group[i], " \\"),
      paste0("RGLB=", fq$LB[i], " \\"),
      paste0("RGPL=", fq$PL[i], " \\"),
      paste0("RGPU=", fq$PU[i], " \\"),
      paste0("RGSM=", fq$sample[i]),
      #paste("RGPL=illumina \"),
      paste0(""),
      "",
      #paste("rm", sorted_bam),
      file=shid, sep="\n", append=TRUE)
  message("###>>> set up add RG using picard-tools!")
  return(dedupRG_bam)
}


##########
set_realignInDels <- function(fq, inputbam, i, indels.vcf, ref.fa, gatkpwd, run, shid){
  dir.create("$HOME/tmp", showWarnings = FALSE)

  ### input and output files
  realigned_bam <- gsub("bam", "indelrealigned.bam", inputbam)
  intervals <- paste0(fq$out[i], ".forIndelRealigner.intervals")

  cat("### Define intervals to target for local realignment",
      paste0("java -Xmx", floor(as.numeric(run[4])/1.024), "g ", "-jar ", gatkpwd, " \\"),
      paste0("-T RealignerTargetCreator \\"),
      paste0("-R ", ref.fa, " \\"),
      paste0("-I ", inputbam, " \\"),
      paste0("--known ", indels.vcf, " \\"),
      paste0("-o ", intervals),
      paste(""),
      file=shid, sep="\n", append=TRUE)

  cat("### Perform local realignment around indels",
      paste0("### link: https://www.broadinstitute.org/gatk/guide/article?id=2800"),
      paste0("java -Xmx", floor(as.numeric(run[4])/1.024), "g -Djava.io.tmpdir=$HOME/tmp \\"),
      paste0("-jar ", gatkpwd, " \\"),
      paste0("-I ", inputbam, " \\"),
      paste0("-R ", ref.fa, " \\"),
      paste0("-T IndelRealigner \\"),
      paste0("-targetIntervals ", intervals, " \\"),
      paste0("-o ", realigned_bam, " \\"),
      paste0("-known ", indels.vcf),
      paste0("--consensusDeterminationModel KNOWNS_ONLY \\"),
      paste0("-LOD 0.4"),
      paste(""),
      file=shid, sep="\n", append=TRUE)
  message("###>>> set up Realign InDels using GATK!")
  return(realigned_bam)

}

##########
set_recalBases <- function(fq, inputbam, i, indels.vcf, dbsnp.vcf, ref.fa, gatkpwd, run, shid){

  ### input and output files
  recal_table <- paste0(fq$out[i], ".recal_data.table")
  post_recal_table <- paste0(fq$out[i], ".post_recal_data.table")
  plotpdf <- paste0(fq$out[i], ".recalibration_plots.pdf")
  recal_bam <- gsub("bam$", "recal.bam", inputbam)

  #1. Analyze patterns of covariation in the sequence dataset
  cat("### Recalibrate base quality scores = run BQSR",
      "### link: https://www.broadinstitute.org/gatk/guide/article?id=2801",
      paste0("java -Xmx", floor(as.numeric(run[4])/1.024), "g ", "-jar ", gatkpwd, " \\"),
      paste0("-T BaseRecalibrator \\"),
      paste0("-R ", ref.fa, "\\"),
      paste0("-I ", inputbam, " \\"),
      paste0("-knownSites ", dbsnp.vcf, " \\"),
      #paste0("-knownSites ", indels.vcf, " \\"),
      paste0("-o ", recal_table),
      "",
      file=shid, sep="\n", append=TRUE)

  #2. Do a second pass to analyze covariation remaining after recalibration
  cat(paste0("java -Xmx", floor(as.numeric(run[4])/1.024), "g ", "-jar ", gatkpwd, " \\"),
      paste0("-T BaseRecalibrator \\"),
      paste0("-R ", ref.fa, " \\"),
      paste0("-I ", inputbam, " \\"),
      paste0("-knownSites ", dbsnp.vcf, " \\"),
      paste0("-knownSites ", indels.vcf, " \\"),
      paste0("-BQSR ", recal_table, " \\"),
      paste0("-o ", post_recal_table),
      "",
      file=shid, sep="\n", append=TRUE)

  #3. Generate before/after plots
  cat(paste0("java -Xmx", floor(as.numeric(run[4])/1.024), "g ", "-jar ", gatkpwd, " \\"),
      paste0("-T AnalyzeCovariates \\"),
      paste0("-R ", ref.fa, " \\"),
      paste0("-before ", recal_table, " \\"),
      paste0("-after ", post_recal_table, " \\"),
      paste0("-plots ", plotpdf),
      "",
      file=shid, sep="\n", append=TRUE)

  #4. Apply the recalibration to your sequence data
  cat(paste0("java -Xmx", floor(as.numeric(run[4])/1.024), "g ", "-jar ", gatkpwd, " \\"),
      paste0("-T PrintReads \\"),
      paste0("-R ", ref.fa, " \\"),
      paste0("-I ", inputbam, " \\"),
      paste0("-BQSR ", recal_table, " \\"),
      paste0("-o ", recal_bam),
      "",
      file=shid, sep="\n", append=TRUE)
  message("###>>> set up Recalibrate Bases using GATK!")
  return(recal_bam)
}

vcaller <- function(fq, inputbam, i, ref.fa, gatkpwd, run, shid){
  ### input and output files
  gvcf <- paste0(fq$out[i], ".g.vcf")

  cat("### Call variants in your sequence data",
      paste0("java -Xmx", floor(as.numeric(run[4])/1.024), "g ", "-jar ", gatkpwd, " \\"),
      paste0("-T HaplotypeCaller \\"),
      paste0("-R ", ref.fa, " \\"),
      paste0("-I ", inputbam, " \\"),
      #paste0("--genotyping_mode DISCOVERY \\"),
      "-ERC GVCF \\",
      #paste0("-stand_emit_conf 10 \\"),
      #paste0("-stand_call_conf 30 \\"),
      paste0("-o ", gvcf),
      file=shid, sep="\n", append=TRUE)
  message("###>>> set up Variants calling using GATK HaplotypeCaller!")

}





