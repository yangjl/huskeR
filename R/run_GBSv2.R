#' \code{Run TASSEL GBSv2 on slurm-based HPC}
#'
#' a single slurm job.
#' help doc: https://bitbucket.org/tasseladmin/tassel-5-source/wiki/Tassel5GBSv2Pipeline
#'
#' @parm outdir Out put file directory. A batch of subdir will be created under it. [chr, ="largedata/gbs"]
#' @param shfile A shell script in "slurm-script/". [chr, "slurm-script/run_it.sh"]
#' @param fqdir Absolute dir for the fastq files, can be nested within the dir. [chr, "my/fq/filedir"]
#' @param keyfile The key file. The GBSv2 key file has 4 required headers: Flowcell, Lane, Barcode and FullSampleName. [chr, "data/file.txt]
#' @param bt2_idx The idx file for bowtie2. [chr, "/lustre/work/jyanglab/jyang21/dbcenter/AGP/AGPv4/AGPv4_bt2_index/AGPv4_bt2_index"]
#' @param sam The output sam file. [chr, "tagexport/tagsForAlignFullvs_F1s.sam"]
#' @param ref Fastq file of the ref. [chr, "/lustre/work/jyanglab/jyang21/dbcenter/AGP/AGPv4/AGPv4_bt2_index/AGPv4_bt2.fa"]
#' @param snpqc_file SNP qc file. [chr, "largedata/gbs/outputStats.txt"]
#'
#' @param db Relative file name of the db. [chr, =NULL, or "data/gbs.db"]
#' @param mem Memory usage. [num, 50]
#' @param cpu Number of CPUs to use. [num, 50]
#' @param kmerlen Kmer length [num, =64]
#' @param enzyme Enzyme used for cutting. [chr, ="ApeKI"]
#' @param mnlcov mnLCov. [num, =0.05]
#' @param mnmaf Min MAF. [num, =0.001]
#' @param deleteolddata Delete Old Data. [chr, ="true"]
#'
#' @param production SNP production. [logical, =F]
#' @param seq2tag Generated tags from fastq files. [logical, =T]
#' @param tag2fq extract tags from db for alignment to reference. [logical, =T]
#' @param bt2 Run bt2 alignment. [logical, =T]
#' @param snpcall Call SNP. [logical, =T]
#'
#'
#' @return return a batch of shell scripts.
#'
#' @examples
#'
#'
#' @export
run_GBSv2 <- function(outdir="largedata/gbs", shfile, fqdir, keyfile,
                      bt2_idx, sam, ref, snpqc_file, h5,
                      db=NULL, mem, cpu, kmerlen=64, enzyme="ApeKI",
                      mnlcov=0.05, mnmaf=0.01, deleteolddata="true",
                      production=FALSE, seq2tag=TRUE, tag2fq=TRUE, bt2=TRUE, snpcall=TRUE){
    #### create dir if not exist
    dir.create("slurm-script", showWarnings = FALSE)

    #dir.create(paste(outdir, "tagexport", sep="/"), showWarnings = FALSE)
    #dir.create(paste(outdir, "discovery", sep="/"), showWarnings = FALSE)
    #dir.create(paste(outdir, "db", sep="/"), showWarnings = FALSE)

    #if(is.null(db)){
    #    db = paste(outdir, "db", sep="/")
    #    dir.create(db, showWarnings = FALSE)
    #}

    # 1. load modules
    cat(paste("## [run_GBSv2]: Tassel5.2 GBSv2", Sys.time(), sep=" "),
        paste0("module load java/1.8 tassel/5.2 bowtie/2.2"),
        "",
        file=shfile, sep="\n", append=FALSE)

    # 2. generate tags from fastq file
    if(production){
        ProductionSNPCaller(shfile, mem, db, enzyme, fqdir, keyfile, kmerlen, h5)
    }else{
        if(seq2tag){
            GBSSeqToTag(shfile, mem, fqdir, db, keyfile, kmerlen, enzyme)
        }
        tagexport = paste(outdir, "tagexport.fa.gz", sep="/")
        if(tag2fq){
            TagExportToFastq(shfile, mem, db, tagexport, mindepth=1)
        }
        if(bt2){
            run_bowtie2(shfile, cpu, bt2_idx, tagexport, sam)
            SAMToGBSdb(shfile, mem, sam, db)
        }
        if(snpcall){
            DiscoverySNP(shfile, mem, db, ref, mnlcov=0.05, mnmaf=0.001, deleteolddata)
            SNPQuality(shfile, mem, db, snpqc_file)
        }
    }

}


#' \code{generate tags from fastq file}
#'
#'
#' @rdname run_GBSv2
GBSSeqToTag <- function(shfile, mem, fqdir, db, keyfile, kmerlen, enzyme){

  cat(paste("## [run_GBSv2]: set GBSSeqToTag, Tassel5.2 GBSv2"),

      paste0("run_pipeline.pl -Xmx", mem, "g -fork1 -GBSSeqToTagDBPlugin -e ", enzyme, " \\"),
      paste("-i", fqdir, " \\"),
      paste("-db", db, " \\"),
      paste("-k", keyfile, " \\"),

      "-c 10 \\",
      paste("-kmerLength", kmerlen, "\\"),
      "-minKmerL 20 \\",
      "-mnQS 20 \\",
      "-mxKmerNum 100000000 \\",
      "-batchSize 8 \\",
      "-endPlugin -runfork1",
      "",

      file=shfile, sep="\n", append=TRUE)
}


#' \code{extract tags from db for alignment to reference}
#'
#' @param tagexport Fastq.gz file to export tags. [chr, "tagexport/tagsForAlign.fa.gz"]
#' @param mindepth The minimum count of reads for a tag to be output. [num, =1]
#'
#'
#' @rdname run_GBSv2
TagExportToFastq <- function(shfile, mem, db, tagexport, mindepth=1){

    cat(paste("## [run_GBSv2]: set TagExportToFastq, Tassel5.2 GBSv2"),

        paste0("run_pipeline.pl -Xmx", mem, "g -fork1 -TagExportToFastqPlugin \\"),
        paste("-db", db, " \\"),
        paste("-o", tagexport, " \\"),
        paste("-c", mindepth, " \\"),
        "-endPlugin -runfork1",
        "",

        file=shfile, sep="\n", append=TRUE)
}

#' \code{map tags to reference use bt2}
#' #--very-sensitive-local -D 20 -R 3 -N 0 -L 20 -i S,1,0.50
#' pigz -k -d -p 50 Zea_mays.AGPv4.dna.toplevel.fa.gz
#' bowtie2-build Zea_mays.AGPv4.dna.toplevel.fa AGPv4_bt2_index --threads 50
#'
#'
#'
#' @rdname run_GBSv2
run_bowtie2 <- function(shfile, cpu, bt2_idx, tagexport, sam){

    cat(paste("## [run_GBSv2]: set bowtie2, Tassel5.2 GBSv2"),

        paste0("bowtie2 -p ", cpu, " \\"),
        "--very-sensitive-local \\",
        paste("-x", bt2_idx, " \\"),
        paste("-U", tagexport, " \\"),
        paste("-S", sam),
        "",
        file=shfile, sep="\n", append=TRUE)
}


#' \code{store tag physical map position in database}
#'
#'
#' @rdname run_GBSv2
SAMToGBSdb <- function(shfile, mem, sam, db){

    cat(paste("## [run_GBSv2]: set SAMToGBSdb, Tassel5.2 GBSv2"),

        paste0("run_pipeline.pl -Xmx", mem, "g -fork1 -SAMToGBSdbPlugin \\"),
        paste("-i", sam, " \\"),
        paste("-db", db, " \\"),

        paste("-aProp 0 \\"), #Minimum length of aligned base pair to store the SAM entry
        paste("-aLen 0 \\"), #Minimum proportion of sequence that must align to store the SAM entry
        "-endPlugin -runfork1",
        "",
        file=shfile, sep="\n", append=TRUE)
}


#' \code{run SNP discovery plugin}
#'
#'
#' @rdname run_GBSv2
DiscoverySNP <- function(shfile, mem, db, ref, mnlcov=0.05, mnmaf=0.001, deleteolddata ){

    cat(paste("## [run_GBSv2]: set DiscoverySNPCaller, Tassel5.2 GBSv2"),

        paste0("run_pipeline.pl -Xmx", mem, "g -fork1 -DiscoverySNPCallerPluginV2 \\"),
        paste("-db", db, " \\"),
        paste("-ref", ref, " \\"),

        paste("-mnLCov", mnlcov, " \\"),
        paste("-mnMAF", mnmaf, " \\"),
        paste("-deleteOldData", deleteolddata, " \\"),
        "-endPlugin -runfork1",
        "",
        file=shfile, sep="\n", append=TRUE)
}


#' \code{check SNP quality parameters}
#'
#'
#' @rdname run_GBSv2
SNPQuality <- function(shfile, mem, db, snpqc_file){

    cat(paste("## [run_GBSv2]: set SNPQualityProfilerPlugin, Tassel5.2 GBSv2"),

        paste0("run_pipeline.pl -Xmx", mem, "g -fork1 -SNPQualityProfilerPlugin \\"),
        paste("-db", db, " \\"),
        paste0("-statFile \'", snpqc_file, "\' \\"),

        "-endPlugin -runfork1",
        "",
        file=shfile, sep="\n", append=TRUE)
}



#' \code{run production pipeline}
#'
#' @param h5 Output hdf5 format. [chr, "gbs/discovery/output.h5"]
#'
#' @rdname run_GBSv2
ProductionSNPCaller <- function(shfile, mem, db, enzyme, fqdir, keyfile, kmerlen, h5){

    cat(paste("## [run_GBSv2]: set ProductionSNPCallerPluginV2, Tassel5.2 GBSv2"),

        paste0("run_pipeline.pl -Xmx", mem, "g -fork1 -ProductionSNPCallerPluginV2 \\"),
        paste("-db", db, " \\"),
        paste("-e", enzyme, " \\"),
        paste("-i", fqdir, " \\"),
        paste("-k", keyfile, " \\"),
        paste("-kmerLength", kmerlen, "\\"),
        paste("-o", h5, "\\"),

        "-endPlugin -runfork1",
        "",
        file=shfile, sep="\n", append=TRUE)
}










