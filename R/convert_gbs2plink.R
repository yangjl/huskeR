#' \code{Change gbs hmp to plink format}
#'
#' load gbs file, change to tfile format for PLINK, filtering multiple alleles and
#' then change to "A A" coding and change IUPAC Ambiguity Codes
#'
#' @param gbsfile An input gbs in hmp format. Note without # lines. [chr, ="g2f_2014_zeagbsv27v5hmp_nodash.txt"]
#' @param gbs If gbsfile is NULL load this gbs object. [obj, data.table=NULL]
#' @param iupac Change the IUPAC Code. [logical, =TRUE]
#' @param mode GBS code mode. mode=1, A, T, C, G; mode=2, AA, TT, CC, GG. [num, =1]
#' @param outfile Output file name. [chr, ="~/dbcenter/seeds_data/chr10_filetered_unimputed.bed"]
#' @param jmph
#'
#' plink --tped test.tped --tfam test.tfam --make-bed --missing-genotype N --allow-no-sex --make-founders --out ../Data/sim
#' @export
convert_gbs2plink <- function(gbsfile="g2f_2014_zeagbsv27v5hmp_nodash.txt",
                              gbs=NULL, iupac=TRUE, mode=1,
                    outfile="~/dbcenter/seeds_data/chr10_filetered_unimputed.bed"){
    ### GBS
    #library("data.table")

    ### read in GBS file
    if(!is.null(gbsfile)){
        gbs <- fread(gbsfile, data.table=FALSE)
    }

    message(sprintf("[gbs2plink]: Loaded [ %s ] SNPs and [ %s ] cols for file [%s]!", nrow(gbs), ncol(gbs), gbsfile))
    ### sanity check
    sc <- sum(names(gbs)[1:11] == c("rs#", "alleles", "chrom", "pos", "strand",
                                    "assembly#", "center",
                                    "protLSID", "assayLSID", "panelLSID", "QCcode"))
    if(sc != 11){stop("[gbs2plink]: Your files seems not a GBS hmp file!")}


    ### change to BED5+ format
    tped <- gbs[, c(3,1,4,4,2,5, 12:ncol(gbs))]
    names(tped)[1:6] <- c("chr", "snpid", "genetic", "pos", "alleles", "nchar")
    nms <- names(tped)
    nms <- gsub("\\..*$", "", nms)
    names(tped) <- nms
    tped$genetic <- 0
    message(sprintf("[gbs2plink]: Changed to PLINK (tped) format and start filtering ..."))

    ### filter SNPs contain multiple alleles
    tped$nchar <- nchar(as.character(tped$alleles))
    tped <- subset(tped, nchar == 3)

    ### remove "indel"
    idx <- grep("-", tped$alleles)

    tped <- tped[-idx, -5:-6]
    #idx <- grep("-", subg$alleles)
    #subg <- subg[-idx,]
    message(sprintf("[gbs2plink]: Remaining [ %s ] sites with two variants!", nrow(tped)))

    ### TMAP
    tmap <- data.frame(fid= names(tped)[-1:-4], iid=names(tped)[-1:-4], father=0, mother=0, sex=0, pheno=-9)


    message(sprintf("[gbs2plink]: Start to transforming and recoding ..."))
    ### change to two identical haplotypes
    for(i in 5:ncol(tped)){
        if(mode == 1){
            tped[, i] <- paste(tped[,i], tped[,i])
        }
        if(mode == 2){
            tped[, i] <- paste( gsub(".$", "", tped[,i] ),  gsub("^.", "", tped[,i] ))
        }
        #print(i)
    }

    if(iupac){
        message(sprintf("[gbs2plink]: running IUPAC ..."))
        tped <- run_iupac(tped)
    }

    message(sprintf("[gbs2plink]: Start to writing tped and tfam ..."))
    out.tped <- paste(outfile, ".tped", sep="")
    fwrite(tped, out.tped, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

    out.tmap <- paste(outfile, ".tfam", sep="")
    fwrite(tmap, out.tmap, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

    message(sprintf("[gbs2plink]: DONE!"))
}

#' @rdname convert_gbs2plink
run_iupac <- function(subg){
  ###change IUPAC Ambiguity Codes
  #M    A or C    K
  #R	A or G	Y
  #W	A or T	W
  #S	C or G	S
  #Y	C or T	R
  #K	G or T	M
  subg[subg=="M M"] <- "A C"
  subg[subg=="R R"] <- "A G"
  subg[subg=="W W"] <- "A T"
  subg[subg=="S S"] <- "C G"
  subg[subg=="Y Y"] <- "C T"
  subg[subg=="K K"] <- "G T"
  subg[subg=="0 0"] <- "N N"
  #subg[subg=="- -"] <- "N N"
  return(subg)
}
