### load gbs file, change to tfile format for PLINK, filtering multiple alleles and 
### then change to "A A" coding and change IUPAC Ambiguity Codes
gbs2plink <- function(gbsfile="largedata/genotypes/SB.imputed.hmp",
                    outfile="~/dbcenter/seeds_data/chr10_filetered_unimputed.bed", ...){
    ### GBS  
  
    ### read in GBS file
    gbs <- fread(gbsfile, data.table=FALSE, ...)
    message(sprintf("[gbs2plink]: Loaded [ %s ] SNPs and [ %s ] cols for file [%s]!", nrow(gbs), ncol(gbs), gbsfile))
    
    ### sanity check
    sc <- sum(names(gbs)[1:11] == c("rs.", "alleles", "chrom", "pos", "strand", 
                                    "assembly.", "center",
                                    "protLSID", "assayLSID", "panelLSID", "QCcode"))
    if(sc != 11){stop("[gbs2plink]: Your files seems not GBS hmp file!")}
    
    
    ### change to BED5+ format
    tped <- gbs[, c(3,1,4,4,2,5, 12:ncol(gbs))]
    names(tped)[1:6] <- c("chr", "snpid", "genetic", "pos", "alleles", "nchar")
    nms <- names(tped)
    nms <- gsub("\\..*$", "", nms)
    names(tped) <- nms
    tped$genetic <- 0
    message(sprintf("[gbs2plink]: Changed to BED5+ format and start filtering ..."))
    
    ### filter SNPs contain multiple alleles
    tped$nchar <- nchar(as.character(tped$alleles))
    tped <- subset(tped, nchar == 3)
    tped <- tped[, -5:-6]
    #idx <- grep("-", subg$alleles)
    #subg <- subg[-idx,]
    message(sprintf("[gbs2plink]: Remaining [ %s ] sites with two variations!", nrow(tped)))
    
    ### TMAP
    tmap <- data.frame(fid= names(tped)[-1:-4], iid=names(tped)[-1:-4], father=0, mother=0, sex=0, pheno=-9)
    
    
    message(sprintf("[gbs2plink]: Start to transforming, recoding and writing ..."))
    ### change to two identical haplotypes
    for(i in 5:ncol(tped)){
        tped[, i] <- paste( gsub(".$", "", tped[,i] ),  gsub("^.", "", tped[,i] ))
        #print(i)
    }
    
    out.tped <- paste(outfile, ".tped", sep="")
    write.table(tped, out.tped, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
    
    out.tmap <- paste(outfile, ".tfam", sep="")
    write.table(tmap, out.tmap, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
    
    message(sprintf("[gbs2plink]: DONE!"))
}

iupac <- function(subg){
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
  return(subg)
}