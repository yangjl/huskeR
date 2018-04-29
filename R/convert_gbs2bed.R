### load gbs file, change to bed5+ format, filtering multiple alleles and 
### then change to "A A" coding and change IUPAC Ambiguity Codes
gbs2bed <- function(gbsfile="largedata/genotypes/SB.imputed.hmp",
                    outfile="~/dbcenter/seeds_data/chr10_filetered_unimputed.bed", ...){
    ### GBS hapmap file  
  
    ### read in GBS file
    gbs <- fread(gbsfile, data.table=FALSE, ...)
    message(sprintf("[gbs2bed]: Loaded [ %s ] SNPs and [ %s ] cols for file [%s]!", nrow(gbs), ncol(gbs), gbsfile))
    
    ### sanity check
    sc <- sum(names(gbs)[1:11] == c("rs.", "alleles", "chrom", "pos", "strand", 
                                    "assembly.", "center",
                                    "protLSID", "assayLSID", "panelLSID", "QCcode"))
    if(sc != 11){stop("[gbs2bed]: Your files seems not GBS hmp file!")}
    
    
    ### change to BED5+ format
    gbs <- gbs[, c(3,4,4,1,2,5, 12:ncol(gbs))]
    names(gbs)[1:6] <- c("chr", "start", "end", "snpid", "alleles", "nchar")
    nms <- names(gbs)
    nms <- gsub("\\..*$", "", nms)
    names(gbs) <- nms
    gbs$start <- gbs$start -1
    message(sprintf("[gbs2bed]: Changed to BED5+ format and start filtering ..."))
    
    ### filter SNPs contain multiple alleles
    gbs$nchar <- nchar(as.character(gbs$alleles))
    subg <- subset(gbs, nchar == 3)
    subg <- subg[, -6]
    #idx <- grep("-", subg$alleles)
    #subg <- subg[-idx,]
    message(sprintf("[gbs2bed]: Remaining [ %s ] sites with two variations!", nrow(subg)))
    
    message(sprintf("[gbs2bed]: Start to transforming, recoding and writing ..."))
    ### change to two identical haplotypes
    for(i in 6:ncol(subg)){
        subg[, i] <- paste(subg[, i], subg[, i])
        #print(i)
    }
    
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
    
    write.table(subg, outfile, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    message(sprintf("[gbs2bed]: DONE!"))
}
