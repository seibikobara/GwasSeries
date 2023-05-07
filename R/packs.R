#' @title 
#' Gwas catalog manipulation
#' @description 
#' return grange
#' @aliases gwas_grange
#' @param csvFilepath path of csv file
#' @param gap the gap bps across the position of SNP
#' @export
#' @import tidyverse GenomicRanges


gwas_grange <- function(csvFilePath, gap){
    path = file.path(csvFilePath)
    gwas = readr::read_csv(path)
    chr_loc = stringr::str_split(gwas$Location, pattern = stringr::regex(":"))
    temp = lapply(chr_loc, as.numeric)
    temp1 = as.data.frame(do.call(rbind, temp))
    names(temp1) = c("chr","pos") 
    # rmeove na
    temp1 = temp1[complete.cases(temp1),]

    gr = GenomicRanges::GRanges(seqnames = temp1$chr,
            range = IRanges::IRanges(start = temp1$pos-gap, end = temp1$pos+gap) 
     )
     return(gr)
}



#' @title 
#' Chromatin states manipulation
#' @description 
#' return grange
#' @aliases chrome_state_grange
#' @param csvFilepath path of csv file
#' @param chr chrmosome (string) ex, chr1, chr2...
#' @export
#' @import tidyverse GenomicRanges



# chromatin state to grange
chrom_state_grange = function(csvFilePath,chr){
    path = file.path(csvFilePath)
    temp = readr::read_csv(path)
    states = GenomicRanges::GRanges(seqnames = IRanges::Rle(chr, nrow(temp)),
                range = IRanges::IRanges(start = temp$start, end = temp$end) 
     )
    # add meta (states)
    values(states) = data.frame(states = temp$state)
    return(states)
}




