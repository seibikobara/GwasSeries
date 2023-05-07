#' @title 
#' Gwas catalog manipulation
#' @description 
#' return grange
#' @aliases gwas_grange
#' @param csvFilepath path of csv file
#' @param gap the gap bps across the position of SNP
#' @export
#' @import tidyverse magrittr GenomicFeature GenomicRanges janitor


gwas_grange <- function(csvFilePath, gap){
    path = file.path(csvFilePath)
    gwas = read_csv(path)
    chr_loc = str_split(gwas$Location, pattern = regex(":"))
    temp = lapply(chr_loc, as.numeric)
    temp1 = as.data.frame(do.call(rbind, temp))
    names(temp1) = c("chr","pos") 
    # rmeove na
    temp1 = temp1[complete.cases(temp1),]

    gr = GRanges(seqnames = temp1$chr,
            range = IRanges(start = temp1$pos-gap, end = temp1$pos+gap) 
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
#' @import tidyverse magrittr GenomicFeature GenomicRanges janitor



# chromatin state to grange
chrom_state_grange = function(csvFilePath,chr){
    path = file.path(csvFilePath)
    temp = read_csv(path)
    states = GRanges(seqnames = Rle(chr, nrow(temp)),
                range = IRanges(start = temp$start, end = temp$end) 
     )
    # add meta (states)
    values(states) = data.frame(states = temp$state)
    return(states)
}




