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
    # add chr to chr numeric
    temp1 = temp1 |> dplyr::rowwise() |> dplyr::mutate(newchr = paste0("chr",chr))
    gr = GenomicRanges::GRanges(seqnames = temp1$newchr,
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
#' @export
#' @import tidyverse GenomicRanges



# chromatin state to grange
chrom_state_grange = function(csvFilePath){
    path = file.path(csvFilePath)
    temp = readr::read_csv(path)
    states = GenomicRanges::GRanges(seqnames = temp$chrom,
                range = IRanges::IRanges(start = temp$start, end = temp$end) 
     )
    # add meta (states)
    ## convert state to interpretable state name
    temp = temp |> dplyr::mutate(new_state=
                    ifelse(state == 1, "1: Active TSS",
                    ifelse(state == 2, "2: Flanking Active TSS",
                    ifelse(state == 3, "3: Transcr.at gene 5' and 3'",
                    ifelse(state == 4, "4: Strong transcription",
                    ifelse(state == 5, "5: Weak transcription",
                    ifelse(state == 6, "6: Genic enhancers",
                    ifelse(state == 7, "7: Enhancers",
                    ifelse(state == 8, "8: ZNF genes & repeats",
                    ifelse(state == 9, "9: Heterochromatin",
                    ifelse(state == 10, "Bivalent/Poised TSS",
                    ifelse(state == 11, "11: Flanking Bivalent TSS/Enh",
                    ifelse(state == 12, "12: Bivalent Enhancer",
                    ifelse(state == 13, "13: Repressed PolyComb",
                    ifelse(state == 14, "14: Weak Repressed PolyComb",
                    ifelse(state == 15, "15: Quiescent/Low", NA))))))))))))))))
    values(states) = data.frame(states = temp$new_state)
    return(states)
}





