# Read vcf file
my.read.vcf <- function(file, special.char="##", ...) {
    # Making a search term that looks like: "##.*", tells R 
    # to find anything containing the pattern "##" followed 
    # by anything (* is wildcard)   
    my.search.term <- paste0(special.char, ".*")  
    # Replace any line containing the search term with nothing (in other words remove it)
    clean.lines <- sub(my.search.term, "", readLines(file)) 
    # Replace the #CHROM term in the header with CHROM, so 
    # R doesn't treat it as a special character   
    clean.lines2 <- sub("#CHROM", "CHROM", clean.lines) 
    # Pass the cleaned up lines to read.table
    read.table(..., text=paste(clean.lines2, collapse="\n")) 
}

GBS=my.read.vcf(file= "/01_Data/SAP_GS.vcf", header=TRUE, stringsAsFactors = TRUE, as.is=TRUE)

# Function Parse vcf file to convert to -1,0,1 format
parse.GBS <- function(x) {
    unique.x <- unique(x)
    alleles <- setdiff(unique.x,union("H","N"))
    y <- rep(0,length(x))
    y[which(x==alleles[1])] <- -1
    y[which(x==alleles[2])] <- 1
    y[which(x=="N")] <- NA
    return(y)
}

X <- apply(GBS[, -c(1:f.column)],1,parse.GBS) # genotype matrix
