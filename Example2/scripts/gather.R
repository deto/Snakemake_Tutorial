samples <- snakemake@input
out_file <- snakemake@output[["expression"]]

data <- lapply(samples, function(s){

    sname <- strsplit(s, "/")[[1]][1]

    sample_data <- read.table(s, sep="\t", row.names=1)
    colnames(sample_data) <- sname
    return(sample_data)
})

data <- do.call(cbind, data)

write.table(data, out_file, sep="\t", col.names=NA)
