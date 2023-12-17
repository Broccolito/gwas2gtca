library(dplyr)
library(purrr)
library(data.table)
library(GenomicRanges)
library(rtracklayer)

args = commandArgs(trailingOnly=TRUE)

filename = args[1]
chr = args[2]
pos = args[3]
build = args[4]
A1 = args[5]
A2 = args[6]
freq = args[7]
b = args[8]
se = args[9]
N = args[10]
output_name = args[11]

cat("Loading GWAS summary statistics...\n")
d = fread(filename)

cat("Parsing chromosome names...\n")
d[[chr]] = gsub(pattern = "chr", replacement = "", d[[chr]])
d[[chr]] = paste0("chr", d[[chr]])

cat("Generating GRange object...\n")
drange = GRanges(
  seqnames = d[[chr]],
  ranges = IRanges(start = d[[pos]], end = d[[pos]] + nchar(d[[A2]]) - 1),
  strand = "*",
  A1 = d[[A1]],
  A2 = d[[A2]],
  freq = d[[freq]],
  b = d[[b]],
  se = d[[se]],
  N = d[[N]]
)

cat("Liftover to HG38 if needed...\n")
if(build != "HG38"){
  ch = import.chain("hg19ToHg38.over.chain")
  drange = liftOver(drange, ch)
}

cat("Converting to dataframe...\n")
d_hg38 = as.data.frame(drange)

cat("Generating standard SNP names...\n")
SNP = paste(gsub(pattern = "chr", replacement = "", d_hg38[["seqnames"]]), 
            d_hg38[["start"]], d_hg38[["A2"]], d_hg38[["A1"]], sep = ":")

cat("Generating .ma result file and filtering out INDELs...\n")
d_ma = data.frame(
  chr = d_hg38[["seqnames"]],
  SNP,
  A1 = d_hg38[["A1"]],
  A2 = d_hg38[["A2"]],
  freq = d_hg38[["freq"]],
  b = d_hg38[["b"]],
  se = d_hg38[["se"]],
  N = d_hg38[["N"]]
) %>%
  filter(nchar(A1) == 1) %>%
  filter(nchar(A2) == 1)

null_return = d_ma %>%
  split(.$chr) %>%
  map(function(x){
    chromosome = x$chr[1]
    cat(paste0("Writing .ma file for ", chromosome, "...\n"))
    x = x[,-1]
    fwrite(x, file = paste0(output_name, "_", chromosome, ".ma"), sep = " ")
  })

cat("GWAS formatting completed...\n\n")




