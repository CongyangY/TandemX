options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 2, timeout = 600)
cran <- c("BiocManager", "remotes", "stringr", "stringdist", "circlize", "seqinr", "doParallel", "getopt", "ape")
install.packages(cran, destdir = "/opt/source_archives")
BiocManager::install(c("Biostrings", "pwalign", "msa"), version = "3.20",
                     ask = FALSE, update = FALSE, destdir = "/opt/source_archives")
required <- c(cran, "Biostrings", "pwalign", "msa")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop(paste("Missing comparator dependencies:", paste(missing, collapse = ", ")))
packages <- installed.packages()
write.table(packages[, c("Package", "Version", "Built", "LibPath")],
            "/opt/provenance/R_packages.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), "/opt/provenance/R_sessionInfo.txt")
