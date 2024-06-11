setwd(here::here())
setwd("doc/paper-jss/")

system('R -e \'knitr::knit("sdmTMB-paper.Rnw");tinytex::latexmk("sdmTMB-paper.tex", engine = "xelatex", clean = FALSE)\'')

system("gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dFirstPage=51 -dLastPage=54 -sOUTPUTFILE=reviewer-response2.pdf sdmTMB-paper.pdf")

d <- readLines("sdmTMB-paper.Rnw")
i <- grep("response2}", d)
d[seq(i-1, i+1)] <- paste0("% ", d[seq(i-1, i+1)])
writeLines(d, "sdmTMB-paper.Rnw")

system('R -e \'knitr::knit("sdmTMB-paper.Rnw");tinytex::latexmk("sdmTMB-paper.tex", engine = "xelatex", clean = FALSE)\'')

system("mv sdmTMB-paper.pdf sdmTMB-paper-V3-changes.pdf")

d <- readLines("sdmTMB-paper.Rnw")
i <- grep("\\[final\\]", d)
d[i] <- "\\usepackage[final]{changes}"
d[i-1] <- "% \\usepackage{changes}"
writeLines(d, "sdmTMB-paper.Rnw")

system('R -e \'knitr::knit("sdmTMB-paper.Rnw");tinytex::latexmk("sdmTMB-paper.tex", engine = "xelatex", clean = FALSE)\'')

system("mv sdmTMB-paper.pdf sdmTMB-paper-V3.pdf")

## revert state:
d <- readLines("sdmTMB-paper.Rnw")
i <- grep("response2}", d)
d[seq(i-1, i+1)] <- gsub("\\% ", "", d[seq(i-1, i+1)])
i <- grep("\\[final\\]", d)
d[i] <- "% \\usepackage[final]{changes}"
d[i-1] <- "\\usepackage{changes}"
writeLines(d, "sdmTMB-paper.Rnw")

# system('R -e \'knitr::knit("sdmTMB-paper.Rnw");tinytex::latexmk("sdmTMB-paper.tex", engine = "xelatex", clean = FALSE)\'')

setwd(here::here())
