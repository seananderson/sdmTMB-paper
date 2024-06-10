old_folder <- "~/Downloads/sdmTMB-paper-53264baf950427e13209950484b881712f1d74c4/doc/paper-jss/"
new_folder <- "~/src/sdmTMB-paper/doc/paper-jss/"

wd <- getwd()
setwd(new_folder)

f2 <- readLines(paste0(new_folder, "sdmTMB-paper.tex"))
for (i in seq_along(f2)) {
  x <- f2[i]
  x <- gsub("\\\\Rev\\{[a-zA-Z0-9]+\\}", "", x) # reviewer reference tags
  f2[i] <- x
}

# remove response letter at end:

resp <- grep("response2", f2)
f2 <- c(f2[seq(1, resp-1)], "\\end{document}")

writeLines(f2, paste0(new_folder, "/main-clean.tex"))
system(paste0("latexdiff ", old_folder, "/sdmTMB-paper.tex main-clean.tex > diff.tex"))

system("pdflatex diff.tex")
system("bibtex diff")
system("pdflatex diff.tex")
system("pdflatex diff.tex")
system("pdflatex diff.tex")
system("pdflatex diff.tex")

system("pdflatex main-clean.tex")
system("bibtex main-clean")
system("pdflatex main-clean.tex")
system("pdflatex main-clean.tex")
system("pdflatex main-clean.tex")
system("pdflatex main-clean.tex")

system("cp main-clean.pdf anderson-etal-main.pdf")
system("cp diff.pdf anderson-etal-diff.pdf")

system("pdflatex main.tex")
system("bibtex main")
system("pdflatex main.tex")
system("pdflatex main.tex")
system("pdflatex main.tex")
system("pdflatex main.tex")

system("open .")

setwd(wd)
