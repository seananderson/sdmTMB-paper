all: paper.pdf sdmTMB-cover.pdf

paper.pdf: paper.Rmd refs.bib
	Rscript -e "rmarkdown::render('paper.Rmd')"

sdmTMB-cover.pdf: sdmTMB-cover.md template.tex
	Rscript -e "rmarkdown::render('sdmTMB-cover.md')"
