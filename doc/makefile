TEXT = paper
REFS = refs

all: $(TEXT).pdf appendix-model.pdf

%.pdf: %.Rmd $(REFS).bib
	Rscript -e "rmarkdown::render('paper.Rmd')"
continuous:
	while true; do make --silent; sleep 3; done

appendix-model.pdf: appendix-model.Rmd refs.bib
	Rscript -e "rmarkdown::render('appendix-model.Rmd')"
