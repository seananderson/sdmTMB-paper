TEXT = paper
REFS = refs

all: $(TEXT).pdf

%.pdf: %.Rmd $(REFS).bib
	Rscript -e "rmarkdown::render('paper.Rmd')"
continuous:
	while true; do make --silent; sleep 3; done
