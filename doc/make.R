files <- c(
  "sdmTMB-paper",
  "appendix-model",
  "appendix-pcod-sdmTMB",
  "appendix-snow-sdmTMB",
  "appendix-binomial-inla",
  "appendix-pcod-VAST-tweedie",
  "appendix-validation"
)

# https://github.com/rstudio/rmarkdown/issues/1673
render_separately <- function(...) callr::r(
  function(...) rmarkdown::render(..., envir = globalenv()), args = list(...), show = TRUE)

for (i in seq_along(files)) {
  render_separately(paste0(here::here("doc", files[i]), ".Rmd"))
}

names_to <- c(
  "sdmTMB-paper.pdf",
  "appendix1-model.pdf",
  "appendix2-index-sdmTMB.pdf",
  "appendix3-owl-sdmTMB.pdf",
  "appendix4-INLA-comparison.pdf",
  "appendix5-VAST-comparison.pdf",
  "appendix6-speed-validation.pdf"
)

loc <- "~/Dropbox/sdmTMB-paper/"
dir.create(loc, showWarnings = FALSE)

for (i in seq_along(files)) {
  from <- paste0(here::here("doc", files[i]), ".pdf")
  to <- paste0(loc, names_to[i])
  cat(to, "\n")
  file.copy(from, to, overwrite = TRUE)
}

system("cd doc; ./clean-bib.sh")


# f <- paste(names_to, collapse = " ")

# system(paste0("cd ~/Documents/sdmTMB-paper/;'/System/Library/Automator/Combine PDF Pages.action/Contents/Resources/join.py' -o sdmTMB-paper-combined.pdf ", f))

# system(paste0("cd ~/Dropbox/sdmTMB-paper/;/opt/homebrew/bin/gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=sdmTMB-paper-combined.pdf ", f))

system("open ~/Dropbox/sdmTMB-paper/")
