files <- c(
  "paper.pdf",
  "appendix-model.pdf",
  "appendix-pcod-sdmTMB.pdf",
  "appendix-binomial-inla.pdf",
  "appendix-pcod-VAST-tweedie.pdf",
  "appendix-validation.pdf"
)

names <- c(
  "sdmTMB-paper.pdf",
  "appendix1-model.pdf",
  "appendix2-index-sdmTMB.pdf",
  "appendix3-INLA-comparison.pdf",
  "appendix4-VAST-comparison.pdf",
  "appendix5-speed-validation.pdf"
)

loc <- "~/Dropbox/sdmTMB-paper/"

for (i in seq_along(files)) {
  from <- here::here("doc", files[i])
  to <- paste0(loc, names[i])
  cat(to, "\n")
  file.copy(from, to, overwrite = TRUE)
}
