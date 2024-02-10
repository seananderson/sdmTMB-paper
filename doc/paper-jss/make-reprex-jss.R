RUN_SPIN <- TRUE

setwd(here::here())
setwd("doc/paper-jss/")
dir.create("reprex", showWarnings = FALSE)

knitr::purl("sdmTMB-paper.Rnw", documentation = 1L)

system("mv sdmTMB-paper.R reprex/sdmTMB-paper.R")

d <- readLines("reprex/sdmTMB-paper.R")
# d[grepl("## ----child=", d)]
# d <- d[!grepl("## ----child=", d)]

do_replace <- function(find, replacement) {
  cat("found:\n")
  i <- grep(find, d)
  print(d[i])
  d[grepl(find, d)] <- replacement
  cat("after replacement:\n")
  print(d[i])
  d
}

d[grep("knitr::render_sweave()", d)]
d <- do_replace("knitr::render_sweave()", "# knitr::render_sweave()")

d[grep("dpi = 140,", d)]
d <- do_replace("dpi = 140,", "  dpi = 140,")

d[grep(" cache = TRUE,", d)]
d <- do_replace(" cache = TRUE,", " # cache = TRUE,")
d <- do_replace(" autodep = TRUE,", " # autodep = TRUE,")

d[grep("-strip all", d)]
d <- do_replace("-strip all", "  # optipng = \"-strip all\"")

d[grep("optipng = knitr", d)]
d <- do_replace("optipng = knitr", "# knitr::knit_hooks$set(optipng = knitr::hook_optipng)")

# i <- grep("## ----compare-table", d)
# j <- grep("## ----setup-pcod", d)
# d <- d[-seq(i, j-2)]

# i <- grep("owl-knitr-setup", d)
# j <- grep("setup-owls", d)
# d <- d[-seq(i, j-1)]


d[grep("snow <- read", d)]
d <- do_replace("snow <- read", "snow <- readRDS(\"snow-data.rds\")")

# d[grep("owl <- here::", d)]
# d <- do_replace("owl <- here::", "# owl <- here::here(\"figs\", \"owl-nao-effect.png\")")

# d[grep("knitr::include_graphics\\(owl", d)]
# d <- do_replace("knitr::include_graphics\\(owl", "# knitr::include_graphics(owl)")

d[grep("data/ne_10m_lakes", d)]
d <- gsub("data/ne_10m_lakes", "ne_10m_lakes", d)

d[grep("lakes <- sf::st_read", d)]
d <- do_replace("lakes <- sf::st_read", "lakes <- sf::st_read(\"ne_10m_lakes\", quiet = TRUE)")

d[grep("file.exists\\(here::here\\(\"ne_10m_lakes\\\"", d)]
d <- do_replace("file.exists\\(here::here\\(\"ne_10m_lakes\\\"", "if (!file.exists(\"ne_10m_lakes\")) {")

d[grep("here::here\\(\"ne_10m_lakes.zip\"\\)", d)]
d <- gsub("here::here\\(\"ne_10m_lakes.zip\"\\)", "\"ne_10m_lakes.zip\"", d)
d <- gsub("here::here\\(\"ne_10m_lakes\"\\)", "\"ne_10m_lakes\"", d)

d <- c(d, "sessionInfo()")

d[grepl("here::", d)]

# i <- grep("## ----setup-owls", d)
# d <- c(d[seq(1, i-1)], "\n# Snowy Owl appendix -----------------------------------------------------------------------\n", d[seq(i, length(d))])

# i <- grep("## ----setup-pcod", d)
# d <- c(d[seq(1, i-1)], "\n# Pacific Cod appendix -----------------------------------------------------------------------\n", d[seq(i, length(d))])

# i <- grep("## ----inla-knitr-setup", d)
# d <- c(d[seq(1, i-1)], "\n# INLA comparison appendix -----------------------------------------------------------------------\n", d[seq(i, length(d))])

i <- grep("dog-update-old, ", d) + 1
j <- grep("## fit_dpl <- update", d)
for (s in i:j) {
  d[s] <- gsub("^## ", "", d[s])
}

i <- grep("dog-update,", d) - 1
j <- grep("dog-update-old,", d) - 2
d <- c(d[1:i], "if (new_deltas) {", d[(i+1):(j-1)], "}", d[(j):length(d)])

i <- grep("dog-update-old,", d) - 1
j <- grep("----dog-aic,", d) - 2
d <- c(d[1:i], "if (!new_deltas) {", d[(i+1):(j-1)], "}", d[(j):length(d)])

i <- grep("inlabrufit-timing", d) + 1
j <- grep("bru_max_iter = 1, num.threads = \"1:1\"", d) + 2
for (s in i:j) {
  d[s] <- gsub("^## ", "", d[s])
}

# these make figs not show up in html:
d <- gsub('dev = \\"pdf\\"', 'dev = \\"png\\"', d)
d <- gsub('out.width=\\"[0-9]+in\\", ', "", d)

# check:
d[grepl("out.width=\\'[0-9]+in\\'", d)]

writeLines(d, "reprex/replication-code-main.R")

system("cp ~/src/sdmTMB-paper/data/SNOW_data.rds ~/src/sdmTMB-paper/doc/paper-jss/reprex/snow-data.rds")

system("rm reprex/sdmTMB-paper.R")

setwd("reprex")
if (RUN_SPIN) system("R -e 'knitr::spin(\"replication-code-main.R\")'")

d <- readLines(here::here("analysis/timing.R"))
writeLines(d, "replication-code-timing.R")
if (RUN_SPIN) system("R -e 'knitr::spin(\"replication-code-timing.R\")'")
#
# d <- readLines(here::here("analysis/pcod-fig.R"))
# # d <- gsub('silent = TRUE', 'silent = FALSE', d)
# writeLines(d, "replication-code-pcod.R")
# if (RUN_SPIN) system("R -e 'knitr::spin(\"replication-code-pcod.R\")'")

setwd(here::here())
