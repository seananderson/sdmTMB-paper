knitr::purl("paper-jss.Rmd", documentation = 1L)

# fix this:


if (!file.exists(here::here("data/ne_10m_lakes"))) {
  zip_file <- paste0("https://www.naturalearthdata.com/http//www.naturalearthdata.com/",
    "download/10m/physical/ne_10m_lakes.zip")
  download.file(zip_file, destfile = here::here("data/ne_10m_lakes.zip"))
  unzip(here::here("data/ne_10m_lakes.zip"), exdir = here::here("data/ne_10m_lakes"))
}


## ----shapes-read, echo=TRUE-------------------------------------------------------------------------------------------------
coast <- rnaturalearth::ne_coastline(scale = "medium", returnclass = "sf") %>%
  sf::st_transform(crs = Albers)
lakes <- sf::st_read(here::here("data/ne_10m_lakes"), quiet = TRUE)
lakes <- lakes[lakes$scalerank == 0, ] %>% sf::st_transform(crs = Albers)



# trash this:

## ----compare-table, results='asis', eval=FALS
# ....

# put this in root folder:

snow <- readRDS(here::here("data/SNOW_data.rds"))


# look for all here::here



knitr::spin("pcod-fig.R")
knitr::spin("timing.R")




# add timing script onto bottom

# go change:
saveRDS(out, file = "analysis/timing-cache-parallel-2022-11-01.rds")
out <- readRDS("analysis/timing-cache-parallel-2022-11-01.rds")

# at bottom of added timing

# add sessionInfo()  at end!

# remove all child chunks
# remove knitr::opts lines except first?
# remove   optipng = "-strip all"
# remove knitr::knit_hooks$set(optipng = knitr::hook_optipng)

# add pcod-fig.R stuff to bottom
