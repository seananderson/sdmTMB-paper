
library(neonstore)
library(tidyverse)
library(ISOweek)
library(neonUtilities)
# access data using EW's token
#neon_download("DP1.10022.001", type = "expanded",.token=NEON_TOKEN)

beetles <- neon_download("DP1.10022.001", type = "expanded",.token=NEON_TOKEN)
neon_store(product = "DP1.10022.001")

# Update `scientificName`, `taxonID`, `taxonRank` and `morphospeciesID` using assignments from parataxonomy and expert taxonomy.
#
library(dplyr)
library(stringi)

clean_names <- function (x)
{
  s <- stringi::stri_split_regex(x, "/", simplify = TRUE)[,1]
  s <- stringi::stri_extract_all_words(s, simplify = TRUE)
  if (dim(s)[2] > 1)
    stringi::stri_trim(paste(s[, 1], s[, 2]))
  else stringi::stri_trim(s[, 1])
}

resolve_taxonomy <- function(sorting, para, expert){

  taxonomy <-
    left_join(sorting,
              select(para, subsampleID, individualID, scientificName, taxonRank, taxonID, morphospeciesID),
              by = "subsampleID")  %>%
    ## why are there so many other shared columns (siteID, collectDate, etc?  and why don't they match!?)
    ## we use `select` to avoid these
    left_join(
      select(expert, -uid, -namedLocation, -domainID, -siteID, -collectDate, -plotID, -setDate, -collectDate),
      by = "individualID") %>%
    distinct() %>%
    ## Prefer the para table cols over the sorting table cols only for sampleType=="other carabid"
    mutate(taxonRank.x = ifelse(is.na(taxonRank.y) | sampleType != "other carabid", taxonRank.x, taxonRank.y),
           scientificName.x = ifelse(is.na(scientificName.y) | sampleType != "other carabid", scientificName.x, scientificName.y),
           taxonID.x = ifelse(is.na(taxonID.y) | sampleType != "other carabid", taxonID.x, taxonID.y),
           morphospeciesID.x =  ifelse(is.na(morphospeciesID.y) | sampleType != "other carabid", morphospeciesID.x, morphospeciesID.y)) %>%
    ## Prefer expert values where available
    mutate(taxonRank = ifelse(is.na(taxonRank), taxonRank.x, taxonRank),
           scientificName = ifelse(is.na(scientificName), scientificName.x, scientificName),
           taxonID = ifelse(is.na(taxonID), taxonID.x, taxonID),
           morphospeciesID =  ifelse(is.na(morphospeciesID), morphospeciesID.x, morphospeciesID),
           nativeStatusCode = ifelse(is.na(nativeStatusCode.y), nativeStatusCode.x, nativeStatusCode.y),
           sampleCondition = ifelse(is.na(sampleCondition.y), sampleCondition.x, sampleCondition.y)
    ) %>%
    select(-ends_with(".x"), -ends_with(".y")) # %>%
  #    select(-individualCount)
  ## WARNING: if the subsample is split into separate taxa by experts, we do not know
  ## how many of the total count should go to each taxon in the the split
  ## since only part of that subsample have been pinned.
  ## We should flag these cases in some manner.


  #### Should we add a "species" column, using morphospecies or the best available?
  ## Use morphospecies if available for higher-rank-only classifications,
  ## Otherwise, binomialize the scientific name:
  taxonomy <- taxonomy %>%
    mutate(morphospecies =
             ifelse(taxonRank %in% c("subgenus", "genus", "family", "order") & !is.na(morphospeciesID),
                    morphospeciesID,
                    clean_names(scientificName)
             )
    )

  ## Beetles must be identified as carabids by both sorting table and the taxonomists (~3 non-Carabidae slip through in sorting)
  beetles <- taxonomy %>%
    filter(grepl("carabid", sampleType)) %>%
    filter(family == "Carabidae" | is.na(family))

  beetles
}

## Load data from raw files
sorting <- neon_table("bet_sorting-expanded")
para <- neon_table("bet_parataxonomistID-expanded")
expert <- neon_table("bet_expertTaxonomistIDProcessed-expanded")
field <- neon_table("bet_fielddata-expanded")


#### Generate derived richness table  ####################
beetles <- resolve_taxonomy(sorting, para, expert) %>%
  mutate(iso_week = ISOweek::ISOweek(collectDate),
         time = ISOweek::ISOweek2date(paste0(iso_week, "-1"))) %>%
  as_tibble()

richness <- beetles %>%
  select(taxonID, siteID, collectDate, time) %>%
  distinct() %>%
  count(siteID, time) %>%
  rename(richness = n)  %>%
  ungroup()



#### Generate derived abundance table ####################

## Does not reflect taxonomic corrections!
## Allows for some counts even when richness is NA

effort <- field %>%
  mutate(iso_week = ISOweek::ISOweek(collectDate),
         time = ISOweek::ISOweek2date(paste0(iso_week, "-1"))) %>%
  group_by(siteID, time) %>%
  summarise(trapnights = as.integer(sum(collectDate - setDate)),
            .groups = "drop")

counts <- beetles %>%
  mutate(iso_week = ISOweek::ISOweek(collectDate),
         time = ISOweek::ISOweek2date(paste0(iso_week, "-1"))) %>%
  group_by(siteID, time) %>%
  summarise(count = sum(as.numeric(individualCount), na.rm = TRUE),
            .groups = "drop")

abund <- counts %>%
  left_join(effort) %>%
  arrange(time) %>%
  mutate(abundance = count / trapnights) %>%
  select(siteID, time, abundance) %>%
  ungroup()

targets_na <- full_join(abund, richness)
# join back in site info
fields = dplyr::group_by(field, siteID) %>%
  dplyr::summarize(decimalLatitude = decimalLatitude[1],
                   decimalLongitude = decimalLongitude[1])

targets_na = dplyr::left_join(targets_na, fields)

targets_na = dplyr::select(targets_na, decimalLatitude,
                           decimalLongitude,
                           abundance,
                           richness,
                           time)
saveRDS(targets_na, "beetles.rds")
