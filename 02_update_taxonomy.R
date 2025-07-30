# Set the Environment ####
rm(list = ls()) # Reset R's brain

# Load necessary libraries
library(readxl)   # to read Excel files
library(dplyr)    # fro data manipulations
library(tidyr)    # fro data manipulations
library(ggplot2)  # for data visualization
library(sf)       # for spatial data and mapping
library(viridis)  # for color palette

# Load Ukraine admin boundaries from a file
load("./data/Ukraine_poly_adm1.Rdata")


# Red data ####
# OTU table
otu_table <- readxl::read_excel("Phyloseq_Tables.xlsx",
                                sheet = "OTU_table")
# Sample data
samples <- readxl::read_excel("Phyloseq_Tables.xlsx",
                              sheet = "Samples")
# Taxonomy assignment
taxonomy <- readxl::read_excel("Phyloseq_Tables.xlsx",
                              sheet = "Taxonomy") %>% 
  select(-Reference_Sequence)

unite_id <- read.csv("./data/unite_blastresult.csv") %>% 
  rename(id = occurrenceId)


# Re-calculate richness ####

# Merge taxonomy with UNITE blast results
sh_abundance <- otu_table %>%
  # Join OTU table with SH annotations
  left_join(unite_id %>% select(id, scientificName), by = "id") %>% 
  # Group by SH
  group_by(scientificName) %>%
  # and sum abundances across OTUs
  summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE))) %>%
  ungroup()

sh_abundance

# Drop non-identified OTUs
sh_abundance <- slice(sh_abundance, -1)


# Convert to long format and compute richness & read depth
sample_metrics <- sh_abundance %>%
  pivot_longer(-scientificName, names_to = "sample", values_to = "abundance") %>%
  group_by(sample) %>%
  summarise(
    richness_detected = sum(abundance > 0),
    Read_depth = sum(abundance),
    .groups = "drop"
  )

ggplot(sample_metrics, aes(x = Read_depth, y = richness_detected)) +
  geom_point() +
  geom_abline(slope = model$coefficients[2], intercept = model$coefficients[1])
# Fit linear model and compute residuals
model <- lm(richness_detected ~ Read_depth, data = sample_metrics)
sample_metrics$richness_adj <- resid(model)

# Visualize adjusted richness ####

# Join with metadata and convert to sf object
samples_sf <- samples %>%
  left_join(sample_metrics, by = c("id" = "sample")) %>%
  filter(!is.na(richness_adj)) %>%  # remove samples without abundance info
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

# Make a plot
ggplot(samples_sf) +
  geom_sf(data = ukraine_poly) +
  geom_sf(aes(colour = richness_adj), size = 3) + 
  scale_colour_viridis(option = "D", name = "Read Depth -\nAdjusted Richness") +
  theme_minimal() +
  labs(title = "Adjusted EcM Fungal Richness per Sample",
       subtitle = "Residuals of linear model: richness ~ read depth",
       caption = "SH-based richness, UNITE v9.0")

# Make an interactive plot
library(mapview)
mapview(samples_sf, zcol = "richness_adj")

# End of the script ####

