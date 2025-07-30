# Prepare environment ####
rm(list = ls()) # Reset R`s brain

library(vegan)
library(tidyverse)

# Load data, using species-sample matrix (in wide format) 

data_spe <- read_csv("./nuts/data_wide_spe.csv") %>% 
  as.data.frame()

# Assign the first column as row names
rownames(data_spe) <- data_spe[, 1]
# Remove the first column
data_spe <- data_spe[, -1]

data_spe

# Load environmental data
env_data <- read_csv("./nuts/nuts_data_env.csv") %>% 
  unite("sample_id", region, sample, remove = F) %>% 
  arrange(region, sample) %>% 
  column_to_rownames("sample_id")

env_data

## NMDS (Non-metric Multi-Dimensional Scaling) ####
set.seed(123)
nmds.nuts <- metaMDS(data_spe,           # Our community-by-species matrix
                     distance = "bray",  # Ecological distance
                     k = 2,              # The number of reduced dimensions
                     trymax = 100)       # The number of default iterations

# metaMDS has automatically applied a square root transformation and calculated the 
# Bray-Curtis distances for our community-by-site matrix.

# Shepard plot, which shows scatter around the regression between the interpoint distances 
# in the final configuration (i.e., the distances between each pair of communities) against 
# their original dissimilarities.

stressplot(nmds.nuts) # Large scatter around the line suggests that original dissimilarities are not well preserved in the reduced number of dimensions

# Basic plot function
plot(nmds.nuts)
plot(nmds.nuts, type = "t")

# Filter and fit environmental data to the ordination
envfit_vars <- env_data %>%
  select(energy, proteins, fats, carbohydrates, price)

fit <- envfit(nmds.nuts, envfit_vars, permutations = 999)

# Plot environmental data vectors
plot(fit)

# Step-by-step plotting
ordiplot(nmds.nuts, type="n")                                      # coordinate plot
orditorp(nmds.nuts, display="species",col="red", air=0.01)            # species with names
orditorp(nmds.nuts, display="sites",cex=1,air=0.01)               # groups


# STEP 6: Plot with ggplot2
library(ggplot2)
library(ggrepel)

# Extract NMDS scores
site_scores <- as.data.frame(scores(nmds.nuts, display = "sites"))
site_scores$sample_id <- rownames(site_scores)

species_scores <- as.data.frame(scores(nmds.nuts, display = "species"))
species_scores$species <- rownames(species_scores)

vectors <- as.data.frame(scores(fit, display = "vectors"))
vectors$variable <- rownames(vectors)

# Combine with environmental metadata
# Assumes rows are in same order!
nmds_df <- bind_cols(site_scores, env_data)

# Plot ordination with samples as points and environmental metadata as vectors
ggplot(nmds_df, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(colour = region, shape = origin), size = 4) +
  scale_size_continuous(range = c(2, 6)) +
  theme_minimal() +
  labs(title = "NMDS Ordination of Nut & Fruit mix data",
       colour = "Region",
       shape = "Origin") + 
  geom_segment(data = vectors,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  geom_text(data = vectors,
            aes(x = NMDS1, y = NMDS2, label = variable),
            hjust = 1.1, vjust = 1.1)


# We can also add species scores:
ggplot(nmds_df, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(colour = region, shape = origin), size = 5) +
  labs(title = "NMDS Ordination of Nut & Fruit mix data",
       colour = "Region",
       shape = "Origin") +
  geom_text_repel(data = species_scores,
                  aes(x = NMDS1, y = NMDS2, label = species),
                  colour = "darkred", size = 3.5) +
  geom_segment(data = vectors,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.2, "cm")),
               colour = "black") +
  geom_text(data = vectors,
            aes(x = NMDS1, y = NMDS2, label = variable),
            hjust = 1.1, vjust = 1.1, size = 4.5) +
  theme_bw() 

ggsave("./figures/nuts_nmds.png", width = 8.3, height = 5, units = "in", dpi = 300)
