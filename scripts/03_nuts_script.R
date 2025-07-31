# Set the Environment ####
rm(list = ls()) # Reset R's brain

# Load necessary libraries
library("tidyverse")

# Load from csv file
data <- read_csv("./nuts/nuts_data_raw.csv")

# Check data
head(data)
summary(data)

# Check for misspellings
unique(data$species)

data$species[
  data$species == "peanuts"
] <- "peanut"


# Transform raw data to abundance table and wide format
data %>%
  group_by(region, sample) %>%
  summarise(n_ids = n(), .groups = "drop")

data_abund <- data %>% 
  group_by(region, sample, species) %>% 
  summarise(abundance = n(), .groups = "drop")

head(data_abund)

# Convert to wide format
data_wide <- data_abund %>% 
  pivot_wider(names_from = species,
              values_from = abundance,
              values_fill = 0) %>% 
  unite("sample_id", region, sample, remove = T) %>% 
  column_to_rownames("sample_id")

write.csv(data_wide, "./nuts/data_wide_spe.csv")


# Let's explore data!
# Look at species richness

data_regions <- data_abund %>% 
  group_by(region, species) %>% 
  summarise(abundance = sum(abundance))

library(vegan) # vegan = "vegetation analysis"

# How many species per region?
SR_region <- data_regions %>% 
  group_by(region) %>% 
  summarise(SR = specnumber(abundance))

SR_region

# How many species per sample?
SR_sample <- data_abund %>% 
  group_by(region, sample) %>% 
  summarise(SR = specnumber(abundance))
  
SR_sample

# Mean species richness per region
mean_SR_sample <- SR_sample %>%
  summarise(mSR = mean(SR))
  
mean_SR_sample

# Add mean species richness as a variable
SR_all <- SR_region %>% 
  mutate(mSR_sample = mean_SR_sample$mSR)

# Create simple plot
library(ggplot2)

ggplot(SR_sample, aes(x = region, y = SR, colour = region)) +
  geom_boxplot() +
  geom_point(size = 3) +
  labs(
    title = "Species richness",
    y = "species richness"
  ) +
  theme_bw()

?specnumber


# Abundance curves
dat_landscapes_sort <- data_regions |>  
  arrange(region, desc(abundance))

dat_landscapes_sort <- dat_landscapes_sort |> 
  group_by(region) |> 
  mutate(rank = row_number())

ggplot(dat_landscapes_sort, 
       aes(x = rank, y = abundance, color = region, shape = region)) +
  geom_line() +
  geom_point(size = 3) +
  scale_y_log10() +
  ylab("abundance") +
  ggtitle("Rank-abundance-curves all regions") +
  theme_bw()


# Rarefaction-based richness estimation ####
# More: https://sites.google.com/view/chao-lab-website/software/inext?authuser=0

# install.packages("iNEXT")
library(iNEXT)

# Abundance data
# Transpose the data so species are rows, samples are columns
data_t <- t(data_wide)

# Convert to a list of abundance vectors per sample
abundance_list <- apply(data_t, 2, as.integer)

out <- iNEXT(abundance_list, q = 0, datatype = "abundance")

ggiNEXT(out, type = 1)  # sample-size-based rarefaction/extrapolation


# calculate estimated richness based on rarefaction
chao1 <- iNEXT::ChaoRichness(t(data_wide), datatype = "abundance")
richness_est <- cbind(rownames(chao1), chao1$Estimator) %>% 
  as.data.frame() %>% 
  rename(plot = 1,
         richness = 2) %>% 
  mutate(richness = as.numeric(richness))

richness_est


# Incidence approach - sample coverage
data_t_pa <- ifelse(data_t>0, 1, 0)

m2 <- list(data_t_pa)
names(m2) <- "IncidenceData - Nuts"
out.pf <- iNEXT(m2, q = c(0,1,2), datatype = "incidence_raw")

ggiNEXT(out.pf, type = 2)

chao1 <- iNEXT::ChaoRichness(m2, datatype = "incidence_raw")

chao1



# Incidence approach - species richness
# Convert counts to presence/absence
data_inc <- data_wide > 0

# Count in how many samples each species occurs
species_incidence_freq <- colSums(data_inc)

# Total number of sampling units
n_units <- nrow(data_wide)

# Create the input for iNEXT
incidence_input <- list(c(n_units, species_incidence_freq))
names(incidence_input) <- "IncidenceData - Nuts"

out_inc <- iNEXT(incidence_input, q = 0, datatype = "incidence_freq")

ggiNEXT(out_inc, type = 1)


# End of the script ####

