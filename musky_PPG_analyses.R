#########################################################
#                       musky PPGs                      #
#########################################################
# Put files in directory /PPGs/for_manuscript.


## libraries required ####
options(contrasts = c("contr.sum", "contr.poly"))

# libraries required
library(tidyverse)
library(forcats)
library(broom)
library(vegan)
library(car)
library(geosphere)
library(RColorBrewer)
library(patchwork)
library(emmeans)
library(ggrepel)
library(scales)
library(ggpattern)
library(gridExtra)
library(cowplot)

################### plotting aesthetics ########################################


# aesthetics w/ legend
ng1 <- theme(aspect.ratio=0.7,panel.background = element_blank(),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             panel.border=element_blank(),
             axis.line.x = element_line(color="black",size=1),
             axis.line.y = element_line(color="black",size=1),
             axis.ticks=element_line(color="black"),
             axis.text=element_text(color="black",size=15),
             axis.title=element_text(color="black",size=1),
             axis.title.y=element_text(vjust=2,size=17),
             axis.title.x=element_text(vjust=0.1,size=17),
             axis.text.x=element_text(size=15),
             axis.text.y=element_text(size=15),
             strip.text.x = element_text(size = 10, colour = "black",face = "bold"),
             strip.background = element_rect(colour="black"),
             legend.position = "right", legend.direction="vertical",
             legend.text=element_text(size=17), legend.key = element_rect(fill = "white"),
             legend.title = element_text(size=17, face = "bold", hjust = 0.5),
             legend.key.size = unit(1.0, "cm"),
             legend.background = element_blank(),
             legend.box.background = element_rect(colour = "black"))

# aesthetics w/out legend
ng2 <- theme(aspect.ratio=0.7,panel.background = element_blank(),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             panel.border = element_blank(),
             axis.line.x = element_line(color="black",size=1),
             axis.line.y = element_line(color="black",size=1),
             axis.ticks=element_line(color="black"),
             axis.text=element_text(color="black",size=15),
             axis.title=element_text(color="black",size=1),
             axis.title.y=element_text(vjust=2,size=17),
             axis.title.x=element_text(vjust=0.1,size=17),
             axis.text.x=element_text(size=15),
             axis.text.y=element_text(size=15),
             strip.text.x = element_text(size = 17, colour = "black",face = "bold"),
             strip.background = element_rect(colour="black"),
             legend.position = "none")

########

## NMDS for chemical defense analyses ####
## import data

# concentrations of PPGs
chems <- read_csv("PPGs/for_manuscript/Musky_PPG_conc.csv") # inds are from different planting than cg1 & cg2
# each line and rep combination is a unique maternal line

chem_vals <- select(chems, `PPG b`:`Unkn 16`) %>% 
  as.matrix()

# env data
musky_monk_master <- read_csv("PPGs/for_manuscript/musky_monk_env.csv") %>% 
  mutate(pop_num = as.character(pop_num))

# dissimilarity matrix
chem_dis <- vegdist(chem_vals, method = "bray")

# checking for optimal # of dims
NMDS.scree <- function(x) { #where x is the name of the data frame variable
  plot(rep(1, 10), replicate(10, metaMDS(x, autotransform = F, k = 1)$stress), 
       xlim = c(1, 10),ylim = c(0, 0.35), 
       xlab = "# of Dimensions", ylab = "Stress", main = "NMDS stress plot")
  for (i in 1:10) {
    points(rep(i + 1,10),replicate(10, metaMDS(x, autotransform = F, k = i + 1)$stress))
  }
}

NMDS.scree(chem_dis) # two dimensions appear to be sufficient

# adding env vars
chem_env <- chems %>% 
  separate(`Sample Name`, c("pop_num", "ind", "rep", "leaf_pair")) %>% 
  unite("ID", c("pop_num", "ind", "rep"), remove = FALSE) %>% 
  left_join(., select(musky_monk_master, pop_name:range), by = "pop_num")



# make NMDS
set.seed(321)

NMDS_PPG <- metaMDS(chem_dis, k = 2)

ef <- envfit(NMDS_PPG, select(chem_env, `PPG b`:`Unkn 16`), perm = 999)

# plot NMDS
plot(NMDS_PPG, display = "sites")
plot(ef, p.max = 0.05, col = "black")


##

NMDS_scores <- vegan::scores(ef, display = "vectors")


### add NMDS scores ####
axes_scores <- vegan::scores(NMDS_PPG, choices = c(1, 2), tidy = TRUE)

chem_df <- cbind(chem_env, axes_scores)



##### R color brewer palette colours ####
display.brewer.pal(7, "Set2")

brewer.pal(n = 7, "Set2")

# HEX codes for colours for each region

# Willamette #66C2A5

# Sierra #FFD92F

# Sequoia #E5C494



########

## How do total PPGs change through development and by population? ####

#### what predicts total PPGs ####
PPG_pop_develop <- lm(`Total PPGs` ~ leaf_pair*pop_name, data = chem_df)

# check residuals
plot(PPG_pop_develop) # looks good enough 

plot(density(resid(PPG_pop_develop))) # pretty normal actually

summary(PPG_pop_develop)

ANOVA_pop_develop <- Anova(PPG_pop_develop, type = "3") 

ANOVA_pop_develop

# post hoc analysis by pop
PPG_pop_develop_pop_name <- emmeans(PPG_pop_develop, tukey ~ pop_name) 
PPG_pop_develop_pop_name

# for more precise values
as.data.frame(PPG_pop_develop_pop_name)

# Willamette (TAB) - Sequoia (CJC) (p = 0.628, t = 0.925) ∆ = 3.43 ± 3.71
# Willamette (TAB) - Sierra (SRN) (p < 0.0001, t = -5.870) ∆ = -21.78 ± 3.71
# Sequoia (CJC) - Sierra (SRN) (p < 0.0001, t = -6.795) ∆ = -25.22 ± 3.71

# post hoc analysis by leaf pair
PPG_pop_develop_leaf_pair <- emmeans(PPG_pop_develop, tukey ~ leaf_pair)
PPG_pop_develop_leaf_pair

# for more precise values
as.data.frame(PPG_pop_develop_leaf_pair)

# early set - late set (p = 0.002; t = 3.293) ∆ = 9.98 ± 3.03 mg/g

# or 19.55 %
(51.04563-41.06867)/51.04563*100

#### Total PPGs (Figure 2) ####
chem_df %>%
  mutate(pop_name = factor(pop_name, levels = c("CJC", "SRN", "TAB")),
         pop_name_legend = factor(pop_name, levels = c("TAB", "SRN", "CJC")),
         leaf_pair = factor(dplyr::recode(leaf_pair, "3" = "Early", "8" = "Late"), levels = c("Early", "Late"))) %>%
  ggplot(aes(x = pop_name, y = `Total PPGs`, fill = pop_name_legend, linetype = leaf_pair)) + 
  geom_boxplot(width = 0.5) +
  geom_point(position = position_dodge(width = 0.5), alpha = 0.5, size = 2.5, aes(shape = leaf_pair)) +
  labs(x = "Population",
       y = "Total PPGs (mg/g)",
       fill = "Population",
       linetype = "Ontogeny",
       shape = "Ontogeny") +
  scale_fill_manual(values = c("TAB" = "#66C2A5", "SRN" = "#FFD92F", "CJC" = "#E5C494")) +
  ng1 +
    theme(legend.text = element_text(hjust = 0),
          legend.title = element_text(hjust = 0))


### how do NMDS axes change through development and population ####

#NMDS1
NMDS1_lm <- lm(NMDS1 ~ leaf_pair*pop_name, data = chem_df)
summary(NMDS1_lm)

Anova(NMDS1_lm, type = 3)

#NMDS2
NMDS2_lm <- lm(NMDS2 ~ leaf_pair*pop_name, data = chem_df)
summary(NMDS2_lm)

Anova(NMDS2_lm, type = 3)

## PPG arsenals through development and by population ####

#### what predicts dissimilarity ####

# both terms
set.seed(321)
PPG_arsenal_PERM <- adonis2(chem_dis ~ leaf_pair*pop_name, data = chem_env,
                            permutations = 999, method = "bray", by = "terms")

PPG_arsenal_PERM
# p = 0.001 leaf_pair
# p = 0.001 range
# p = 0.952 leaf_pair*range




#### PPG arsenals (Figure 3) ####

chem_df %>%
  mutate(pop_name = fct_reorder(pop_name, -latitude),
         leaf_pair = dplyr::recode(leaf_pair, "3" = "Early", "8" = "Late")) %>%
  ggplot(aes(x = NMDS1, y= NMDS2)) + 
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.75) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.75) + 
  geom_point(aes(colour = pop_name, shape = leaf_pair), size = 2.5) + 
  geom_line(aes(group = ID, colour = pop_name)) +
  labs(x = "NMDS1", y = "NMDS2", colour = "Population", shape = "Ontogeny") +
  stat_ellipse(aes(colour = pop_name), size = 0.75) +
  scale_colour_manual(values = c("#66C2A5", "#FFD92F", "#E5C494")) +
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  ng1 +
  theme(legend.text = element_text(hjust = 0),
        legend.title = element_text(hjust = 0))


## congeneric comparison ####

### find closest mos and gutt pops for each region ####

# Load df of coordinates of current study M. moschatus pops and Kooyers et al 2017 M. guttatus pops
mos_gutt_pops <- read_csv("PPGs/for_manuscript/mos_gutt_PPG_pops.csv")


# Step 1: Separate the *mos* and *gutt* populations
mos_pops <- mos_gutt_pops %>% filter(species == "mos")
gutt_pops <- mos_gutt_pops %>% filter(species == "gutt")

# Step 2: Create a matrix of coordinates for each group
mos_coords <- mos_pops %>% select(longitude, latitude) %>% as.matrix()
gutt_coords <- gutt_pops %>% select(longitude, latitude) %>% as.matrix()

# Step 3: Calculate the distances between each *mos* population and all *gutt* populations
# Each row represents a *gutt* population, and each column represents a *mos* population
dist_matrix <- distm(gutt_coords, mos_coords, fun = distHaversine)

# Step 4: Find the top 3 closest *gutt* populations for each *mos* population
closest_three_gutt <- data.frame()  # Initialize an empty data frame

# Loop over each column (each `mos` population)
for (i in 1:ncol(dist_matrix)) {
  # Get distances for this `mos` population
  distances <- dist_matrix[, i]
  
  # Get the indices of the top 3 closest `gutt` populations
  closest_indices <- order(distances)[1:3]
  
  # Get elevation for mos and corresponding gutt populations 
  mos_elev <- mos_pops$elevation_m[i]
  gutt_elevs <- gutt_pops$elevation_m[closest_indices]
  
  # Create a data frame with the results for this `mos` population
  temp_df <- data.frame(
    mos_population = rep(mos_pops$pop_name[i], 3),  # Current `mos` population
    closest_gutt_population = gutt_pops$pop_name[closest_indices],  # Closest `gutt` populations
    Distance = distances[closest_indices] / 1000,  # Corresponding distances in km
    Elevation_diff_m = abs(gutt_elevs - mos_elev) 
  )
  
  # Add to the final data frame
  closest_three_gutt <- rbind(closest_three_gutt, temp_df)
}

# Step 5: View results
as_tibble(closest_three_gutt)

# went with HAC for TAB (17.9 km), YVO for SRN (89.7 km) and NAD for CJC (53.1 km)
# elevations - HAC: 1285, TAB : 1440; YVO: 1495, SRN: 1232; NAD: 1689, CJC: 1807


#### stats tests for congeneric comparison ####

# import data and manipulate
# Kooyers et al. 2017 (ran with full data frame NOT included in data repository)
gutt_PPG_df <- read_csv("PPGs/ppg_r_noDKR6.csv") %>% 
  filter(induction == "c") %>%  # remove induced defenses
  filter(Population %in% c("NAD", "YVO", "HAC")) %>% 
  mutate(species = "gutt") %>% 
  rename("Calc A" = calceolariosideA, 
         Conand = conandroside,
         Verb = verbascoside,
         "Calc B" = CalceolariosideB,
         Mimulo = Mimuloside,
         "Unkn 16" = Unkn16,
         "Total PPGs" = Totalppgs,
         "Unkn 10" = Unkn10,
         latitude = Latitude,
         longitude = Longitude,
         elevation_m = Altitude,
         pop_name = Population) %>% 
  group_by(pop_name, Line) %>% 
  mutate(across(`Unkn 10`:`Total PPGs`, ~mean(.))) %>%
  ungroup() %>% 
  distinct(pop_name, Line, .keep_all = TRUE) %>% 
  rename(ind = Line,
         rep = Ind) %>%
  select(-Flat:-Pos, -LineNum, -Flowerdate:-dryweight, -Cline:-`C stable isotope ratio`, -july_pet:-transect) %>% 
  mutate(range = case_when(pop_name == "NAD" ~ "Sequoia",
                           pop_name == "YVO" ~ "Sierra",
                           pop_name == "HAC" ~ "Willamette"))

# writing out simplified data frame for data repository ### i kept uknk 10 here so i can use the same df below for the figure and to make the means file without need to have an ambiguous df of means
write_csv(gutt_PPG_df, file = "PPGs/for_manuscript/Kooyers_et_al_2017_simplified.csv")

# import simplified data frame from Kooyers et al. 2017 included in data repository
gutt_con_df <- read_csv("PPGs/for_manuscript/Kooyers_et_al_2017_simplified.csv") %>% 
  select(-`Unkn 10`)

# current study
mos_con_df <- chem_df %>% 
  filter(leaf_pair == 3) %>%  # remove older leaves
  mutate(species = "mos") %>% 
  select(-ID:-pop_num, -leaf_pair:-`PPG h`, -NMDS1:-label) %>% 
  select(pop_name, ind:species)

# bring em together
mos_gutt_con_df <- rbind(gutt_con_df, mos_con_df) %>% 
  mutate(rel_CalcA = `Calc A` / `Total PPGs`) %>% 
  mutate(rel_Conand = Conand / `Total PPGs`)


## start some stats

# for total PPGs
totalPPGs_mosgutt <- lm(`Total PPGs` ~ range*species, data = mos_gutt_con_df)

Anova(totalPPGs_mosgutt, type = 3)

# for relative calceolarioside A
rel_calcA_mosgutt <- lm(rel_CalcA ~ range*species, data = mos_gutt_con_df)

Anova(rel_calcA_mosgutt, type = 3) 

# for relative conandroside
reL_conand_mosgutt <- lm(rel_Conand ~ range*species, data = mos_gutt_con_df)

Anova(reL_conand_mosgutt, type = 3)


### to make a plot ####
##### make df with means of both species ####

# Calculate means and standard deviations for gutt
gutt_PPG_means <- gutt_PPG_df %>% 
group_by(pop_name) %>%
  summarise(across(`Unkn 10`:`Total PPGs`,
                   list(mean = ~ mean(., na.rm = TRUE),
                        sd = ~ sd(., na.rm = TRUE)),
                   .names = "{col}_{fn}"),
            N = n()) %>%
  ungroup() %>%
  mutate(species = "gutt")

# remove later leaves and remove unecessary columns
filtered_mos_PPG_raw <- chem_df %>%
  filter(leaf_pair == 3) %>%
  select(-ID, -pop_num, -ind, -rep, -leaf_pair)

# Calculate means and standard deviations for mos
mos_PPG_means <- filtered_mos_PPG_raw %>%
  group_by(pop_name) %>%
  summarise(across(`PPG b`:`Total PPGs`,
                   list(mean = ~ mean(., na.rm = TRUE),
                        sd = ~ sd(., na.rm = TRUE)),
                   .names = "{col}_{fn}"),
            N = n()) %>%
  ungroup() %>%
  mutate(species = "mos")

# piece species dfs together
combined_PPG_df <- bind_rows(mos_PPG_means, gutt_PPG_means) %>%
  mutate(range = case_when(
    pop_name %in% c("CJC", "NAD") ~ "Sequoia",
    pop_name %in% c("YVO", "SRN") ~ "Sierra",
    pop_name %in% c("HAC", "TAB") ~ "Willamette",
    TRUE ~ NA_character_)) %>%
  mutate(across(where(is.numeric), ~ replace_na(., 0)))


#### more manipulation and aesthetics ####
combined_PPG_df_se <- combined_PPG_df %>%
  mutate(across(ends_with("_sd"), 
                ~ . / sqrt(N), 
                .names = "{sub('_sd$', '_se', col)}"))

# List of mean columns from specified range
mean_columns <- c("PPG b_mean", "PPG c_mean", "PPG d_mean", "PPG e_mean",
                  "PPG f_mean", "PPG g_mean", "PPG h_mean", "Calc A_mean", 
                  "Conand_mean", "Verb_mean", "Unkn 10_mean", "Calc B_mean", 
                  "Mimulo_mean", "Unkn 16_mean", "Total PPGs_mean")

# Reshape the data for plotting, including mean, se, and range
combined_PPG_df_long <- combined_PPG_df_se %>%
  select(pop_name, range, all_of(mean_columns), ends_with("_se")) %>%  # Select mean and SE columns along with range
  pivot_longer(cols = -c(pop_name, range),                             # Pivot all but "pop_name" and "range"
               names_to = c("measure", "stat"),                        # Split names into "measure" and "stat"
               names_pattern = "(.*)_(mean|se)",                       # Capture suffixes "_mean" or "_se"
               values_to = "value") %>%                                # Store these values in "value"
  pivot_wider(names_from = "stat", values_from = "value") %>%          # Separate into "mean" and "se" columns
  group_by(pop_name, range) %>%                                        # Keep "range" intact for each pop_name
  slice_max(mean, n = 3) %>%                                           # Get top 3 mean values per population 
  ungroup()

# Define custom colors for the populations
custom_colors <- c("TAB" = "#66C2A5",
                   "HAC" = "#A6D8D0",  # Lighter shade of #66C2A5
                   "SRN" = "#FFD92F",
                   "YVO" = "#FFE476",  # Lighter shade of #FFD92F
                   "CJC" = "#E5C494",
                   "NAD" = "#F0DDC5")  # Lighter shade of #E5C494


#### Create plot ####

###### Congeneric comparison (Figure 4) ####
ggplot(combined_PPG_df_long %>% 
         mutate(range = factor(range, levels = c("Willamette", "Sierra", "Sequoia"))) %>%
         mutate(pop_name = factor(pop_name, levels = c("TAB", "HAC", "SRN", "YVO", "CJC", "NAD"))) %>%
         arrange(range, desc(measure)), 
       aes(x = mean, y = pop_name, fill = pop_name, pattern = measure)) +
  # Bar plot with patterns and error bars
  geom_bar_pattern(stat = "identity", 
                   position = position_dodge(width = 0.8), 
                   width = 0.7, 
                   pattern_density = 0.01, 
                   pattern_spacing = 0.07, 
                   pattern_color = "black", 
                   color = "black", 
                   aes(pattern_fill = pop_name)) + 
  
  # Error bars
  geom_errorbar(aes(xmin = mean - se, xmax = mean + se), 
                position = position_dodge(width = 0.8), 
                width = 0.2) +
  
  # Custom color and pattern scales
  scale_fill_manual(values = custom_colors) +
  scale_pattern_manual(values = c("crosshatch", "stripe", "none")) +
  
  # Adjust x-axis to control spacing
  scale_x_continuous(expand = c(0, 1)) +
  
  # Adjust y-axis and theme
  theme_minimal() +
  theme(
    legend.position = "none",                 
    legend.key.size = unit(0.8, "cm"),       
    legend.text = element_text(size = 12),   
    legend.title = element_text(size = 13),  
    strip.text.y = element_text(face = "bold", size = 15),  
    axis.text.y.left = element_text(angle = 90, hjust = 0.5),  
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black", size = 1),
    axis.line.y = element_line(color = "black", size = 1),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black", size = 15),
    axis.title = element_text(color = "black", size = 1),
    axis.title.y = element_text(vjust = 2, size = 17),
    axis.title.x = element_text(vjust = 0.1, size = 17),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    strip.text.x = element_text(size = 10, colour = "black", face = "bold")
  ) +
  
  # Labels and facets
  labs(x = "PPG concentration (mg/g)", y = "Population", fill = "Population", pattern = "Measure") +
  facet_grid(rows = vars(range), scales = "free_y", space = "free")


###### Create legend ####

# Create a data frame with the measure labels and a dummy variable for positioning
legend_data <- data.frame(
  measure = c("Total PPGs", "Conand", "Calc A"),
  x = c(1, 1, 1),  # Single x value to position boxes in a column
  y = c(3, 2, 1)   # y values to control vertical position
)

# Custom pattern styles for each measure
pattern_values <- c("crosshatch", "stripe", "none")

# Create a custom legend plot
legend_plot <- ggplot(legend_data, aes(x = x, y = y, pattern = measure)) +
  # Create square tiles with patterns
  geom_tile_pattern(
    width = 0.4, height = 0.4,               # Control tile size
    fill = "white",                          # White fill inside tiles
    color = "black",                         # Black border for tiles
    pattern_color = "black",                 # Pattern lines in black
    pattern_fill = "black",                  # Pattern fill color
    pattern_spacing = 0.07,                  # Customize pattern spacing
    pattern_density = 0.01                   # Customize pattern density
  ) +
  # Add text labels beside the patterned tiles
  geom_text(aes(label = measure), hjust = -0.5, size = 5) +  # Position text labels
  # Customize the pattern scale to match the patterns
  scale_pattern_manual(values = pattern_values) +
  # Adjust plot limits to provide spacing for labels
  xlim(0.5, 1.8) +
  ylim(0.5, 3.5) +
  # Simplify the theme for the legend-only plot
  theme_void() +  # Use a minimal theme
  theme(
    legend.position = "none",  # Disable default legends
    plot.margin = margin(1, 1, 1, 1, "cm")  # Add padding around the plot
  )

# Display the custom legend plot
print(legend_plot)



########

## misc. ####

#### how do total PPGs change along NMDS axes? ####

#NMDS1
NMDS1_totalPPGlm <- lm(`Total PPGs` ~ NMDS1, data = chem_df)

summary(NMDS1_totalPPGlm)

# Create the scatter plot with regression line
ggplot(chem_df, aes(x = NMDS1, y = `Total PPGs`, color = pop_name)) +
  geom_point(aes(shape = leaf_pair), size = 3, alpha = 0.8) +  # Scatter plot with population-based colors
  # geom_line(aes(group = ID, colour = pop_name)) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +  # Regression line
  ng2 +  # Clean theme
  theme(legend.key = element_blank()) +
  labs(x = "NMDS1", y = "Total PPGs (mg/g)", color = "Population") +  # Axis labels
  scale_color_manual(values = c("TAB" = "#66C2A5", "SRN" = "#FFD92F", "CJC" = "#E5C494")) 

#NMDS2
NMDS2_totalPPGlm <- lm(`Total PPGs` ~ NMDS2, data = chem_df)

summary(NMDS2_totalPPGlm)

# Create the scatter plot with regression line
ggplot(chem_df, aes(x = NMDS2, y = `Total PPGs`, color = pop_name)) +
  geom_point(aes(shape = leaf_pair),size = 3, alpha = 0.8) +  # Scatter plot with population-based colors
  geom_smooth(method = "lm", color = "black", se = FALSE) +  # Regression line
  ng2 +  # Clean theme
  theme(legend.key = element_blank()) +
  labs(x = "NMDS2", y = "Total PPGs (mg/g)", color = "Population") +  # Axis labels
  scale_color_manual(values = c("TAB" = "#66C2A5", "SRN" = "#FFD92F", "CJC" = "#E5C494")) 


#### developmental changes for each PPG by population ####

# note: when normalizing, log(x + 1) OR (1 - x)^1/2

# what about changes for specific PPGs

# Calceolarioside A
calcA_lm_notrans <- lm(`Calc A` ~ leaf_pair*pop_num, data = chem_df)

# residuals
plot(density(resid(calcA_lm_notrans))) # pretty right-skewed

# with log (1 + x) trans
calcA_lm <- lm(log1p(`Calc A`) ~ leaf_pair*pop_num, data = chem_df)

# residuals
plot(density(resid(calcA_lm))) # much better

summary(calcA_lm)

Anova(calcA_lm, type = 3) 


# Conandroside
conand_lm_notrans <- lm(`Conand` ~ leaf_pair*pop_num, data = chem_df)

# residuals
plot(density(resid(conand_lm_notrans))) # pretty right-skewed

# with log (1 + x) trans
conand_lm <- lm(log1p(`Conand`) ~ leaf_pair*pop_num, data = chem_df)

# residuals
plot(density(resid(conand_lm))) # better

summary(conand_lm)

Anova(conand_lm, type = 3)


# PPG h
PPGh_lm_notrans <- lm(`PPG h` ~ leaf_pair*pop_num, data = chem_df)

# residuals
plot(density(resid(PPGh_lm_notrans))) # a little left-skewed

# with log (1 + x) trans
PPGh_lm <- lm(log1p(`PPG h`) ~ leaf_pair*pop_num, data = chem_df)

# residuals
plot(density(resid(PPGh_lm))) # a little right-skewed

summary(PPGh_lm)

Anova(PPGh_lm, type = 3) 


# Verbascoside
verb_lm_notrans <- lm(`Verb` ~ leaf_pair*pop_num, data = chem_df)

# residuals
plot(density(resid(verb_lm_notrans))) # right-skewed

# with log (1 + x) trans
verb_lm <- lm(log1p(`Verb`) ~ leaf_pair*pop_num, data = chem_df)

# residuals
plot(density(resid(verb_lm))) # a little better

summary(verb_lm)

Anova(verb_lm, type = 3) 


# Unknown 16
unkn16_lm_notrans <- lm(`Unkn 16` ~ leaf_pair*pop_num, data = chem_df)

# residuals
plot(density(resid(unkn16_lm_notrans))) # right-skewed

# with log (1 + x) trans
unkn16_lm <- lm(log1p(`Unkn 16`) ~ leaf_pair*pop_num, data = chem_df)

# residuals
plot(density(resid(unkn16_lm))) # much better

summary(unkn16_lm)

Anova(unkn16_lm, type = 3) 


# PPG g
PPGg_lm_notrans <- lm(`PPG g` ~ leaf_pair*pop_num, data = chem_df)

# residuals
plot(density(resid(PPGg_lm_notrans))) # left-skewed

# with log (1 + x) trans
PPGg_lm <- lm(log1p(`PPG g`) ~ leaf_pair*pop_num, data = chem_df)

# residuals
plot(density(resid(PPGg_lm))) # better

summary(PPGg_lm)

Anova(PPGg_lm, type = 3) 


# Calceilarioside B
calcB_lm_notrans <- lm(`Calc B` ~ leaf_pair*pop_num, data = chem_df)

# residuals
plot(density(resid(calcB_lm_notrans))) # really right-skewed

# with log (1 + x) trans
calcB_lm <- lm(log1p(`Calc B`) ~ leaf_pair*pop_num, data = chem_df)

# residuals
plot(density(resid(calcB_lm))) # a little better

summary(calcB_lm)

Anova(calcB_lm, type = 3) 


# Mimuloside
mimulo_lm_notrans <- lm(`Mimulo` ~ leaf_pair*pop_num, data = chem_df)

# residuals
plot(density(resid(mimulo_lm_notrans))) # right-skewed

# with log (1 + x) trans
mimulo_lm <- lm(log1p(`Mimulo`) ~ leaf_pair*pop_num, data = chem_df)

# residuals
plot(density(resid(mimulo_lm))) # a little better

summary(mimulo_lm)

Anova(mimulo_lm, type = 3) 


# PPG d
PPGd_lm_notrans <- lm(`PPG d` ~ leaf_pair*pop_num, data = chem_df)

# residuals
plot(density(resid(PPGd_lm_notrans))) # really right-skewed

# with log (1 + x) trans
PPGd_lm <- lm(log1p(`PPG d`) ~ leaf_pair*pop_num, data = chem_df)

# residuals
plot(density(resid(PPGd_lm))) # a little better

summary(PPGd_lm)

Anova(PPGd_lm, type = 3) 


# PPG b
PPGb_lm_notrans <- lm(`PPG b` ~ leaf_pair*pop_num, data = chem_df)

# residuals
plot(density(resid(PPGb_lm_notrans))) # right-skewed

# with log (1 + x) trans
PPGb_lm <- lm(log1p(`PPG b`) ~ leaf_pair*pop_num, data = chem_df)

# residuals
plot(density(resid(PPGb_lm))) # a little better

summary(PPGb_lm)

Anova(PPGb_lm, type = 3)


# PPG c
PPGc_lm_notrans <- lm(`PPG c` ~ leaf_pair*pop_num, data = chem_df)

# residuals
plot(density(resid(PPGc_lm_notrans))) # really right-skewed

# with log (1 + x) trans
PPGc_lm <- lm(log1p(`PPG c`) ~ leaf_pair*pop_num, data = chem_df)

# residuals
plot(density(resid(PPGc_lm))) # a little better

summary(PPGc_lm)

Anova(PPGc_lm, type = 3)


# PPG e
PPGe_lm_notrans <- lm(`PPG e` ~ leaf_pair*pop_num, data = chem_df)

# residuals
plot(density(resid(PPGe_lm_notrans))) # really right-skewed

# with log (1 + x) trans
PPGe_lm <- lm(log1p(`PPG e`) ~ leaf_pair*pop_num, data = chem_df)

# residuals
plot(density(resid(PPGe_lm))) # a little better

summary(PPGe_lm)

Anova(PPGe_lm, type = 3) 


# PPG f
PPGf_lm_notrans <- lm(`PPG f` ~ leaf_pair*pop_num, data = chem_df)

# residuals
plot(density(resid(PPGf_lm_notrans))) # really right-skewed

# with log (1 + x) trans
PPGf_lm <- lm(log1p(`PPG f`) ~ leaf_pair*pop_num, data = chem_df)

# residuals
plot(density(resid(PPGf_lm))) # a little better

summary(PPGf_lm)

Anova(PPGf_lm, type = 3) #%>% tidy(.) %>% write_csv(., file = "PPGs/PPGf.csv")


##
# Let's see it
chem_df %>%
  mutate(range = fct_reorder(range, -latitude)) %>%
  mutate(Development = if_else(leaf_pair == 3, "Early", "Late")) %>%
  gather(., PPG, value, `PPG b`: `Unkn 16`, `Total PPGs`) %>%  # Include Total PPGs in gather
  group_by(PPG, Development, range) %>% 
  summarise(mean_value = mean(value, na.rm = TRUE),
            se_value = sd(value, na.rm = TRUE) / sqrt(n()), .groups = 'drop') %>%  # Calculate standard error
  mutate(PPG = factor(PPG, levels = c("Total PPGs", "Calc A", "Conand", "PPG h", 
                                      "Verb",  "Unkn 16", "PPG g", "Calc B", 
                                      "Mimulo", "PPG d", "PPG b", "PPG c", 
                                      "PPG e", "PPG f"))) %>%  # Set factor levels for PPG
  arrange(desc(PPG)) %>%  # Arrange panels in descending order of levels
  ggplot(aes(x = Development, y = log(mean_value + 1))) +  # Apply log transformation to mean_value
  geom_point(aes(colour = range), size = 2.5) +
  geom_line(aes(group = range, colour = range), size = 0.8) +
  geom_errorbar(aes(ymin = log(mean_value - se_value + 1), ymax = log(mean_value + se_value + 1), colour = range),
                width = 0.2) +  # Add error bars with log transformation
  facet_wrap(~PPG, scales = "free") +  # Each panel with its own x-axis and y-axis
  ng2 +  # custom
  theme(strip.background = element_blank()) +
  labs(y = "log[Concentration (mg/g) + 1]") +  # Label for y-axis
  scale_colour_manual(values = c("#66C2A5", "#FFD92F", "#E5C494"))


# aspect ratio = 12" x 7" for PDF portait



########

## supplement ####
#### Table of means and SDs for each PPG by population (Table S1) ####
# calculate means and SD for each PPG across pops and ontogeny
PPG_mean_SD <- chem_env %>%
  mutate(
    Population = case_when(
      pop_num == 195 ~ "TAB",
      pop_num == 206 ~ "SRN",
      pop_num == 201 ~ "CJC",
      TRUE ~ as.character(pop_num)
    ),
    Ontogeny = case_when(
      leaf_pair == 3 ~ "younger",
      leaf_pair == 8 ~ "older",
      TRUE ~ as.character(leaf_pair)
    )
  ) %>%
  group_by(Population, Ontogeny) %>%
  summarise(across(`PPG b`:`Total PPGs`, list(mean = mean, sd = sd), 
                   .names = "{.col}_{.fn}"), .groups = "drop")

write_csv(PPG_mean_SD, file = "PPGs/for_manuscript/PPG_mean_SD.csv")


#### Table of correlations between PPGs (Table S2) ####

# Define the desired order of PPGs
PPG_order_TableS2 <- c("Total PPGs",
                       "Calc A",
                       "Conand",
                       "Verb",
                       "Calc B",
                       "Mimulo",
                       "Unkn 16",
                       "PPG b",
                       "PPG c",
                       "PPG d",
                       "PPG e",
                       "PPG f",
                       "PPG g",
                       "PPG h")

# Select only PPG-related columns and reorder them
ppg_vars <- chem_env %>%
  select(all_of(PPG_order_TableS2))

# Compute correlation matrix
cor_matrix <- cor(ppg_vars, use = "pairwise.complete.obs", method = "pearson")

# Compute p-value matrix
p_matrix <- matrix(NA, ncol = ncol(ppg_vars), nrow = ncol(ppg_vars))
colnames(p_matrix) <- colnames(ppg_vars)
rownames(p_matrix) <- colnames(ppg_vars)

for (i in 1:ncol(ppg_vars)) {
  for (j in 1:ncol(ppg_vars)) {
    if (i < j) {  # Only compute for upper triangle
      test <- cor.test(ppg_vars[[i]], ppg_vars[[j]], method = "pearson")
      p_matrix[i, j] <- test$p.value
    }
  }
}

# Create a table with correlation coefficients below the diagonal and p-values above
formatted_table <- matrix("", ncol = ncol(ppg_vars), nrow = ncol(ppg_vars))
colnames(formatted_table) <- colnames(ppg_vars)
rownames(formatted_table) <- colnames(ppg_vars)

for (i in 1:ncol(ppg_vars)) {
  for (j in 1:ncol(ppg_vars)) {
    if (i < j) {
      formatted_table[i, j] <- sprintf("%.3f", p_matrix[i, j])  # P-values above diagonal
    } else if (i > j) {
      formatted_table[i, j] <- sprintf("%.3f", cor_matrix[i, j])  # Correlations below diagonal
    } else {
      formatted_table[i, j] <- "—"  # Diagonal empty or placeholder
    }
  }
}

# Convert to a data frame for better readability
formatted_df <- as.data.frame(formatted_table)

# Ensure the order is preserved
formatted_df <- formatted_df %>%
  rownames_to_column(var = "PPG") %>%
  arrange(factor(PPG, levels = ppg_order))  # Reorder rows

# Reorder columns as well
formatted_df <- formatted_df %>%
  select(PPG, all_of(ppg_order))

# Print the formatted table
print(formatted_df, row.names = FALSE)

# Write to CSV
write_csv(formatted_df, file = "PPGs/for_manuscript/PPG_cor_table.csv")


#### plot arrows of NMDS with ggplot (Figure S1) ####

# Define the desired order of PPGs
PPG_order_FigS1 <- c("Calc A",
                     "Conand",
                     "Verb",
                     "Calc B",
                     "Mimulo",
                     "Unkn 16",
                     "PPG b",
                     "PPG c",
                     "PPG d",
                     "PPG e",
                     "PPG f",
                     "PPG g",
                     "PPG h")

# make the df from vegan output
NMDS_sum <- cbind(PPGs = rownames(NMDS_scores),
                  NMDS1 = ef$vectors$arrows[ ,1], 
                  NMDS2 = ef$vectors$arrows[ ,2],
                  r2 = ef$vectors$r, 
                  p = ef$vectors$pvals) %>% 
  as.data.frame() %>%
  mutate(across(NMDS1:p, .fns = as.numeric)) %>%
  mutate(NMDS1 = NMDS1*sqrt(r2)) %>%
  mutate(NMDS2 = NMDS2*sqrt(r2)) %>% 
  mutate(PPGs = factor(PPGs, levels = reorder_PPGs)) %>%
  arrange(PPGs)

# the plot
ggplot(NMDS_sum, aes(x = NMDS1, y = NMDS2)) + 
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.75) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.75) +
  ng1 +
  labs(x = "NMDS1", y = "NMDS2", size = 5) +
  coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1)) +
  geom_segment(data = NMDS_sum, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(1/2, "picas")), color = "black") +
  geom_text_repel(data = NMDS_sum %>% filter(p < 0.05), aes(x = NMDS1, y = NMDS2, label = PPGs), fontface = "bold",
                  size = 5, point.padding = 1.5) +
  geom_text_repel(data = NMDS_sum %>% filter(p > 0.05), aes(x = NMDS1, y = NMDS2, label = PPGs),
                  size = 5, point.padding = 1.5) 


#### get table of how much variance is explained for each NMDS axis (Table S3) ####

# Define the desired order of PPGs
PPG_order_TableS3 <- c("Calc A",
                       "Conand",
                       "Verb",
                       "Calc B",
                       "Mimulo",
                       "Unkn 16",
                       "PPG b",
                       "PPG c",
                       "PPG d",
                       "PPG e",
                       "PPG f",
                       "PPG g",
                       "PPG h")

NMDS_sum <- cbind(PPGs = rownames(NMDS_scores),
                  NMDS1 = ef$vectors$arrows[ ,1], 
                  NMDS2 = ef$vectors$arrows[ ,2],
                  r2 = ef$vectors$r, 
                  p = ef$vectors$pvals) %>% 
  as.data.frame() %>%
  mutate(across(NMDS1:p, .fns = as.numeric)) %>%
  mutate(NMDS1 = NMDS1*sqrt(r2)) %>%
  mutate(NMDS2 = NMDS2*sqrt(r2)) %>% 
  mutate(PPGs = factor(PPGs, levels = PPG_order_TableS3)) %>%
  arrange(PPGs)

write_csv(NMDS_sum, file = "PPGs/for_manuscript/PPG_NMDS_var.csv")









#####















##