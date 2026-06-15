setwd("D:/atacama_project")
library(dplyr)
library(readxl)
library(ggplot2)
library(ggeffects)
library(phyloseq)
library(ggsci)

######################################################
##### 1. Preparation of Data for Diversity Plots #####
######################################################
all.alpha.noW <- readRDS("D:/atacama_project/R_objects/all.alpha.noW.rds")
metadata <- read_xlsx("AtacamaMetadata.xlsx")

#1.1: Add Depth to the data frame 
all.alpha.noW$Depth <- metadata$Depth[match(rownames(all.alpha.noW), metadata$SampleID)]

#1.2: Add Depth Categories 
all.alpha.noW$DepthCategory <- cut(
  all.alpha.noW$Depth,
  breaks = c(0, 5, 30, Inf),
  labels = c("Surface", "Sub-Surface", "Deep-Soil"),
  right = FALSE
)

all.metadata.noW <- all.alpha.noW
saveRDS(all.metadata.noW, "R_Objects/all.metadata.noW.rds")

all.metadata.noW$logDepth <- log(all.metadata.noW$Depth + 0.01)

all.metadata.noW$logDepth.centered <- all.metadata.noW$logDepth - mean(all.metadata.noW$logDepth)
saveRDS(all.metadata.noW, "R_Objects/all.meta.noW.log.centre.rds")
all.meta.noW.log.centre <- all.metadata.noW


############################################
##### 2. Diversity Plots Against Depth #####
############################################
##### 2.1: Shannon Diversity #####
m.shannon_depth <- lm(Shannon.Diversity ~ Depth + Dataset,
                data = all.meta.noW.log.centre)

summary(m.shannon_depth)
# Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  2.792546   0.400213   6.978 2.01e-08 ***
# Depth       -0.006566   0.010151  -0.647  0.52144    
# DatasetF     1.674814   0.541650   3.092  0.00361 ** 
# DatasetN     2.301923   0.432200   5.326 4.17e-06 ***
# DatasetS     1.667336   0.568908   2.931  0.00557 ** 

# Residual standard error: 1.047 on 40 degrees of freedom
# Multiple R-squared:  0.4195,	Adjusted R-squared:  0.3615 
# F-statistic: 7.228 on 4 and 40 DF,  p-value: 0.0001771

#EXtract adjusted R square and p-value of depth 
adj_r2 <- summary(m.shannon_depth)$adj.r.squared
p_depth <- summary(m.shannon_depth)$coefficients["Depth", "Pr(>|t|)"]

#Create label for plotting 
label <- paste0("Adjusted R² = ", sprintf("%.3f", adj_r2),
  "\nDepth p-value = ", format.pval(p_depth, digits = 3))

pdf("Results/m.shannon_depth_plot.pdf", width = 8.5, height = 10)

pred <- ggpredict(
  m.shannon_depth,
  terms = "Depth [all]"
)

ggplot() +
  geom_point(
    data = all.meta.noW.log.centre,
    aes(x = Depth, y = Shannon.Diversity, colour = Dataset),size = 2) +
  geom_ribbon(data = pred, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_line(data = pred, aes(x = x, y = predicted), colour = "black", linewidth = 1.2) +
  labs(
    title = "Relationship Between Shannon Diversity and Soil Depth",
    x = "Soil Depth",
    y = "Shannon Diversity",
    colour = "Dataset") +
  annotate("text", x = Inf, y = Inf, label = label, hjust = 1.1, vjust = 1.5, size = 4) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.background = element_rect(colour = "black", fill = "white"),
        axis.title.x = element_text(face = "bold", margin = margin(t=10)),
        axis.title.y = element_text(face = "bold", margin = margin(r=10)))

dev.off()

#Try Log Centered Transformation as well
m.shannon_depth_log <- lm(Shannon.Diversity ~ logDepth.centered + Dataset,
                      data = all.meta.noW.log.centre)
summary(m.shannon_depth_log)
# Coefficients:
#                    Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        2.72746    0.37132   7.345 6.23e-09 ***
# logDepth.centered -0.07017    0.08660  -0.810  0.42256    
# DatasetF           1.42923    0.53020   2.696  0.01023 *  
# DatasetN           2.30192    0.43093   5.342 3.96e-06 ***
# DatasetS           1.64563    0.56899   2.892  0.00616 ** 

# Residual standard error: 1.044 on 40 degrees of freedom
# Multiple R-squared:  0.4229,	Adjusted R-squared:  0.3652 
# F-statistic: 7.329 on 4 and 40 DF,  p-value: 0.0001586

pdf("Results/m.shannon_depth_log_plot.pdf", width = 8.5, height = 10)

pred <- ggpredict(
  m.shannon_depth_log,
  terms = "logDepth.centered [all]"
)

ggplot() +
  geom_point(
    data = all.meta.noW.log.centre,
    aes(x = logDepth.centered, y = Shannon.Diversity, colour = Dataset),size = 2) +
  geom_ribbon(data = pred, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_line(data = pred, aes(x = x, y = predicted), colour = "black", linewidth = 1.2) +
  labs(
    title = "Relationship Between Shannon Diversity and Soil Log Depth",
    x = "Soil Log Depth",
    y = "Shannon Diversity",
    colour = "Dataset") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.background = element_rect(colour = "black", fill = "white"),
        axis.title.x = element_text(face = "bold", margin = margin(t=10)),
        axis.title.y = element_text(face = "bold", margin = margin(r=10)))

dev.off()

#Conclusion: Spread of Data still looks the same, don't need to transform 


##### 2.2: Observed Richness ##### 
m.richness_depth <- lm(Observed.Richness ~ Depth + Dataset,
                      data = all.meta.noW.log.centre)

summary(m.richness_depth)

# Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  48.7789   168.5777   0.289    0.774    
# Depth        -0.4353     4.2760  -0.102    0.919    
# DatasetF    350.3965   228.1535   1.536    0.132    
# DatasetN    799.7045   182.0511   4.393    8e-05 ***
# DatasetS    341.1963   239.6349   1.424    0.162

# Residual standard error: 440.9 on 40 degrees of freedom
# Multiple R-squared:  0.3585,	Adjusted R-squared:  0.2944 
# F-statistic: 5.589 on 4 and 40 DF,  p-value: 0.001138

#EXtract adjusted R square and p-value of depth 
adj_r2 <- summary(m.richness_depth)$adj.r.squared
p_depth <- summary(m.richness_depth)$coefficients["Depth", "Pr(>|t|)"]

#Create label for plotting 
label <- paste0("Adjusted R² = ", sprintf("%.3f", adj_r2),
                "\nDepth p-value = ", format.pval(p_depth, digits = 3))

pdf("Results/m.richness_depth_plot.pdf", width = 8.5, height = 10)

pred <- ggpredict(
  m.richness_depth,
  terms = "Depth [all]"
)

ggplot() +
  geom_point(
    data = all.meta.noW.log.centre,
    aes(x = Depth, y = Observed.Richness, colour = Dataset),size = 2) +
  geom_ribbon(data = pred, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_line(data = pred, aes(x = x, y = predicted), colour = "black", linewidth = 1.2) +
  labs(
    title = "Relationship Between Observed Richness and Soil Depth",
    x = "Soil Depth",
    y = "Observed Richness",
    colour = "Dataset") +
  annotate("text", x = Inf, y = Inf, label = label, hjust = 1.1, vjust = 1.5, size = 4) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.background = element_rect(colour = "black", fill = "white"),
        axis.title.x = element_text(face = "bold", margin = margin(t=10)),
        axis.title.y = element_text(face = "bold", margin = margin(r=10)))

dev.off()


##### 2.3: Simpson Diversity #####
m.simpson_depth <- lm(Simpson.Diversity ~ Depth + Dataset,
                       data = all.meta.noW.log.centre)

summary(m.simpson_depth)
# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.8763737  0.0182128  48.119  < 2e-16 ***
# Depth       -0.0003228  0.0004620  -0.699 0.488726    
# DatasetF     0.0883804  0.0246493   3.586 0.000905 ***
# DatasetN     0.0845020  0.0196684   4.296 0.000108 ***
# DatasetS     0.0853757  0.0258897   3.298 0.002052 **

# Residual standard error: 0.04764 on 40 degrees of freedom
# Multiple R-squared:  0.3465,	Adjusted R-squared:  0.2812 
# F-statistic: 5.302 on 4 and 40 DF,  p-value: 0.0016

#EXtract adjusted R square and p-value of depth 
adj_r2 <- summary(m.simpson_depth)$adj.r.squared
p_depth <- summary(m.simpson_depth)$coefficients["Depth", "Pr(>|t|)"]

#Create label for plotting 
label <- paste0("Adjusted R² = ", sprintf("%.3f", adj_r2),
                "\nDepth p-value = ", format.pval(p_depth, digits = 3))

pdf("Results/m.simpson_depth_plot.pdf", width = 8.5, height = 10)

pred <- ggpredict(
  m.simpson_depth,
  terms = "Depth [all]"
)

ggplot() +
  geom_point(
    data = all.meta.noW.log.centre,
    aes(x = Depth, y = Simpson.Diversity, colour = Dataset),size = 2) +
  geom_ribbon(data = pred, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_line(data = pred, aes(x = x, y = predicted), colour = "black", linewidth = 1.2) +
  labs(
    title = "Relationship Between Simpson Diveristy and Soil Depth",
    x = "Soil Depth",
    y = "Simpson Diversity",
    colour = "Dataset") +
  annotate("text", x = Inf, y = Inf, label = label, hjust = 1.1, vjust = 1.5, size = 4) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.background = element_rect(colour = "black", fill = "white"),
        axis.title.x = element_text(face = "bold", margin = margin(t=10)),
        axis.title.y = element_text(face = "bold", margin = margin(r=10)))

dev.off()


##############################################################
##### 3. Diversity Plots Against log Moisture (Centered) #####
##############################################################
all.alpha.noW.log.centre <- readRDS("D:/atacama_project/R_objects/all.alpha.noW.log.centre.rds")

##### 3.1: Shannon Diversity #####
m.shannon_moisture <- lm(Shannon.Diversity ~ logMoisture.centered + Dataset,
                      data = all.alpha.noW.log.centre)

summary(m.shannon_moisture)
# Coefficients:
#                       Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           2.48599    0.21514  11.555 2.55e-14 ***
# logMoisture.centered  0.53256    0.05917   9.001 3.67e-11 ***
# DatasetF              2.47960    0.31142   7.962 8.89e-10 ***
# DatasetN              2.36075    0.24986   9.448 9.67e-12 ***
# DatasetS              1.66545    0.32676   5.097 8.69e-06 ***

# Residual standard error: 0.605 on 40 degrees of freedom
# Multiple R-squared:  0.8061,	Adjusted R-squared:  0.7867 
# F-statistic: 41.58 on 4 and 40 DF,  p-value: 9.628e-14

#EXtract adjusted R square and p-value of depth 
adj_r2 <- summary(m.shannon_moisture)$adj.r.squared
p_moisture <- summary(m.shannon_moisture)$coefficients["logMoisture.centered", "Pr(>|t|)"]

#Create label for plotting 
label <- paste0("Adjusted R² = ", sprintf("%.3f", adj_r2),
                "\nlog Moisture p-value = ", format.pval(p_moisture, digits = 3))

pdf("Results/m.shannon_moisture_plot.pdf", width = 8.5, height = 10)

pred <- ggpredict(
  m.shannon_moisture,
  terms = "logMoisture.centered [all]")

ggplot() +
  geom_point(
    data = all.alpha.noW.log.centre,
    aes(x = logMoisture.centered, y = Shannon.Diversity, colour = Dataset),size = 2) +
  geom_ribbon(data = pred, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_line(data = pred, aes(x = x, y = predicted), colour = "black", linewidth = 1.2) +
  labs(
    title = "Relationship Between Shannon Diversity and Soil log Moisture",
    x = "Soil log Moisture",
    y = "Shannon Diversity",
    colour = "Dataset") +
  annotate("text", x = Inf, y = Inf, label = label, hjust = 1.1, vjust = 1.5, size = 4) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.background = element_rect(colour = "black", fill = "white"),
        axis.title.x = element_text(face = "bold", margin = margin(t=10)),
        axis.title.y = element_text(face = "bold", margin = margin(r=10)))

dev.off()

#For w/o values 
pdf("Results/m.shannon_moisture_plot(wo.values).pdf", width = 8.5, height = 10)

pred <- ggpredict(
  m.shannon_moisture,
  terms = "logMoisture.centered [all]")

ggplot() +
  geom_point(
    data = all.alpha.noW.log.centre,
    aes(x = logMoisture.centered, y = Shannon.Diversity, colour = Dataset),size = 2) +
  geom_ribbon(data = pred, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_line(data = pred, aes(x = x, y = predicted), colour = "black", linewidth = 1.2) +
  labs(
    title = "Relationship Between Shannon Diversity and Soil log Moisture",
    x = "Soil log Moisture",
    y = "Shannon Diversity",
    colour = "Dataset") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.background = element_rect(colour = "black", fill = "white"),
        axis.title.x = element_text(face = "bold", margin = margin(t=10)),
        axis.title.y = element_text(face = "bold", margin = margin(r=10)))

dev.off()


##### 3.2: Observed Richness #####
m.richness_moisture <- lm(Observed.Richness ~ logMoisture.centered + Dataset,
                         data = all.alpha.noW.log.centre)

summary(m.richness_moisture)
# Coefficients:
#                       Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            -33.83     109.81  -0.308 0.759596    
# logMoisture.centered   194.74      30.20   6.448 1.11e-07 ***
# DatasetF               680.72     158.96   4.282 0.000113 ***
# DatasetN               821.22     127.54   6.439 1.14e-07 ***
# DatasetS               328.22     166.79   1.968 0.056044 .  

# Residual standard error: 308.8 on 40 degrees of freedom
# Multiple R-squared:  0.6854,	Adjusted R-squared:  0.6539 
# F-statistic: 21.79 on 4 and 40 DF,  p-value: 1.328e-09

#EXtract adjusted R square and p-value of depth 
adj_r2 <- summary(m.richness_moisture)$adj.r.squared
p_moisture <- summary(m.richness_moisture)$coefficients["logMoisture.centered", "Pr(>|t|)"]

#Create label for plotting 
label <- paste0("Adjusted R² = ", sprintf("%.3f", adj_r2),
                "\nMoisture p-value = ", format.pval(p_moisture, digits = 3))

pdf("Results/m.richness_moisture_plot.pdf", width = 8.5, height = 10)

pred <- ggpredict(
  m.richness_moisture,
  terms = "logMoisture.centered [all]")

ggplot() +
  geom_point(
    data = all.alpha.noW.log.centre,
    aes(x = logMoisture.centered, y = Observed.Richness, colour = Dataset),size = 2) +
  geom_ribbon(data = pred, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_line(data = pred, aes(x = x, y = predicted), colour = "black", linewidth = 1.2) +
  labs(
    title = "Relationship Between Observed Richness and Soil log Moisture",
    x = "Soil log Moisture",
    y = "Observed Richness",
    colour = "Dataset") +
  annotate("text", x = Inf, y = Inf, label = label, hjust = 1.1, vjust = 1.5, size = 4) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.background = element_rect(colour = "black", fill = "white"),
        axis.title.x = element_text(face = "bold", margin = margin(t=10)),
        axis.title.y = element_text(face = "bold", margin = margin(r=10)))

dev.off()

#For w/o values 
pdf("Results/m.richness_moisture_plot(wo.values).pdf", width = 8.5, height = 10)

pred <- ggpredict(
  m.richness_moisture,
  terms = "logMoisture.centered [all]")

ggplot() +
  geom_point(
    data = all.alpha.noW.log.centre,
    aes(x = logMoisture.centered, y = Observed.Richness, colour = Dataset),size = 2) +
  geom_ribbon(data = pred, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_line(data = pred, aes(x = x, y = predicted), colour = "black", linewidth = 1.2) +
  labs(
    title = "Relationship Between Observed Richness and Soil log Moisture",
    x = "Soil log Moisture",
    y = "Observed Richness",
    colour = "Dataset") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.background = element_rect(colour = "black", fill = "white"),
        axis.title.x = element_text(face = "bold", margin = margin(t=10)),
        axis.title.y = element_text(face = "bold", margin = margin(r=10)))

dev.off()

##### 3.3: Simpson Index #####
m.simpson_moisture <- lm(Simpson.Diversity ~ logMoisture.centered + Dataset,
                          data = all.alpha.noW.log.centre)

summary(m.simpson_moisture)
#Coefficients:
#                       Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          0.863425   0.012173  70.927  < 2e-16 ***
# logMoisture.centered 0.020749   0.003348   6.198 2.49e-07 ***
# DatasetF             0.118508   0.017621   6.725 4.53e-08 ***
# DatasetN             0.086794   0.014138   6.139 3.01e-07 ***
# DatasetS             0.085721   0.018490   4.636 3.75e-05 ***  

# Residual standard error: 0.03423 on 40 degrees of freedom
# Multiple R-squared:  0.6626,	Adjusted R-squared:  0.6288 
# F-statistic: 19.64 on 4 and 40 DF,  p-value: 5.22e-09

#EXtract adjusted R square and p-value of depth 
adj_r2 <- summary(m.simpson_moisture)$adj.r.squared
p_moisture <- summary(m.simpson_moisture)$coefficients["logMoisture.centered", "Pr(>|t|)"]

#Create label for plotting 
label <- paste0("Adjusted R² = ", sprintf("%.3f", adj_r2),
                "\nMoisture p-value = ", format.pval(p_moisture, digits = 3))

pdf("Results/m.simpson_moisture_plot.pdf", width = 8.5, height = 10)

pred <- ggpredict(
  m.simpson_moisture,
  terms = "logMoisture.centered [all]")

ggplot() +
  geom_point(
    data = all.alpha.noW.log.centre,
    aes(x = logMoisture.centered, y = Simpson.Diversity, colour = Dataset),size = 2) +
  geom_ribbon(data = pred, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_line(data = pred, aes(x = x, y = predicted), colour = "black", linewidth = 1.2) +
  labs(
    title = "Relationship Between Simpson Diversity and Soil log Moisture",
    x = "Soil log Moisture",
    y = "Observed Richness",
    colour = "Dataset") +
  annotate("text", x = Inf, y = Inf, label = label, hjust = 1.1, vjust = 1.5, size = 4) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.background = element_rect(colour = "black", fill = "white"),
        axis.title.x = element_text(face = "bold", margin = margin(t=10)),
        axis.title.y = element_text(face = "bold", margin = margin(r=10)))

dev.off()

#wo values
pdf("Results/m.simpson_moisture_plot(wo.values).pdf", width = 8.5, height = 10)

pred <- ggpredict(
  m.simpson_moisture,
  terms = "logMoisture.centered [all]")

ggplot() +
  geom_point(
    data = all.alpha.noW.log.centre,
    aes(x = logMoisture.centered, y = Simpson.Diversity, colour = Dataset),size = 2) +
  geom_ribbon(data = pred, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_line(data = pred, aes(x = x, y = predicted), colour = "black", linewidth = 1.2) +
  labs(
    title = "Relationship Between Simpson Diversity and Soil log Moisture",
    x = "Soil log Moisture",
    y = "Observed Richness",
    colour = "Dataset") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.background = element_rect(colour = "black", fill = "white"),
        axis.title.x = element_text(face = "bold", margin = margin(t=10)),
        axis.title.y = element_text(face = "bold", margin = margin(r=10)))

dev.off()


###############################################
##### 4. Data Preparation for Phyla Plots #####
###############################################
rare.min.no.contam <- readRDS("R_Objects/rare.min.no.contam.rds")
rare.min.no.contam.noW <- rare.min.no.contam[names(rare.min.no.contam) != "W"]

###### Combine the ps objects ######
combined.phy <- do.call(merge_phyloseq, rare.min.no.contam.noW)

###### Add Moisture Groups ######
sample_data(combined.phy)$MoistureGroup <- cut(
  sample_data(combined.phy)$Moisture,
  breaks = c(-Inf, 2.5, 7.5, Inf),
  labels = c("Low", "Medium", "High"),
  right = FALSE
)

##### Relative Abundance ######
phy.rel <- transform_sample_counts(combined.phy, function(x) x / sum(x))

##### Aggregate at Phylum Level ######
phy.phylum <- tax_glom(phy.rel, taxrank = "Phylum")
phy.phylum.df <- psmelt(phy.phylum)

##### Handle Unassigned ######
phy.phylum.df$Phylum <- ifelse(
  is.na(phy.phylum.df$Phylum) |
    grepl("Unassigned", phy.phylum.df$Phylum),
  "Unassigned",
  as.character(phy.phylum.df$Phylum)
)

##### Top 10 Phyla Overall ######
top10.phyla <- phy.phylum.df %>%
  group_by(Phylum) %>%
  summarise(TotalAbundance = sum(Abundance), .groups = "drop") %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 10) %>%
  pull(Phylum)

##### Create 'Other' Category ######
phy.phylum.df$Phylum <- ifelse(phy.phylum.df$Phylum %in% c(top10.phyla, "Unassigned"),
                               phy.phylum.df$Phylum, "Other")

##### Mean Relative Abundance ######
plot.df <- phy.phylum.df %>%
  group_by(MoistureGroup, Phylum) %>%
  summarise(Abundance = mean(Abundance), .groups = "drop")

##### Force Bars to Sum to 100% ######
plot.df <- plot.df %>%
  group_by(MoistureGroup) %>%
  mutate(Abundance = Abundance / sum(Abundance)) %>%
  ungroup()

##### Put Other and Unassigned at Bottom, Rest by Abundance #####
phylum.order <- plot.df %>%
  group_by(Phylum) %>%
  summarise(TotalAbundance = sum(Abundance)) %>%
  arrange(desc(TotalAbundance)) %>%
  pull(Phylum)

phylum.order <- c(setdiff(phylum.order, c("Other", "Unassigned")), 
                  "Other", "Unassigned")

plot.df$Phylum <- factor(plot.df$Phylum, levels = phylum.order)


##################################################
##### 5. Top 10 Phyla Stacked Bar Plot ###########
##################################################
p.phylum <- ggplot(plot.df, aes(x = MoistureGroup, y = Abundance * 100, fill = Phylum)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_d3(palette = "category20") +
  theme_bw() +
  theme(axis.title.x = element_text(face = "bold", margin = margin(t=10)),
        axis.title.y = element_text(face = "bold", margin = margin(r=10)), 
        axis.text = element_text(colour = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold",),
        legend.background = element_rect(colour = "black", fill = "white", linewidth = 0.3),
        legend.title = element_text(face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(
    title = "Top 10 Phyla Across Moisture Groups",
    x = "Moisture Group",
    y = "Mean Relative Abundance (%)")

pdf("Results/Top10Phylum_MoistureGroups.pdf", width = 8.5, height = 10)

p.phylum

dev.off()


#############################################
##### 6. Top 10 Phyla Bubble Plot ###########
#############################################
##### Create Bubble Plot Data Frame ######
bubble.df <- phy.phylum.df %>%
  filter(Phylum %in% top10.phyla) %>%
  group_by(MoistureGroup, Phylum) %>%
  summarise(MeanAbundance = mean(Abundance), .groups = "drop")

##### Order Phyla by Overall Abundance ######
bubble.df$Phylum <- factor(bubble.df$Phylum, levels = rev(phylum.order))

##### Plot Bubble Plot ######
pdf("Results/Top10Phylum_bubble_plot.pdf", width = 8.5, height = 10)

ggplot(
  bubble.df,
  aes(x = MoistureGroup,y = Phylum, size = MeanAbundance * 100)) +
  geom_point(colour = "#1F4E79", alpha = 0.8) +
  scale_size_continuous(name = "Mean Abundance (%)", range = c(1, 25)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(colour = "black", fill = "white", linewidth = 0.5),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(size=12),
        axis.title.x = element_text(face = "bold", margin = margin(t=10)),
        axis.title.y = element_text(face = "bold", margin = margin(r=10)),
        axis.text = element_text(colour = "black")) +
  labs(title = "Top 10 Phyla Within Each Moisture Group",
       x = "Moisture Group",
       y = "Phylum")

dev.off()


##########################################
##### 7. Top 10 Phyla Heat Map ###########
##########################################
plot.df$Phylum <- factor(plot.df$Phylum, levels = rev(phylum.order))

p.heatmap <- ggplot(plot.df,
                    aes(x = MoistureGroup, y = Phylum, fill = Abundance * 100)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  scale_fill_gradient(low = "white", high ="green", name = "Mean \nAbundance (%)") +
  theme_bw() +
  theme(strip.placement = "outside", 
        legend.background = element_rect(colour = "black",fill = "white", linewidth = 0.3), 
        legend.title = element_text(size = 9, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 9),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.3, "cm"),
        plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b = 10)),
        strip.text.y = element_text(face = "bold", angle = 0),
        panel.grid = element_blank(), 
        axis.title.x = element_text(face = "bold", margin = margin(t = 10)),
        axis.title.y = element_text(face = "bold", margin = margin(r = 10)), 
        axis.text.y = element_text(colour = "black"), 
        axis.text.x = element_text(colour = "black")) +
  labs(
    title = "Top 10 Phyla Across Moisture Groups",
    x = "Moisture Group",
    y = "Phylum")

pdf("Results/Top10Phylum_HeatMap.pdf", width = 8.5, height = 10)

p.heatmap

dev.off()


###############################################
##### 8. Data Preparation for Genus Plots #####
###############################################
##### Aggregate at Genus Level #####
phy.genus <- tax_glom(phy.rel, taxrank = "Genus")
phy.genus.df <- psmelt(phy.genus)

##### Handle Unassigned #####
#For Genus
phy.genus.df$Genus <- ifelse(
  is.na(phy.genus.df$Genus) | grepl("Unassigned", phy.genus.df$Genus),
  "Unassigned", as.character(phy.genus.df$Genus))

#For Phylum
phy.genus.df$Phylum <- ifelse(
  is.na(phy.genus.df$Phylum) | grepl("Unassigned", phy.genus.df$Phylum),
  "Unassigned", as.character(phy.genus.df$Phylum))

#Keep the phylum information
genus.phylum <- phy.genus.df %>%
  select(Genus, Phylum) %>%
  distinct()

##### Get Top 10 Actual Genera #####
top10.genus <- phy.genus.df %>%
  filter(Genus != "Unassigned") %>%
  group_by(Genus) %>%
  summarise(TotalAbundance = sum(Abundance), .groups = "drop") %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 10) %>%
  pull(Genus)

##### Create 'Other' Category #####
phy.genus.df$Genus <- ifelse(phy.genus.df$Genus %in% c(top10.genus, "Unassigned"),
  phy.genus.df$Genus, "Other")

##### Mean Relative Abundance #####
plot.df.genus <- phy.genus.df %>%
  group_by(MoistureGroup, Genus) %>%
  summarise(Abundance = mean(Abundance), .groups = "drop")

##### Add Phylum Back #####
genus.phylum <- phy.genus.df %>%
  select(Genus, Phylum) %>%
  distinct(Genus, .keep_all = TRUE)

plot.df.genus <- left_join(
  plot.df.genus,
  genus.phylum,
  by = "Genus"
)

##### Remove Other and Unassigned #####
plot.df.genus <- plot.df.genus %>%
  filter(!Genus %in% c("Other", "Unassigned"))

##### Order Genus According to Abundance #####
genus.order <- plot.df.genus %>%
  group_by(Genus) %>%
  summarise(TotalAbundance = sum(Abundance), .groups = "drop") %>%
  arrange(desc(TotalAbundance)) %>%
  pull(Genus)


####################################
##### 9. Top 10 Genus Heat Map #####
####################################
plot.df.genus$Genus <- factor(plot.df.genus$Genus, levels = rev(genus.order))

p.heatmap <- ggplot(
  plot.df.genus, aes(x = MoistureGroup, y = Genus, fill = Abundance * 100)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  facet_grid(Phylum ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient(low = "white", high = "red",
                      name = "Mean \nAbundance (%)") +
  theme_bw() +
  theme(strip.placement = "outside", 
        legend.background = element_rect(colour = "black", fill = "white",
                                         linewidth = 0.3), 
        legend.title = element_text(size = 9, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 9),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.3, "cm"),
        plot.title = element_text(hjust = 0.5, face = "bold",
                                  margin = margin(b = 10)),
        strip.text.y = element_text(face = "bold", angle = 0),
        panel.grid = element_blank(), 
        axis.title.x = element_text(face = "bold",
                                    margin = margin(t = 15)),
        axis.title.y = element_text(face = "bold",
                                    margin = margin(r = 15)), 
        axis.text.y = element_text(size = 10, face = "italic",
                                   colour = "black"), 
        axis.text.x = element_text(colour = "black")) +
  labs(
    title = "Top 10 Genera Across Moisture Groups",
    x = "Moisture Group",
    y = "Genus")

pdf("Results/Top10Genus_Heatmap.pdf", width = 8.5, height = 10)

p.heatmap

dev.off()

#######################################################
##### 10. Top 10 Genus Faceted by Phylum Bar Plot #####
#######################################################
plot.df.genus$Genus <- factor(plot.df.genus$Genus, levels = genus.order)

p.genus.faceted <- ggplot(plot.df.genus,
                          aes(x = MoistureGroup, y = Abundance * 100, fill = Genus)) +
  geom_bar(stat = "identity", width = 0.8) +
  facet_wrap(~ Phylum, scales = "fixed", nrow = 1) +
  scale_fill_d3(palette = "category20") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        strip.text = element_text(face = "bold", size = 12),
        legend.title = element_text(hjust = 0.5, face = "bold"),
        legend.background = element_rect(colour = "black",fill = "white", linewidth = 0.3),
        legend.text = element_text(face = "italic"),
        axis.text = element_text(colour = "black"),
        axis.title.x = element_text(face = "bold", margin = margin(t = 10)),
        axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
        panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  labs(
    title = "Top 10 Genera Across Moisture Groups",
    x = "Moisture Group",
    y = "Mean Relative Abundance (%)",
    fill = "Genus")

#Save Plot 
pdf("Results/Top10Genus_FacetedByPhylum.pdf", width = 8.5, height = 10)

p.genus.faceted

dev.off()

########################################
##### 11. Top 10 Genus Bubble Plot #####
########################################
##### Ordering of Phyla and Genus by Abundance #####
phylum.order <- plot.df.genus %>%
  group_by(Phylum) %>%
  summarise(TotalAbundance = sum(Abundance), .groups = "drop") %>%
  arrange(desc(TotalAbundance)) %>%
  pull(Phylum)

plot.df.genus$Phylum <- factor(plot.df.genus$Phylum, levels = phylum.order)

genus.order <- plot.df.genus %>%
  group_by(Phylum, Genus) %>%
  summarise(TotalAbundance = sum(Abundance), .groups = "drop") %>%
  arrange(factor(Phylum, levels = phylum.order), desc(TotalAbundance)) %>%
  pull(Genus)

plot.df.genus$Genus <- factor(plot.df.genus$Genus, levels = rev(genus.order))

pdf("Results/Top10Genus_bubble_plot.pdf", width = 8.5, height = 10)

ggplot(plot.df.genus,
  aes(x = MoistureGroup, y = Genus, size = Abundance * 100, fill = Phylum)) +
  geom_point(shape = 21, colour = "black", alpha = 0.8) +
  scale_size_continuous(name = "Mean Abundance (%)", range = c(2,25)) +
  scale_fill_d3(palette = "category20") +
  theme_bw() +
  labs(
    title = "Top 10 Genera Across Moisture Groups",
    x = "Moisture Group",
    y = "Genus",
    fill = "Phylum"
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.background = element_rect(colour = "black", fill = "white"), 
    legend.title = element_text(face = "bold", hjust = 0.5), 
    axis.text.y = element_text(face = "italic", colour = "black"), 
    axis.text.x = element_text(colour = "black"), 
    axis.title.x = element_text(face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", margin = margin(r = 10))
    )

dev.off()
