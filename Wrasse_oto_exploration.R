# Lauren wrasse experiment draft code
# length weight from Cefas Silva et al Silva J. F., Ellis J. R. and Ayers R. A. 2013.
#Length-weight relationships of marine fish collected from around the
#British Isles. Sci. Ser. Tech. Rep., Cefas Lowestoft, 150: 109 pp

library(plyr)
library(tidyverse)
library(patchwork)
library(ggsci)
library(ggrepel)
library(lme4)
library(INLA)
library(gridExtra)
library(cowplot)


color_palette <- colorRampPalette(c("darkturquoise", "tomato"))(6)

color_common <-c ( "#999999", "#d55e00", "#009e73", "#cc79a7", "#e69f00", "#56b4e9")
color_species <-c ("#56b4e9", "#e69f00","#d55e00","#cc79a7", "#999999",  "#009e73")
# Baillons S.bailloni #999999
# Ballan L.bergylta #d55e00
# Corkwing S.melops #009e73
# Cuckoo L.mixtus #cc79a7
# Goldsinny C.rupestris #e69f00
# Rock cook C.exoletus #56b4e9

# Rock cook C.exoletus #56b4e9
# Goldsinny C.rupestris #e69f00
# Ballan L.bergylta #d55e00
# Cuckoo L.mixtus #cc79a7
# Baillons S.bailloni #999999
# Corkwing S.melops #009e73


#Boltzmann constant
K<-8.62*10^-5

#### z-score function 

Z.score <- function(data, colname){
  # Get mean and SD
  mu <- mean(data, na.rm = TRUE)
  st_dev <- sd(data, na.rm = TRUE)
  Z_score <- data.frame()
  for(i in 1:length(data)){
    z <- (data[i] - mu) / st_dev
    Z_score <- rbind(Z_score, z)
  }
  colnames(Z_score) <- paste("z", colname, sep = "_")
  return(Z_score)
}

Wrasse_data<-read.csv("Wrasse_Data_Clive.csv")
Dorset_data<-Wrasse_data%>% filter(Location=="Dorset")
Skye_data<-Wrasse_data%>% filter(Location=="Skye")

d18O_data<-read.csv("LOCEAN.csv")

Dorset_d18O <- d18O_data %>%
  filter(between(lat, 49, 52),
         between(lon, -3, 1),
         between(salinity, 30, 36))
Dorset_d18Ow<-mean(Dorset_d18O$d18O)


Skye_d18O <- d18O_data %>%
  filter(between(lat, 56, 59),
         between(lon, -4.5, -1.5),
         between(salinity, 30, 36))
Skye_d18O

Skye_d18Ow<-mean(Skye_d18O$d18O)


#Yuen_2023_Ballan<-read.csv("/Users/trueman/Desktop/R scripts/OtoMet/Lauren_wrasse/Yuen Ballan met data.csv", header=TRUE)
names(Wrasse_data)
str(Wrasse_data)

Wrasse_data$log10Mass<-log(Wrasse_data$Mass,10)
Wrasse_data$lnMass<-log(Wrasse_data$Mass)


plot1<-ggplot(Wrasse_data, aes(d18Ooto, d13Coto, color = Species))+
  #scale_color_npg()+
  geom_point(size=3)+
  theme_classic()+
  scale_color_manual(values = color_species)

plot1+facet_wrap( ~ Location, ncol=2)



d18Ow<-0 # global average 

d13C_DIC<-1 #first guess water / DIC_ Burt et al 1990. d13C DIC North Sea


# temp 
Wrasse_data <- Wrasse_data %>%
  mutate(
    Temp = if_else(
      Location == "Dorset",
      (d18Ooto - Dorset_d18Ow - 3.465) / -0.209,
      (d18Ooto - Skye_d18Ow   - 3.465) / -0.209
    )
  )


#Wrasse_data$Temp<-((Wrasse_data$d18Ooto-d18Ow)-3.465)/-0.209# update Morissette et al 2023 general marine fishes


Wrasse_data$InverseTemp<-1/((K)*(Wrasse_data$Temp+273))

# Cresp

Wrasse_data$Cresp<-(Wrasse_data$d13Coto-Wrasse_data$d13CDIC) / (Wrasse_data$d13Cdiet- Wrasse_data$d13CDIC)# not lipid corrected - Chung et al 2019

## convert Cresp to mg O2/Kg/hr Using linear mean equation from Chung et al 2019  Cresp<- 4.01+(0.000971*O2)  (in mg O2/Kg/hr) 
Wrasse_data$O2_consumption<-(Wrasse_data$Cresp-0.041)/0.000971
Wrasse_data$Org_O2<-Wrasse_data$O2_consumption*(Wrasse_data$Mass/1000)# going from mass specific Oc consumption to whole organism
Wrasse_data$ln_O2_consumption<-log(Wrasse_data$O2_consumption)

summary_FMR_output<-Wrasse_data%>%
  dplyr::group_by(Common_name, Location)%>%
  dplyr::summarise(
 number = dplyr::n(),
    d13Cotomed =  median(d13Coto, na.rm = TRUE),
    d13CotoIQR =  IQR(d13Coto, na.rm = TRUE),
    d18Ootomed =  median(d18Ooto, na.rm = TRUE),
    d18OotoIQR =  IQR(d18Ooto, na.rm = TRUE),
    FMR_median=median(O2_consumption, na.rm = TRUE),
 FMR_min=min(O2_consumption, na.rm = TRUE),
 FMR_max=max(O2_consumption, na.rm = TRUE),
  #  FMR_IQR=IQR(O2_consumption, na.rm = TRUE),
    Temp_median=median(Temp, na.rm = TRUE),
 Temp_min=min(Temp, na.rm = TRUE),
 Temp_max=max(Temp, na.rm = TRUE)
   # TempIQR=IQR(Temp, na.rm = TRUE)
  )


summary_FMR_output_location<-Wrasse_data%>%
  dplyr::group_by(Location)%>%
  dplyr::summarise(
    #   number = dplyr::n(),
    FMR_min=min(O2_consumption, na.rm = TRUE),
    FMR_median=median(O2_consumption, na.rm = TRUE),
    FMR_max=max(O2_consumption, na.rm = TRUE),
    Temp_min=min(Temp, na.rm = TRUE),
    Temp_median=median(Temp, na.rm = TRUE),
    Temp=max(Temp, na.rm = TRUE)
  )

ln_median_o2 <- Wrasse_data %>%
  dplyr::group_by(Species) %>%
  dplyr::summarise(
    ln_median_O2_consumption = median(ln_O2_consumption, na.rm = TRUE)
  )

median_InvTemp <- Wrasse_data %>%
  dplyr::group_by(Species) %>%
  dplyr::summarise(
    median_Inv_Temp = median(InverseTemp, na.rm = TRUE)
  )

median_lnMass <- Wrasse_data %>%
  dplyr::group_by(Species) %>%
  dplyr::summarise(
    median_lnMass = median(lnMass, na.rm = TRUE)
  )


combined_table <- Wrasse_data %>%
  left_join(ln_median_o2, by = "Species")%>%
  left_join(median_InvTemp, by = "Species")%>%
  left_join(median_lnMass, by = "Species")
 




### compare to literature respirometry
Lit_wrasse<-read.csv("Lab_wrasse.csv")

# extract to compare with published SMR

Goldsinny_FMR <- combined_table %>%
  filter(Species == "Ctenolabrus_rupestris" & Temp >= 8.5 & Temp <= 11.5) %>%
  dplyr::select(Mass, Temp, O2_consumption)

Ballan_FMR <- combined_table %>%
  filter(Species == "Labrus_bergylta" & Mass >= 200 & Mass <= 350) %>%
  dplyr::select(Mass, Temp, O2_consumption)

Corkwing_FMR <- combined_table %>%
  filter(Species == "Symphodus_melops" & Temp >= 9 & Mass >= 8 & Mass <= 20) %>%
  dplyr::select(Mass, Temp, O2_consumption)

mean(Goldsinny_FMR$Mass)
mean(Goldsinny_FMR$O2_consumption)

sd(Goldsinny_FMR$Mass)
sd(Goldsinny_FMR$O2_consumption)


## set threshold for 'largest' section to restrict for INLA modellingm (in grammes)




Lab_field_ballan <- ggplot() +
#  geom_line(data = Lit_wrasse[Lit_wrasse$Species=="Labrus_bergylta",], aes(x = Temp, y = SMR), lwd=1) +
  geom_smooth(data = Lit_wrasse[Lit_wrasse$Species=="Labrus_bergylta",], aes(x = Temp, y = MMR), method="loess", lwd=1) +
  geom_smooth(data = Lit_wrasse[Lit_wrasse$Species=="Labrus_bergylta",], aes(x = Temp, y = SMR), method="loess", lwd=1)+
 # geom_point(data = Goldsinny_FMR, aes(x = Temp, y = O2_consumption, color="green"), shape = 16) +
  geom_point(data = Ballan_FMR, aes(x = Temp, y = O2_consumption), color="#d55e00", shape = 16, cex=3) +
#  geom_point(data = Corkwing_FMR, aes(x = Temp, y = O2_consumption, color="orange"), shape = 16) +
  
  
 # geom_point(data = combined_table, aes(x = Temp, y = SMR, color = Species), shape = 17) +
  labs(
       x = "Temperature (°C)",
       y = expression(mg~O[2]~kg^{-1}~h^{-1})) +
  ggtitle("L. bergylta") +
  theme(
    panel.background = element_rect(colour = "black", fill = 'transparent'), # transparent panel bg
    plot.background = element_rect(fill = 'transparent', color = NA), # transparent plot bg
    panel.grid.major = element_blank(), # remove major gridlines
    panel.grid.minor = element_blank(), # remove minor gridlines
    axis.text.y = element_text(size = 18, colour = "black"),
    axis.text.x = element_text(size = 12, colour = "black"), 
    axis.title.x = element_text(size = 20, colour = "black"),
    axis.title.y = element_text(size = 20, colour = "black"),
    plot.margin = unit(c(1, 1, 0, 1), "lines") # Adjust margins for spacing
  )

#create summary table for error bars for corkwing plot
# SMR summary
SMR_sum <- Lit_wrasse %>%
  filter(Species == "Symphodus_melops") %>%
  group_by(Temp) %>%
  summarise(
    mean = mean(SMR, na.rm = TRUE),
    se   = sd(SMR, na.rm = TRUE) / sqrt(n())
  )
# MMR summary
MMR_sum <- Lit_wrasse %>%
  filter(Species == "Symphodus_melops") %>%
  group_by(Temp) %>%
  summarise(
    mean = mean(MMR, na.rm = TRUE),
    se   = sd(MMR, na.rm = TRUE) / sqrt(n())
  )

# corkwing  - Norin single temp data
Lab_field_corkwing <- ggplot() +
  geom_pointrange(data = SMR_sum, 
              aes(x = Temp, y = mean - 0.2, ymin = 119, ymax = 270), linewidth=0.8, width = 0.1, colour = "blue") +  
  geom_pointrange(data = MMR_sum, 
              aes(x = Temp, y= mean + 0.2, ymin = 563, ymax = 940), linewidth=0.8, width = 0.1, colour = "blue") +  
  #geom_jitter(data = Lit_wrasse[Lit_wrasse$Species=="Symphodus_melops",], 
              #aes(x = Temp, y = SMR), lwd=2, width = 0.1) +  # Add jitter with width for x-axis
  #geom_jitter(data = Lit_wrasse[Lit_wrasse$Species=="Symphodus_melops",], 
              #aes(x = Temp, y = MMR), lwd=2, lty=2, width = 0.1) +  # Add jitter with width for x-axis
  
  # geom_point(data = Goldsinny_FMR, aes(x = Temp, y = O2_consumption, color="green"), shape = 16) +
 geom_point(data = Corkwing_FMR, aes(x = Temp, y = O2_consumption), color="#009e73", shape = 16, cex=3) +
  #  geom_point(data = Corkwing_FMR, aes(x = Temp, y = O2_consumption, color="orange"), shape = 16) +
  
  
  # geom_point(data = combined_table, aes(x = Temp, y = SMR, color = Species), shape = 17) +
  labs(title = "S. melops",
       x = "Temperature (°C)",
       y = expression(mg~O[2]~kg^{-1}~h^{-1})) +
  theme(
    panel.background = element_rect(colour = "black", fill = 'transparent'), # transparent panel bg
    plot.background = element_rect(fill = 'transparent', color = NA), # transparent plot bg
    panel.grid.major = element_blank(), # remove major gridlines
    panel.grid.minor = element_blank(), # remove minor gridlines
    axis.text.y = element_text(size = 18, colour = "black"),
    axis.text.x = element_text(size = 12, colour = "black"), 
    axis.title.x = element_text(size = 20, colour = "black"),
    axis.title.y = element_text(size = 20, colour = "black"),
    plot.margin = unit(c(1, 1, 0, 1), "lines") # Adjust margins for spacing
  )

# goldsinny  - Perry, Sayer Davenport
Lab_field_goldsinny <- ggplot() +
  geom_smooth(data = Lit_wrasse[Lit_wrasse$Species=="Ctenolabrus_rupestris",], aes(x = Temp, y = SMR), method="loess", lwd=1) +
  geom_smooth(data = Lit_wrasse[Lit_wrasse$Species=="Ctenolabrus_rupestris",], aes(x = Temp, y = MMR), method="lm",lwd=1) +
  
  geom_point(data = Goldsinny_FMR, aes(x = Temp, y = O2_consumption), color="#e69f00", shape = 16, cex=3) +

  
  labs(title = "C. rupestris",
       x = "Temperature (°C)",
       y = expression(mg~O[2]~kg^{-1}~h^{-1})) +
  theme(
    panel.background = element_rect(colour = "black", fill = 'transparent'), # transparent panel bg
    plot.background = element_rect(fill = 'transparent', color = NA), # transparent plot bg
    panel.grid.major = element_blank(), # remove major gridlines
    panel.grid.minor = element_blank(), # remove minor gridlines
    axis.text.y = element_text(size = 18, colour = "black"),
    axis.text.x = element_text(size = 12, colour = "black"), 
    axis.title.x = element_text(size = 20, colour = "black"),
    axis.title.y = element_text(size = 20, colour = "black"),
    plot.margin = unit(c(1, 1, 0, 1), "lines") # Adjust margins for spacing
  )

# combine plots
Lab_field_plot <- plot_grid(Lab_field_ballan, Lab_field_corkwing, Lab_field_goldsinny, ncol = 3, rel_heights = c(1, 1))

# Display the plot
print(Lab_field_plot)



#MTE prediction 
# Define constants
B0 <- 4.6  # Example normalization constant
a <- -0.25  # Example scaling exponent
E <- -0.65  # Example activation energy (eV)
k <- 8.617e-5  # Boltzmann constant (eV/K)


combined_table <- combined_table %>%
  mutate(
    # Centered predictors
    M_scaled = lnMass - median_lnMass,
    T_scaled = InverseTemp - median_Inv_Temp,
    
    # Correct MTE equation
    ln_MTE_FMR = ln_median_O2_consumption +
      a * M_scaled +
      E * T_scaled,
    
    MTE_FMR = exp(ln_MTE_FMR)
  )


## plot mass against FMR by species  - overlay MTE predicted temp and mass effects
plot2 <- ggplot(combined_table, aes(Mass, O2_consumption, color = Species)) +
  geom_point(size = 1) +
  geom_smooth(method = "gam", lwd = 1) +
  scale_color_manual(values = color_species)+

  geom_smooth(aes(Mass, MTE_FMR), method = "gam", color = "black", lwd=1) + 
  labs(
       x = "Mass (g)",
       y = "FMR (mgO2 Kg-1 h-1") +
  theme(
    panel.background = element_rect(colour = "black", fill = 'transparent'), # transparent panel bg
    plot.background = element_rect(fill = 'transparent', color = NA), # transparent plot bg
    panel.grid.major = element_blank(), # remove major gridlines
    panel.grid.minor = element_blank(), # remove minor gridlines
    axis.text.y = element_text(size = 18, colour = "black"),
    axis.text.x = element_text(size = 12, colour = "black"), 
    axis.title.x = element_text(size = 20, colour = "black"),
    axis.title.y = element_text(size = 20, colour = "black"),
    plot.margin = unit(c(1, 1, 0, 1), "lines") # Adjust margins for spacing
  )+
  facet_wrap(~ Species, ncol = 2, scales = "free")

plot2



#### attempt to curve fit

# Define constants
k <- 8.617e-5  # Boltzmann constant (eV/K)


# Perform nonlinear least squares fitting


fit_parameters_by_group <- function(data) {
  start_values <- list(a = -0.35, E = -0.65)  # Initial guesses: strict MTE canonical values
  
  # Create an empty list to store results
  fit_results <- list()
  
  # Group the data by Species
  grouped_data <- data %>%
    group_by(Species) %>%
    nest()
  
  # Loop over each group
  for (i in 1:nrow(grouped_data)) {
    species <- grouped_data$Species[i]       # Extract the species name
    group_data <- grouped_data$data[[i]]     # Get the data for the current species
    
    # Perform the nonlinear least squares fitting
    fit <- nls(
      ln_O2_consumption ~ ln_median_O2_consumption + (lnMass - median(lnMass)) * a + 
        E * (1 / InverseTemp),
      data = group_data,
      start = start_values,
      algorithm = "port"
    )
    
    # Get the summary of the fit
    fit_summary <- broom::tidy(fit)
    
    # Add the Species name to the results
    fit_summary <- fit_summary %>%
      mutate(Species = species)
    
    # Append the results to the list
    fit_results[[i]] <- fit_summary
  }
  
  # Combine all results into a single tibble
  fit_results <- bind_rows(fit_results)
  
  return(fit_results)
}


# Fit parameters a and E to the data by Species
grouped_fit_results <-fit_parameters_by_group(combined_table)


ungrouped_results <- grouped_fit_results %>%
  ungroup()

# Extract and reshape estimates for `a` and `E`
fitted_parameters <- ungrouped_results %>%
  pivot_wider(
    id_cols = Species,                # Species as unique identifier
    names_from = term,                # Create columns 'a' and 'E'
    values_from = c(estimate, std.error)            # Fill with 'estimate' values
  )




# View the result
print(fitted_parameters)



combined_table <- combined_table %>%
  left_join(fitted_parameters, by = "Species")


combined_table <- combined_table%>%
  mutate(
    M_scaled = lnMass - median_lnMass,
    ln_MTE_FMR_fitted = ln_median_O2_consumption + (M_scaled * estimate_a) + (estimate_E * (1 / InverseTemp)),
    MTE_FMR_fitted = exp(ln_MTE_FMR_fitted),
    ln_MTE_FMR_strict = ln_median_O2_consumption + (M_scaled * -0.25) + (-0.65 * (1 / InverseTemp)),
    MTE_FMR_strict = exp(ln_MTE_FMR_strict)
  )

# Plot the results
plot2 <- ggplot(combined_table, aes(Mass, O2_consumption, color = Species)) +
  geom_point(size = 1) +
  scale_color_manual(values = color_species)+

  theme(legend.position = "none") +
  geom_smooth(aes(Mass, MTE_FMR_fitted, group = Species), color = "black", lwd = 1) +
  geom_smooth(aes(Mass, MTE_FMR_strict, group = Species), color = "red", lwd = 0.5) +
  labs(
    x = "Mass (g)",
    y = "FMR (mgO2 Kg-1 h-1") +
  theme(
    panel.background = element_rect(colour = "black", fill = 'transparent'), # transparent panel bg
    plot.background = element_rect(fill = 'transparent', color = NA), # transparent plot bg
    panel.grid.major = element_blank(), # remove major gridlines
    panel.grid.minor = element_blank(), # remove minor gridlines
    axis.text.y = element_text(size = 18, colour = "black"),
    axis.text.x = element_text(size = 12, colour = "black"), 
    axis.title.x = element_text(size = 20, colour = "black"),
    axis.title.y = element_text(size = 20, colour = "black"),
    plot.margin = unit(c(1, 1, 0, 1), "lines") # Adjust margins for spacing
  )+
  facet_wrap(~ Species, ncol = 2, scales = "free")

print(plot2)

# Plot the results
plot2b <- ggplot(combined_table, aes(Temp, O2_consumption, color = Species)) +
  geom_point(size = 1) +
  scale_color_manual(values = color_species)+
  #geom_smooth(method = "gam", lwd = 1) +
  theme(legend.position = "none") +
  geom_smooth(aes(Temp, MTE_FMR_fitted, group = Species), color = "black", lwd = 1) +
  geom_smooth(aes(Temp, MTE_FMR_strict, group = Species), color = "red", lwd = 0.5) +
  labs(
    x = "Temperature (°C)",
    y = "FMR (mgO2 Kg-1 h-1") +
  theme(
    panel.background = element_rect(colour = "black", fill = 'transparent'), # transparent panel bg
    plot.background = element_rect(fill = 'transparent', color = NA), # transparent plot bg
    panel.grid.major = element_blank(), # remove major gridlines
    panel.grid.minor = element_blank(), # remove minor gridlines
    axis.text.y = element_text(size = 18, colour = "black"),
    axis.text.x = element_text(size = 12, colour = "black"), 
    axis.title.x = element_text(size = 20, colour = "black"),
    axis.title.y = element_text(size = 20, colour = "black"),
    plot.margin = unit(c(1, 1, 0, 1), "lines") # Adjust margins for spacing
  )+
  facet_wrap(~ Species, ncol = 2, scales = "free")

print(plot2b)
# combine plots
Fitted_curves <- plot_grid(plot2, plot2b, ncol = 2, rel_heights = c(1, 1))

# Display the plot
print(Fitted_curves)


C.exoletus_threshold<-20
C.rupestris_threshold<-30
L.bergylta_threshold<-300
L.mixtus_threshold<-150
S.melops_threshold<-75




Wrasse_data <- Wrasse_data %>%
  
  mutate(Species = ifelse(Species == "Ctenolabrus_exoletus" & Mass < C.exoletus_threshold, 
                          "Ctenolabrus_exoletus_S", 
                          ifelse(Species == "Ctenolabrus_exoletus" & Mass >= C.exoletus_threshold, 
                                 "Ctenolabrus_exoletus_L", 
                                 Species)))%>%
  
  mutate(Species = ifelse(Species == "Labrus_bergylta" & Mass < L.bergylta_threshold, 
                          "Labrus_bergylta_S", 
                          ifelse(Species == "Labrus_bergylta" & Mass >= L.bergylta_threshold, 
                                 "Labrus_bergylta_L", 
                                 Species)))%>%
  
  mutate(Species = ifelse(Species == "Labrus_mixtus" & Mass < L.mixtus_threshold, 
                          "Labrus_mixtus_S", 
                          ifelse(Species == "Labrus_mixtus" & Mass >= L.mixtus_threshold, 
                                 "Labrus_mixtus_L", 
                                 Species)))%>%
  
  mutate(Species = ifelse(Species == "Symphodus_melops" & Mass < S.melops_threshold, 
                          "Symphodus_melops_S", 
                          ifelse(Species == "Symphodus_melops" & Mass >= S.melops_threshold, 
                                 "Symphodus_melops_L", 
                                 Species)))%>%
  
  mutate(Species = ifelse(Species == "Ctenolabrus_rupestris" & Mass < C.rupestris_threshold, 
                          "Ctenolabrus_rupestris_S", 
                          ifelse(Species == "Ctenolabrus_rupestris" & Mass >= C.rupestris_threshold, 
                                 "Ctenolabrus_rupestris_L", 
                                 Species)))


  
#z.score mass and inverse temp
Wrasse_data$z_inv_temp <- Z.score(Wrasse_data$InverseTemp, "inv_temp")$z_inv_temp
Wrasse_data$z_ln_mass <- Z.score(Wrasse_data$lnMass,"ln_mass")$z_ln_mass

#collect std dev of ln inv temp and ln mass for back-correction
sd_Inv.Temp<-sd(Wrasse_data$InverseTemp, na.rm=TRUE)
sd_lnMass<-sd(Wrasse_data$lnMass, na.rm=TRUE)




##### INLA glm modelling - all terms as fixed effects, restrict to large
Wrasse_data_large<-Wrasse_data%>%
  filter(Species%in% c("Ctenolabrus_exoletus_L", "Labrus_bergylta_L", "Labrus_mixtus_L", "Symphodus_melops_L", "Ctenolabrus_rupestris_L"))


Wrasse_data_large$Species<-factor(Wrasse_data_large$Species)

median(Wrasse_data_large$Temp, na.rm=TRUE)
median(Wrasse_data_large$Mass, na.rm=TRUE)

Wrasse_data_large$z_inv_temp <- Z.score(Wrasse_data_large$InverseTemp, "inv_temp")$z_inv_temp
Wrasse_data_large$z_ln_mass <- Z.score(Wrasse_data_large$lnMass,"ln_mass")$z_ln_mass

scale(Wrasse_data_large$lnMass, center = TRUE, scale = TRUE)

#collect std dev of ln inv temp and ln mass for back-correction
sd_Inv.Temp<-sd(Wrasse_data_large$InverseTemp, na.rm=TRUE)
sd_lnMass<-sd(Wrasse_data_large$lnMass, na.rm=TRUE)

mean_Inv.Temp<-mean(Wrasse_data_large$InverseTemp, na.rm=TRUE)
mean_lnMass<-mean(Wrasse_data_large$lnMass, na.rm=TRUE)

summary_large<-Wrasse_data_large%>%
  group_by(Species, Location)%>%
  summarize (
    n(),
    median_mass = median(Mass, na.rm=TRUE),
    IQR_mass = IQR(Mass, na.rm=TRUE),
    median_temp = median(Temp, na.rm=TRUE),
    IQR_temp = IQR(Temp, na.rm=TRUE),
    medianFMR = median(O2_consumption, na.rm=TRUE),
    IQRFMR = IQR(O2_consumption, na.rm=TRUE)
    
  )
summary_large

# following Zuur et al Ch 12.  

# I1<- Species and Location as fixed effects on intercept only
I1<-inla(ln_O2_consumption ~ z_ln_mass+ z_inv_temp + Species + Location,

          family = "gaussian",
          #control.family = list(link = "log"),  # Use control.family to specify the link function
          control.predictor = list(compute = TRUE),
          control.compute = list(config=TRUE, dic = TRUE, waic = TRUE),
          data = Wrasse_data_large)

# I2<- Species effect on mass slope and intercept,  Location as fixed effect on intercept only
I2<-inla(ln_O2_consumption ~ (z_ln_mass*Species)+ z_inv_temp  + Location,
         
         family = "gaussian",
         #control.family = list(link = "log"),  # Use control.family to specify the link function
         control.predictor = list(compute = TRUE),
         control.compute = list(config=TRUE, dic = TRUE, waic = TRUE),
         data = Wrasse_data_large)

# I3<- Species effect on  slopes and intercepts, and Location as  fixed effect on intercept only
I3<-inla(ln_O2_consumption ~ (z_ln_mass*Species)+ (z_inv_temp*Species)  + Location,
         
         family = "gaussian",
         #control.family = list(link = "log"),  # Use control.family to specify the link function
         control.predictor = list(compute = TRUE),
         control.compute = list(config=TRUE, dic = TRUE, waic = TRUE),
         data = Wrasse_data_large)

# I4<- Species effect  effects on intercept only, and Location interaction with  mass slope
I4<-inla(ln_O2_consumption ~ (z_ln_mass*Location)+ z_inv_temp  + Species,
         
         family = "gaussian",
         #control.family = list(link = "log"),  # Use control.family to specify the link function
         control.predictor = list(compute = TRUE),
         control.compute = list(config=TRUE, dic = TRUE, waic = TRUE),
         data = Wrasse_data_large)

# I5<- Species effect on intercept only, and Location effect on  with mass and temp slopes
I5<-inla(ln_O2_consumption ~ (z_ln_mass*Location)+ (z_inv_temp*Location)  + Species,
         
         family = "gaussian",
         #control.family = list(link = "log"),  # Use control.family to specify the link function
         control.predictor = list(compute = TRUE),
         control.compute = list(config=TRUE, dic = TRUE, waic = TRUE),
         data = Wrasse_data_large)


# I6<- location effects on temp intercepts and slopes, interaction with location and species
I6<-inla(ln_O2_consumption ~ z_ln_mass+ z_inv_temp*Location *Species*Location,
         
         family = "gaussian",
         #control.family = list(link = "log"),  # Use control.family to specify the link function
         control.predictor = list(compute = TRUE),
         control.compute = list(config=TRUE, dic = TRUE, waic = TRUE),
         data = Wrasse_data_large)

# I7<- Species and location effects on intercepts and slopes with interactions on location
I7<-inla(ln_O2_consumption ~ (z_ln_mass*Species*Location)+ (z_inv_temp*Species*Location)+ (Species*Location),
         
         family = "gaussian",
         #control.family = list(link = "log"),  # Use control.family to specify the link function
         control.predictor = list(compute = TRUE),
         control.compute = list(config=TRUE, dic = TRUE, waic = TRUE),
         data = Wrasse_data_large)


### Model selection

## compare DIC for each of the  models:
Complete_dic<-c(I1$dic$dic,I2$dic$dic,I3$dic$dic, I4$dic$dic, I5$dic$dic, I6$dic$dic, I7$dic$dic) # I1 best model

# Extract the fixed effects summary
model_summary <- summary(I1)

fixed_effects <- as.data.frame(model_summary$fixed)
print(fixed_effects)

# link to sd and mean to recover original units (unscale)
Scaled_fixed<-fixed_effects%>%
  mutate(
    Scaled_fixed <- fixed_effects %>%
    mutate(across(everything(), ~ ifelse(row_number() %in% c(4:7), .+ mean[1], .))),
    mutate(across(everything(), ~ ifelse(row_number() == 2, . / sd_lnMass, .))),
    mutate(across(everything(), ~ ifelse(row_number() == 3, . / sd_Inv.Temp, .))),
    
    
    mutate(across(everything(), ~ ifelse(row_number() == 1, exp(.), .))),
    mutate(across(everything(), ~ ifelse(row_number() %in% c(4:7), exp(.), .))),
  )




# convert coefficients to original units
Common_temp_slope<-fixed_effects$mean[3]/sd_Inv.Temp
Common_temp_slope_lower<-fixed_effects$"0.025quant"[3]/sd_Inv.Temp
Common_temp_slope_upper<-fixed_effects$"0.975quant"[3]/sd_Inv.Temp

Common_mass_slope<-fixed_effects$mean[2]/sd_lnMass
Common_mass_slope_lower<-fixed_effects$"0.025quant"[2]/sd_lnMass
Common_mass_slope_upper<-fixed_effects$"0.975quant"[2]/sd_lnMass


# residual and fitted - model validation
Fit1 <-I1$summary.fitted.values[,"mean"]
E1<- Wrasse_data_large$ln_O2_consumption-Fit1

variance_fitted<-var(Fit1)
variance_resid<-var(E1)

variance_total<-var(Wrasse_data_large$ln_O2_consumption)
variance_resid/variance_fitted


# Generate posterior predictive samples
posterior_samples <- inla.posterior.sample(1000, I1) 

# Extract fitted values
predicted_values <- sapply(posterior_samples, function(x) x$latent)

# Compare with observed data

# Create calibration data
calibration_data <- data.frame(
  Predicted = rowMeans(predicted_values)[1:102],
  Observed =  Wrasse_data_large$ln_O2_consumption,
  residuals <- Wrasse_data_large$ln_O2_consumption - rowMeans(predicted_values)[1:102]
  
)


# Residual plot
residuals_plot<-ggplot(calibration_data, aes(x = Observed, y = residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, color = "red") +
  labs(title = "Residual Plot", x = "Observed", y = "Residuals")


calibration_plot<-ggplot(calibration_data, aes(x = Predicted, y = Observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(title = "Calibration Plot", x = "Predicted values", y = "Observed values")

Model_test_plot <- plot_grid(residuals_plot, calibration_plot, ncol = 1, rel_heights = c(1, 1))
Model_test_plot # outlier low values (below 3.5) poorly fit - may be effects of increasing error at low FMR or genuine low FMR individuals


# plot FMR at median body mass and common temp by species and location




median_z_InvTemp_large <- Wrasse_data_large %>%
 # group_by(Species) %>%
  summarise(
    median_z_Inv_Temp = median(z_inv_temp, na.rm = TRUE)
  )

median_z_lnMass_large <- Wrasse_data_large %>%
  group_by(Species) %>%
  summarise(
    z_ln_mass = median(z_ln_mass, na.rm = TRUE)
  )

# predict (by hand) for common median temperature and median body sizes for each species
Predict_values<-as.data.frame(rbind(median_z_lnMass_large,median_z_lnMass_large))
Predict_values<-Predict_values%>%
  mutate(z_inv_temp = median_z_InvTemp_large$median_z_Inv_Temp,
  Location = as.factor(c(rep("Dorset", 5), rep("Skye", 5))),
    Predicted_lnFMR = 1,
  predicted_lnFMR_lower = 1,
  predicted_lnFMR_upper =1
  
  )
predictedFMR<-rep(NA,10)
predictedFMR[1]<- fixed_effects$mean[1]+(fixed_effects$mean[2]*Predict_values$z_ln_mass[1])+(fixed_effects$mean[3]*Predict_values$z_inv_temp[1])
predictedFMR[2]<- (fixed_effects$mean[1]+fixed_effects$mean[4])+(fixed_effects$mean[2]*Predict_values$z_ln_mass[2])+(fixed_effects$mean[3]*Predict_values$z_inv_temp[1])
predictedFMR[3]<- (fixed_effects$mean[1]+fixed_effects$mean[5])+(fixed_effects$mean[2]*Predict_values$z_ln_mass[3])+(fixed_effects$mean[3]*Predict_values$z_inv_temp[1])
predictedFMR[4]<- (fixed_effects$mean[1]+fixed_effects$mean[6])+(fixed_effects$mean[2]*Predict_values$z_ln_mass[4])+(fixed_effects$mean[3]*Predict_values$z_inv_temp[1])
predictedFMR[5]<- (fixed_effects$mean[1]+fixed_effects$mean[7])+(fixed_effects$mean[2]*Predict_values$z_ln_mass[5])+(fixed_effects$mean[3]*Predict_values$z_inv_temp[1])
predictedFMR[6]<-predictedFMR[1]+fixed_effects$mean[8]
predictedFMR[7]<-predictedFMR[2]+fixed_effects$mean[8]
predictedFMR[8]<-predictedFMR[3]+fixed_effects$mean[8]
predictedFMR[9]<-predictedFMR[4]+fixed_effects$mean[8]
predictedFMR[10]<-predictedFMR[5]+fixed_effects$mean[8]



predictedFMR_lower<-rep(NA,10)
predictedFMR_lower[1]<- fixed_effects$"0.025quant"[1]+(fixed_effects$"0.025quant"[2]*Predict_values$z_ln_mass[1])+(fixed_effects$"0.025quant"[3]*Predict_values$z_inv_temp[1])
predictedFMR_lower[2]<- (fixed_effects$"0.025quant"[1]+fixed_effects$"0.025quant"[4])+(fixed_effects$"0.025quant"[2]*Predict_values$z_ln_mass[2])+(fixed_effects$"0.025quant"[3]*Predict_values$z_inv_temp[1])
predictedFMR_lower[3]<- (fixed_effects$"0.025quant"[1]+fixed_effects$"0.025quant"[5])+(fixed_effects$"0.025quant"[2]*Predict_values$z_ln_mass[3])+(fixed_effects$"0.025quant"[3]*Predict_values$z_inv_temp[1])
predictedFMR_lower[4]<- (fixed_effects$"0.025quant"[1]+fixed_effects$"0.025quant"[6])+(fixed_effects$"0.025quant"[2]*Predict_values$z_ln_mass[4])+(fixed_effects$"0.025quant"[3]*Predict_values$z_inv_temp[1])
predictedFMR_lower[5]<- (fixed_effects$"0.025quant"[1]+fixed_effects$"0.025quant"[7])+(fixed_effects$"0.025quant"[2]*Predict_values$z_ln_mass[5])+(fixed_effects$"0.025quant"[3]*Predict_values$z_inv_temp[1])
predictedFMR_lower[6]<-predictedFMR[1]+fixed_effects$"0.025quant"[8]
predictedFMR_lower[7]<-predictedFMR[2]+fixed_effects$"0.025quant"[8]
predictedFMR_lower[8]<-predictedFMR[3]+fixed_effects$"0.025quant"[8]
predictedFMR_lower[9]<-predictedFMR[4]+fixed_effects$"0.025quant"[8]
predictedFMR_lower[10]<-predictedFMR[5]+fixed_effects$"0.025quant"[8]


predictedFMR_upper<-rep(NA,10)
predictedFMR_upper[1]<- fixed_effects$"0.975quant"[1]+(fixed_effects$"0.975quant"[2]*Predict_values$z_ln_mass[1])+(fixed_effects$"0.975quant"[3]*Predict_values$z_inv_temp[1])
predictedFMR_upper[2]<- (fixed_effects$"0.975quant"[1]+fixed_effects$"0.975quant"[4])+(fixed_effects$"0.975quant"[2]*Predict_values$z_ln_mass[2])+(fixed_effects$"0.975quant"[3]*Predict_values$z_inv_temp[1])
predictedFMR_upper[3]<- (fixed_effects$"0.975quant"[1]+fixed_effects$"0.975quant"[5])+(fixed_effects$"0.975quant"[2]*Predict_values$z_ln_mass[3])+(fixed_effects$"0.975quant"[3]*Predict_values$z_inv_temp[1])
predictedFMR_upper[4]<- (fixed_effects$"0.975quant"[1]+fixed_effects$"0.975quant"[6])+(fixed_effects$"0.975quant"[2]*Predict_values$z_ln_mass[4])+(fixed_effects$"0.975quant"[3]*Predict_values$z_inv_temp[1])
predictedFMR_upper[5]<- (fixed_effects$"0.975quant"[1]+fixed_effects$"0.975quant"[7])+(fixed_effects$"0.975quant"[2]*Predict_values$z_ln_mass[5])+(fixed_effects$"0.975quant"[3]*Predict_values$z_inv_temp[1])
predictedFMR_upper[6]<-predictedFMR[1]+fixed_effects$"0.975quant"[8]
predictedFMR_upper[7]<-predictedFMR[2]+fixed_effects$"0.975quant"[8]
predictedFMR_upper[8]<-predictedFMR[3]+fixed_effects$"0.975quant"[8]
predictedFMR_upper[9]<-predictedFMR[4]+fixed_effects$"0.975quant"[8]
predictedFMR_upper[10]<-predictedFMR[5]+fixed_effects$"0.975quant"[8]

Predict_values<-Predict_values%>%
  mutate(
         Predicted_lnFMR = predictedFMR,
         predicted_lnFMR_lower = predictedFMR_lower,
         predicted_lnFMR_upper = predictedFMR_upper,
         Predicted_FMR = exp(Predicted_lnFMR),
         Predicted_FMR_lower = exp(predicted_lnFMR_lower),
         Predicted_FMR_upper = exp(predicted_lnFMR_upper),
         
  )


# predict fpr common temp (9degrees) and common mass (81 g)
Predict_values_common_mass<-as.data.frame(rbind(median_z_lnMass_large,median_z_lnMass_large))
Predict_values_common_mass<-Predict_values%>%
  mutate(z_inv_temp = median_z_InvTemp_large$median_z_Inv_Temp,
         z_ln_mass = -0.0866,
         Location = as.factor(c(rep("Dorset", 5), rep("Skye", 5))),
         Predicted_lnFMR = 1,
         predicted_lnFMR_lower = 1,
         predicted_lnFMR_upper =1
         
  )
predictedFMR_commonMass<-rep(NA,10)
predictedFMR_commonMass[1]<- fixed_effects$mean[1]+(fixed_effects$mean[2]*Predict_values_common_mass$z_ln_mass[1])+(fixed_effects$mean[3]*Predict_values$z_inv_temp[1])
predictedFMR_commonMass[2]<- (fixed_effects$mean[1]+fixed_effects$mean[4])+(fixed_effects$mean[2]*Predict_values_common_mass$z_ln_mass[2])+(fixed_effects$mean[3]*Predict_values_common_mass$z_inv_temp[1])
predictedFMR_commonMass[3]<- (fixed_effects$mean[1]+fixed_effects$mean[5])+(fixed_effects$mean[2]*Predict_values_common_mass$z_ln_mass[3])+(fixed_effects$mean[3]*Predict_values_common_mass$z_inv_temp[1])
predictedFMR_commonMass[4]<- (fixed_effects$mean[1]+fixed_effects$mean[6])+(fixed_effects$mean[2]*Predict_values_common_mass$z_ln_mass[4])+(fixed_effects$mean[3]*Predict_values_common_mass$z_inv_temp[1])
predictedFMR_commonMass[5]<- (fixed_effects$mean[1]+fixed_effects$mean[7])+(fixed_effects$mean[2]*Predict_values_common_mass$z_ln_mass[5])+(fixed_effects$mean[3]*Predict_values_common_mass$z_inv_temp[1])
predictedFMR_commonMass[6]<-predictedFMR_commonMass[1]+fixed_effects$mean[8]
predictedFMR_commonMass[7]<-predictedFMR_commonMass[2]+fixed_effects$mean[8]
predictedFMR_commonMass[8]<-predictedFMR_commonMass[3]+fixed_effects$mean[8]
predictedFMR_commonMass[9]<-predictedFMR_commonMass[4]+fixed_effects$mean[8]
predictedFMR_commonMass[10]<-predictedFMR_commonMass[5]+fixed_effects$mean[8]



predictedFMR_commonMass_lower<-rep(NA,10)
predictedFMR_commonMass_lower[1]<- fixed_effects$"0.025quant"[1]+(fixed_effects$"0.025quant"[2]*Predict_values_common_mass$z_ln_mass[1])+(fixed_effects$"0.025quant"[3]*Predict_values_common_mass$z_inv_temp[1])
predictedFMR_commonMass_lower[2]<- (fixed_effects$"0.025quant"[1]+fixed_effects$"0.025quant"[4])+(fixed_effects$"0.025quant"[2]*Predict_values_common_mass$z_ln_mass[2])+(fixed_effects$"0.025quant"[3]*Predict_values_common_mass$z_inv_temp[1])
predictedFMR_commonMass_lower[3]<- (fixed_effects$"0.025quant"[1]+fixed_effects$"0.025quant"[5])+(fixed_effects$"0.025quant"[2]*Predict_values_common_mass$z_ln_mass[3])+(fixed_effects$"0.025quant"[3]*Predict_values_common_mass$z_inv_temp[1])
predictedFMR_commonMass_lower[4]<- (fixed_effects$"0.025quant"[1]+fixed_effects$"0.025quant"[6])+(fixed_effects$"0.025quant"[2]*Predict_values_common_mass$z_ln_mass[4])+(fixed_effects$"0.025quant"[3]*Predict_values_common_mass$z_inv_temp[1])
predictedFMR_commonMass_lower[5]<- (fixed_effects$"0.025quant"[1]+fixed_effects$"0.025quant"[7])+(fixed_effects$"0.025quant"[2]*Predict_values_common_mass$z_ln_mass[5])+(fixed_effects$"0.025quant"[3]*Predict_values_common_mass$z_inv_temp[1])
predictedFMR_commonMass_lower[6]<-predictedFMR_commonMass_lower[1]+fixed_effects$"0.025quant"[8]
predictedFMR_commonMass_lower[7]<-predictedFMR_commonMass_lower[2]+fixed_effects$"0.025quant"[8]
predictedFMR_commonMass_lower[8]<-predictedFMR_commonMass_lower[3]+fixed_effects$"0.025quant"[8]
predictedFMR_commonMass_lower[9]<-predictedFMR_commonMass_lower[4]+fixed_effects$"0.025quant"[8]
predictedFMR_commonMass_lower[10]<-predictedFMR_commonMass_lower[5]+fixed_effects$"0.025quant"[8]


predictedFMR_commonMass_upper<-rep(NA,10)
predictedFMR_commonMass_upper[1]<- fixed_effects$"0.975quant"[1]+(fixed_effects$"0.975quant"[2]*Predict_values_common_mass$z_ln_mass[1])+(fixed_effects$"0.975quant"[3]*Predict_values_common_mass$z_inv_temp[1])
predictedFMR_commonMass_upper[2]<- (fixed_effects$"0.975quant"[1]+fixed_effects$"0.975quant"[4])+(fixed_effects$"0.975quant"[2]*Predict_values_common_mass$z_ln_mass[2])+(fixed_effects$"0.975quant"[3]*Predict_values_common_mass$z_inv_temp[1])
predictedFMR_commonMass_upper[3]<- (fixed_effects$"0.975quant"[1]+fixed_effects$"0.975quant"[5])+(fixed_effects$"0.975quant"[2]*Predict_values_common_mass$z_ln_mass[3])+(fixed_effects$"0.975quant"[3]*Predict_values_common_mass$z_inv_temp[1])
predictedFMR_commonMass_upper[4]<- (fixed_effects$"0.975quant"[1]+fixed_effects$"0.975quant"[6])+(fixed_effects$"0.975quant"[2]*Predict_values_common_mass$z_ln_mass[4])+(fixed_effects$"0.975quant"[3]*Predict_values_common_mass$z_inv_temp[1])
predictedFMR_commonMass_upper[5]<- (fixed_effects$"0.975quant"[1]+fixed_effects$"0.975quant"[7])+(fixed_effects$"0.975quant"[2]*Predict_values_common_mass$z_ln_mass[5])+(fixed_effects$"0.975quant"[3]*Predict_values_common_mass$z_inv_temp[1])
predictedFMR_commonMass_upper[6]<-predictedFMR_commonMass_upper[1]+fixed_effects$"0.975quant"[8]
predictedFMR_commonMass_upper[7]<-predictedFMR_commonMass_upper[2]+fixed_effects$"0.975quant"[8]
predictedFMR_commonMass_upper[8]<-predictedFMR_commonMass_upper[3]+fixed_effects$"0.975quant"[8]
predictedFMR_commonMass_upper[9]<-predictedFMR_commonMass_upper[4]+fixed_effects$"0.975quant"[8]
predictedFMR_commonMass_upper[10]<-predictedFMR_commonMass_upper[5]+fixed_effects$"0.975quant"[8]

Predict_values_common_mass<-Predict_values_common_mass%>%
  mutate(
    Predicted_lnFMR = predictedFMR_commonMass,
    predicted_lnFMR_lower = predictedFMR_commonMass_lower,
    predicted_lnFMR_upper = predictedFMR_commonMass_upper,
    Predicted_FMR = exp(Predicted_lnFMR),
    Predicted_FMR_lower = exp(predicted_lnFMR_lower),
    Predicted_FMR_upper = exp(predicted_lnFMR_upper),
    
  )



Location_predict_plot <- ggplot(Predict_values, aes(x = Species, y = log(Predicted_FMR), color = Location)) +
  geom_point(
    aes(group = Location),
    position = position_dodge(width = 0.5), # Matches dodge position of error bars
    size = 6
  ) +
  geom_errorbar(
    aes(group = Location, ymin = log(Predicted_FMR_lower), ymax = log(Predicted_FMR_upper)),
    position = position_dodge(width = 0.5), # Offsets by location within species
    width = 0.3
  ) + 
  ylab(expression(ln~Predicted~FMR~mg~O[2]~kg^{-1}~h^{-1})) +
  ggtitle("Median Species Mass") +
  scale_x_discrete(labels = c(
    "Ctenolabrus_rupestris_L"= "C. ruperstris",
    "Labrus_bergylta_L"= "L. bergylta",
    "Symphodus_melops_L" = "S. melops",
    "Ctenolabrus_exoletus_L"= "C. exoletus",
    "Labrus_mixtus_L"= "L. mixtus"
  ))+
  theme(
    panel.background = element_rect(colour = "black", fill = 'transparent'), # transparent panel bg
    plot.background = element_rect(fill = 'transparent', color = NA), # transparent plot bg
    panel.grid.major = element_blank(), # remove major gridlines
    panel.grid.minor = element_blank(), # remove minor gridlines
    axis.text.y = element_text(size = 16, colour = "black"),
    axis.text.x = element_text(size = 12, colour = "black", angle = 45, hjust = 1), # Rotate x-axis labels
    axis.title.x = element_text(size = 16, colour = "black"),
    axis.title.y = element_text(size = 16, colour = "black"),
    plot.margin = unit(c(1, 1, 0, 1), "lines") # Adjust margins for spacing
  )



Location_predict_common_plot <- ggplot(Predict_values_common_mass, aes(x = Species, y = log(Predicted_FMR), color = Location)) +
  geom_point(
    aes(group = Location),
    position = position_dodge(width = 0.5), # Matches dodge position of error bars
    size = 6
  ) +
  geom_errorbar(
    aes(group = Location, ymin = log(Predicted_FMR_lower), ymax = log(Predicted_FMR_upper)),
    position = position_dodge(width = 0.5), # Offsets by location within species
    width = 0.3
  ) + 
  scale_x_discrete(labels = c(
    "Ctenolabrus_rupestris_L"= "C. ruperstris",
    "Labrus_bergylta_L"= "L. bergylta",
    "Symphodus_melops_L" = "S. melops",
    "Ctenolabrus_exoletus_L"= "C. exoletus",
    "Labrus_mixtus_L"= "L. mixtus"
  ))+
  ylab(expression(ln~Predicted~FMR~mg~O[2]~kg^{-1}~h^{-1})) +
  ggtitle("Common Species Mass") +
  theme(
    panel.background = element_rect(colour = "black", fill = 'transparent'), # transparent panel bg
    plot.background = element_rect(fill = 'transparent', color = NA), # transparent plot bg
    panel.grid.major = element_blank(), # remove major gridlines
    panel.grid.minor = element_blank(), # remove minor gridlines
    axis.text.y = element_text(size = 16, colour = "black"),
    axis.text.x = element_text(size = 12, colour = "black", angle = 45, hjust = 1), # Rotate x-axis labels
    axis.title.x = element_text(size = 16, colour = "black"),
    axis.title.y = element_text(size = 16, colour = "black"),
    plot.margin = unit(c(1, 1, 0, 1), "lines") # Adjust margins for spacing
  )




FMR_large_predict_plot <- plot_grid(Location_predict_plot, Location_predict_common_plot, ncol = 2, rel_heights = c(1, 1))
FMR_large_predict_plot 











#####################################################





# how does experienced temperature vary with body size?

plotMassTemp<-ggplot(Wrasse_data)+
  geom_point(aes(log(Mass,10), Temp, group=Species,color=Species), size=2)+
  geom_smooth(method=lm, formula = y ~ x, se = FALSE, aes(log(Mass,10), Temp, group=Species, color=Species), lwd=0.5)+
  theme(panel.background = element_rect(fill = "transparent"))+
  
  theme(plot.background = element_rect(fill = "transparent", color = NA))+
  theme(axis.text.y = element_text(size = 18)) +
  theme(axis.text.x = element_text(size = 18)) +
  theme(axis.title.x = element_text(size = 22))  +
  theme(axis.title.y = element_text(size = 22)) +
  # theme(panel.background = element_rect(fill = "transparent"))  +
  theme(panel.background = element_rect(colour = "black"))  +
  #  theme(strip.background = element_rect(fill = "white"))   +
  labs(x = "log 10 Mass (g)", y = "Temperature")
plotMassTemp+ facet_wrap( ~ Location, ncol=2)

plotMassTemp_no_log<-ggplot(Wrasse_data)+
  geom_point(aes(Mass, Temp, color=Species), size=2)+
# geom_smooth(method=lm, formula = y ~ x, se = FALSE, aes(log(Weight_g_JR,10), Temp, group=Species, color=Species), lwd=0.5)+
  theme(panel.background = element_rect(fill = "transparent"))+
  
  theme(plot.background = element_rect(fill = "transparent", color = NA))+
  theme(axis.text.y = element_text(size = 18)) +
  theme(axis.text.x = element_text(size = 18)) +
  theme(axis.title.x = element_text(size = 22))  +
  theme(axis.title.y = element_text(size = 22)) +
  # theme(panel.background = element_rect(fill = "transparent"))  +
  theme(panel.background = element_rect(colour = "black"))  +
  #  theme(strip.background = element_rect(fill = "white"))   +
  labs(x = "log 10 Mass (g)", y = "Temperature")
plotMassTemp_no_log+ facet_wrap( ~ Species*Location, ncol=4, scales = "free")



# Mass_specific  O2 consumption vs Temp
MSO2_Temp<-ggplot(Wrasse_data)+
  geom_point(aes(Temp, Mass_specificFMR_100g, group=Species,color=Species), size=2)+
  geom_smooth(method=lm, formula = y ~ x, se = FALSE, aes(Temp, Mass_specificFMR_100g, group=Species, color=Species), lwd=0.5)+
  theme(panel.background = element_rect(fill = "transparent"))+
  
  theme(plot.background = element_rect(fill = "transparent", color = NA))+
  theme(axis.text.y = element_text(size = 18)) +
  theme(axis.text.x = element_text(size = 18)) +
  theme(axis.title.x = element_text(size = 22))  +
  theme(axis.title.y = element_text(size = 22)) +
  # theme(panel.background = element_rect(fill = "transparent"))  +
  theme(panel.background = element_rect(colour = "black"))  +
  #  theme(strip.background = element_rect(fill = "white"))   +
  labs(x = "Temperature", y = "O2 consumption @ 100g (mgO2/Kg/hr)")
MSO2_Temp+ facet_wrap( ~ Location, ncol=2)
MSO2_Temp+ facet_wrap( ~ Species, ncol=3, scales = "free")


##### INLA model here



# extract the slope and use to correct for mass effect to plot mass effect (model simultaneously, but good to visualise this)
# function to apply mass scale and save slopes
fooMass_sens <- function(z) {
  ## coef and se in a data frame
  mr <- data.frame(coef(summary(lm(log(Org_O2,10)~log10Mass, data=z))))
  ## put row names (predictors/indep variables)
  mr$predictor <- rownames(mr)
  mr
}



Whole_animal_Mass_Scale <- ddply(Wrasse_data,"Species", fooMass_sens)
Mass_slopes<-Whole_animal_Mass_Scale %>% 
  filter(predictor=="log10Mass")


cbind(Mass_slopes$Species, Mass_slopes$Estimate)



slope<-as.numeric(Mass_slopeAll$coefficients[2])
mass_scale= slope-1 # slope of mass vs mass-specific O2 consump
Wrasse_data$InverseTemp<-1/((K)*(Wrasse_data$Temp+273))
Wrasse_data$Mass_scaledO2<-(log(Wrasse_data$O2_consumption * Wrasse_data$Weight_g_JR^mass_scale))

#### similar plot with Mass scaling but not scaled to common mass (i.e. Arrhenius plot) -red line - MTE thermal scaling of 0.65eV
Arrhenius_All<-ggplot(Wrasse_data)+
  geom_smooth(method=glm, formula = y ~ x, se = TRUE, aes(InverseTemp, Mass_scaledO2, group=Species, color=Species))+
  geom_point(aes(InverseTemp, Mass_scaledO2, group=Species), colour="dark blue", size=3, pch=21, stroke=0.2)+
  geom_abline(intercept=31, slope=-0.65, col="red", lwd=1, lty=2)+
  xlim(39.4,41)+
  theme(axis.text.y = element_text(size = 16)) +
  theme(axis.text.x = element_text(size = 16)) +
  theme(axis.title.x = element_text(size = 18))  +
  theme(axis.title.y = element_text(size = 18)) +
  theme(panel.background = element_rect(fill = "white"))  +
  theme(panel.background = element_rect(colour = "black"))  +
  theme(strip.background = element_rect(fill = "white"))   +
  labs(x = "Inverse temperature  (1/kT)", y = "Ln (O2 consumption*M^-0.1)")

Arrhenius_All
Arrhenius_All + facet_wrap( ~ Species, ncol=3)


## function to extract the linear regression coefficients  across species
fooTemp_sens <- function(z) {
  ## coef and se in a data frame
  mr <- data.frame(coef(summary(lm(Mass_scaledO2~InverseTemp, data=z))))
  ## put row names (predictors/indep variables)
  mr$predictor <- rownames(mr)
  mr
}

Whole_animal_Temp_Scale <- ddply(Wrasse_data, "Species",fooTemp_sens)
slopesTemp<-Whole_animal_Temp_Scale %>% 
  filter(predictor=="InverseTemp")

cbind(slopesTemp$Species, slopesTemp$Estimate)






# population slopes label plot

#function to generate predict values to allow labels at end of lines
predVals <- function(data, x, y) {
  model <- lm(y ~ x, data)
  predVal <- predict(model, newdata = data)
  return(predVal)
}


Pop_plot<-Wrasse_data %>%
  group_by(Species) %>%
  group_modify(~ mutate(.x, PredO2 = predVals(.x, Temp, Mass_specificFMR_10g))) %>%
  mutate(label = if_else(Temp == max(Temp, na.rm = TRUE), as.character(Species), NA_character_)) %>%
  ggplot(aes(x = Temp, y = Mass_specificFMR_10g, group = Species, color=Species)) + 
  geom_point(cex = 1.5) +
  stat_smooth(
    method = lm,
    se = TRUE,
    aes(Temp, PredO2, group = Species, color = Species),
    lwd = 1.5
  ) +
  guides(color = "none")+
  geom_label_repel(aes(label = label),
                   nudge_x = 2,
                   nudge_y = 1.5,
                   cex=3.5,
                   na.rm = TRUE,
                   label.size = 0.1, 
                   max.overlaps = 25) +  # Increase the max.overlaps value
  #scale_colour_npg()+
  geom_smooth(
    data = Wrasse_data,
    # filter(PopID %in% north),
    # filter(Data_set=="SA_Deep"),
    aes(x = Temp, y = Mass_specificFMR_10g),
    method = "lm",
    color = "black",
    lty=2,
    se = TRUE,
    inherit.aes = FALSE
  ) +
  theme(plot.margin = unit(c(1, 3, 1, 1), "lines"))+
  theme(axis.text.y = element_text(size = 16)) +
  theme(axis.text.x = element_text(size = 16)) +
  theme(axis.title.x = element_text(size = 18))  +
  theme(axis.title.y = element_text(size = 18)) +
  theme(panel.background = element_rect(fill = "white"))  +
  theme(panel.background = element_rect(colour = "black"))  +
  theme(strip.background = element_rect(fill = "white"))   +
  labs(x = "Temperature", y = "FMR[10g] (mgO2/Kg/hr)")
Pop_plot

library(lme4)
library(merTools)
## glmm non Bayesian (split by mass class and temp)
LME_1<-lmer(formula = ln_O2_consumption ~z_ln_Mass + z_inv_temp+(1| Month) + (1|Year),  data   = Wrasse_data)
summary(LME_1)
r.squaredGLMM(LME_1)


# Extract variance components
variance_components <- as.data.frame(VarCorr(LME_1))
random_effect_variance <- variance_components$vcov[variance_components$grp == "Year"]
residual_variance <- attr(VarCorr(LME_1), "sc")^2
total_variance <- random_effect_variance + residual_variance
proportion_random_variance <- random_effect_variance / total_variance


# Calculate R-squared values
r2 <- r.squaredGLMM(LME_1)
marginal_r2 <- r2[1]
conditional_r2 <- r2[2]

# Output the results
cat("Proportion of variance explained by random effects:", proportion_random_variance, "\n")
cat("Marginal R-squared (variance explained by fixed effects):", marginal_r2, "\n")
cat("Conditional R-squared (variance explained by fixed + random effects):", conditional_r2, "\n")


proportion_fixed<-marginal_r2/conditional_r2 # prportion of explained variance attributable to fixed effects (relative effect size)

# Extract fixed effects 
fix_eff <- data.frame(fixef(LME_1))

# Extract random effects 
rand_eff <- data.frame(ranef(LME_1))


