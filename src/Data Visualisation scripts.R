------------------------
#M206 PCoA visualisation
------------------------
  
#Plotting MDS M206 data based on compartment, season, and plant age
#Line to save code:
ggsave("graph_title", plot = , width = 8, height = 6, dpi = 300)

#Compartment
gM206Comp<-ggplot(mds_merged, aes(x = V1, y = V2, color = Compartment)) +
  geom_point(size = 1.8, alpha = 0.8) +  # Slightly larger points, some transparency
  scale_color_viridis_d(option = "D") +  # Viridis discrete color scale, colorblind-friendly
  labs(
    x = "PC1",
    y = "PC2",
    title = "PCoA of M206 Samples by Compartment",
    color = "Compartment"
  ) +
  theme_minimal(base_size = 13) +  # Improves text readability
  theme(
    plot.title = element_text(hjust = 0, vjust = 1.8,),
    legend.position = "right"
  )
ggsave("PCoA_M206_by_Compartment.png", plot = gM206Comp, width = 8, height = 6, dpi = 300)

#Season
mds_merged$Season <- factor(mds_merged$Season, levels = c("2014", "2015", "2016"))
gM206Seas <- ggplot(mds_merged, aes(x = V2, y = V3, color = Season)) +
  geom_point(size = 1.8, alpha = 0.5) +  # Slightly larger points, some transparency
  scale_color_manual(values = season_colors) +
  scale_fill_manual(values = season_colors) +
  labs(
    x = "PC1",
    y = "PC2",
    title = "PCoA of M206 Samples by Season",
    color = "Season"
  ) +
  theme_minimal(base_size = 13) +  # Improves text readability
  theme(
    plot.title = element_text(hjust = 0, vjust = 1.8,),
    legend.position = "bottom",                   # Move legend to bottom
    legend.box = "horizontal",                    # Layout legend items horizontally
    legend.direction = "horizontal",              # Align entries in a row
    legend.title = element_text(size = 11),       # Optional: adjust legend text
    legend.text = element_text(size = 10)
  )
ggsave("PCoA_M206_by_Season.png", plot = gM206Seas, width = 8, height = 6, dpi = 300)


#Age
gM206Age <- ggplot(mds_merged, aes(x = V1, y = V2, color = Age)) +
  geom_point(size = 1.5) +
  scale_color_gradient(low = "blue", high = "orange") +
  labs(x = "PC1", y = "PC2", color = "Plant Age (days)", 
       title = "PCoA of M206 Samples by Plant Age") +
  theme_classic()+
  theme(
    plot.title = element_text(
      hjust = 0.1,  
      vjust = 3,      
      size = 14,        
    )
  )




------------------------
#2016 PCoA visualisation
------------------------

#Age
g2016Age <- ggplot(mds_merged2, aes(x = V1, y = V3, color = Age)) +
  geom_point(size = 1.5) +
  scale_color_gradient(low = "blue", high = "orange") +
  labs(x = "PC1", y = "PC3", color = "Plant Age (days)", 
       title = "PCoA of 2016 Samples by Plant Age") +
  theme_classic()+
  theme(
    plot.title = element_text(
      hjust = 0.1,  
      vjust = 3,      
      size = 14,        
    )
  )

#Compartment
g2016Comp<-ggplot(mds_merged2, aes(x = V1, y = V2, color = Compartment)) +
  geom_point(size = 1.8, alpha = 0.8) +  # Slightly larger points, some transparency
  scale_color_viridis_d(option = "D") +  # Viridis discrete color scale, colorblind-friendly
  labs(
    x = "PC1",
    y = "PC2",
    title = "PCoA of 2016 Samples by Compartment",
    color = "Compartment"
  ) +
  theme_minimal(base_size = 13) +  # Improves text readability
  theme(
    plot.title = element_text(hjust = 0, vjust = 1.8,),
    legend.position = "right"
  )
ggsave("PCoA_2016_by_Compartment.png", plot = g2016Comp, width = 8, height = 6, dpi = 300)


#Cultivar
g2016Cult <- ggplot(mds_merged2, aes(x = V1, y = V2, color = Cultivar)) +
  geom_point(size = 1) +
  scale_color_manual(values = cultivar_colors) +
  scale_fill_manual(values = cultivar_colors) +
  labs(x = "PC1", y = "PC2", title = "PCoA of 2016 Samples by Cultivar") +
  theme_minimal(base_size = 13) +  # Improves text readability
  theme(
    plot.title = element_text(hjust = 0, vjust = 1.8,),
    legend.position = "bottom",                   # Move legend to bottom
    legend.box = "horizontal",                    # Layout legend items horizontally
    legend.direction = "horizontal",              # Align entries in a row
    legend.title = element_text(size = 11),       # Optional: adjust legend text
    legend.text = element_text(size = 10)
  )
ggsave("PCoA_2016_by_Cultivar.png", plot = g2016Cult, width = 8, height = 6, dpi = 300)



------------------------------------------------
#Fishers Alpha diversity calculation visualtions
------------------------------------------------
  
-----------------------------
#M206 diversity visualisation
-----------------------------
  
gAFComAge <- ggplot(alpha_merged, aes(x = Age, y = FisherAlpha, color = Compartment)) +
  geom_point(alpha = 0.6) +
  scale_color_viridis_d(option = "D") +  # Viridis discrete color scale, colorblind-friendly
  geom_smooth(method = "lm", se = TRUE, color = "black") +  # keeps a neutral trend line
  labs(
    title = "Alpha Diversity Across Age by Compartment for M206",
    x = "Plant Age (Days)",
    y = "Fisher’s Alpha Diversity",
    color = "Compartment"
  ) +
  theme_minimal(base_size = 13)
ggsave("AFdiversity_M206_by_Compartment.png", plot = gAFComAge, width = 8, height = 6, dpi = 300)



gAFComAgeVio <- ggplot(alpha_merged, aes(x = AgeBin, y = FisherAlpha, fill = Compartment)) +
  geom_violin(scale = "width", trim = FALSE, alpha = 0.7) +
  geom_jitter(aes(color = Compartment), width = 0.15, size = 0.8, alpha = 0.5) +
  scale_fill_viridis_d(option = "D") +   # Apply Viridis to fill
  scale_color_viridis_d(option = "D") +  # Apply same Viridis to color
  labs(
    title = "AlphaF Diversity by Age Bin and Compartment for M206",
    x = "Age Bin (Days)",
    y = "Fisher’s Alpha Diversity",
    fill = "Compartment",
    color = "Compartment"
  ) +
  theme_minimal(base_size = 13)
ggsave("AFdiversity_M206_by_Compartment2.png", plot = gAFComAgeVio, width = 8, height = 6, dpi = 300)


# Summarise alpha diversity by Age and Season
alpha_summary <- alpha_merged %>%
  group_by(Age, Season) %>%
  summarise(
    mean_alpha = mean(FisherAlpha, na.rm = TRUE),
    sd = sd(FisherAlpha, na.rm = TRUE),
    .groups = 'drop'
  )

# Ensure Season is treated as a categorical factor
alpha_summary$Season <- factor(alpha_summary$Season, levels = c("2014", "2015", "2016"))

# Define better color palette
season_colors <- c("2014" = "#E69F00", "2015" = "#56B4E9", "2016" = "#009E73")

# Plot
gAFSeasAge <- ggplot(alpha_summary, aes(x = Age, y = mean_alpha, color = Season, group = Season)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = season_colors) +
  scale_fill_manual(values = season_colors) +
  labs(
    title = "Mean Fisher’s Alpha Diversity Over Age by Season",
    x = "Plant Age (days)",
    y = "Mean Fisher’s Alpha Diversity"
  ) +
  theme_minimal(base_size = 13)
ggsave("AFdiversity_M206_by_Season.png", plot = gAFSeasAge, width = 8, height = 6, dpi = 300)





-----------------------------
#2016 diversity visualisation
-----------------------------

  
  gAFComAge2 <- ggplot(alpha_merged2, aes(x = Age, y = FisherAlpha, color = Compartment)) +
  geom_point(alpha = 0.6) +
  scale_color_viridis_d(option = "D") +  # Viridis discrete color scale, colorblind-friendly
  geom_smooth(method = "lm", se = TRUE, color = "black") +  # keeps a neutral trend line
  labs(
    title = "Alpha Diversity Across Age by Compartment in 2016",
    x = "Plant Age (Days)",
    y = "Fisher’s Alpha Diversity",
    color = "Compartment"
  ) +
  theme_minimal(base_size = 13)
ggsave("AFdiversity_2016_by_Compartment.png", plot = gAFComAge2, width = 8, height = 6, dpi = 300)




gAFComAgeVio2 <- ggplot(alpha_merged2, aes(x = AgeBin, y = FisherAlpha, fill = Compartment)) +
  geom_violin(scale = "width", trim = FALSE, alpha = 0.7) +
  geom_jitter(aes(color = Compartment), width = 0.15, size = 0.8, alpha = 0.5) +
  scale_fill_viridis_d(option = "D") +   # Apply Viridis to fill
  scale_color_viridis_d(option = "D") +  # Apply same Viridis to color
  labs(
    title = "AlphaF Diversity by Age Bin and Compartment in 2016",
    x = "Age Bin (Days)",
    y = "Fisher’s Alpha Diversity",
    fill = "Compartment",
    color = "Compartment"
  ) +
  theme_minimal(base_size = 13)
ggsave("AFdiversity_2016_by_Compartment2.png", plot = gAFComAgeVio2, width = 8, height = 6, dpi = 300)




# Summarise mean ± SD for each age × cultivar
alpha_summary2 <- alpha_merged2 %>%
  group_by(Age, Cultivar) %>%
  summarise(
    mean_alpha = mean(FisherAlpha, na.rm = TRUE),
    sd = sd(FisherAlpha, na.rm = TRUE),
    .groups = 'drop'
  )

# Assign colors for each cultivar
cultivar_colors <- RColorBrewer::brewer.pal(length(unique(alpha_summary2$Cultivar)), "Set1")
names(cultivar_colors) <- levels(alpha_summary2$Cultivar)

# Plot
gAFCultAge <- ggplot(alpha_summary2, aes(x = Age, y = mean_alpha, color = Cultivar, group = Cultivar)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = cultivar_colors) +
  scale_fill_manual(values = cultivar_colors) +
  labs(
    title = "Mean AlphaF Diversity Over Age by Cultivar in 2016",
    x = "Plant Age (days)",
    y = "Mean Fisher’s Alpha Diversity"
  ) +
  theme_minimal(base_size = 13)
ggsave("AFdiversity_2016_by_Cultivar.png", plot = gAFCultAge, width = 8, height = 6, dpi = 300)


