setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")
# packages ####

library(tidyverse)
library(lme4)

# load data ####

data <- read.csv("./01_RAW_DATA/PYT52_TRAIT_DATA.csv")

data = data %>%
  mutate(
    plot_number = paste0("P",plot_number)
  ) %>%
  gather(
    key = "Trait",
    value = "y",
    -accession_name,-year_harvest,-plot_number
  )


# Manual BLUEs ####

dir.create("./XX_FIGURES/01_BLUEs/", recursive = TRUE)
dir.create("./02_BLUEs_CLUSTERS/", recursive = TRUE)

BLUEs = list()

for(i in unique(data$Trait)) {
  
  print(i)
  
  input = data[
    data$Trait == i,
  ]
  
  
  fixed_model = lmer(
    data = input,
    REML = TRUE,
    formula = y ~ 0 + (1|year_harvest) + accession_name
  )
  
  
  png(
    filename = paste0("./XX_FIGURES/01_BLUEs/",i,".png"),
    units = "cm",
    width = 18,
    height = 9,
    res = 300
  )
  
  Residuals = resid(fixed_model)
  Fitted = fitted(fixed_model)
  
  par(mfrow=c(1,4))
  
  hist(Residuals, main = "")
  
  qqnorm(Residuals, main ="")
  qqline(Residuals)
  
  plot(Residuals~Fitted)
  abline(a=mean(Residuals),b=0)
  
  plot(Residuals)
  
  dev.off()
  
  BLUEs[[i]] = data.frame(
    accession_name = gsub(names(fixef(fixed_model)), pattern = "accession_name", replacement = ""),
    Trait = i,
    BLUE = summary(fixed_model)$coefficients[,1],
    Std_error = summary(fixed_model)$coefficients[,2],
    row.names = NULL
  )
  
  
  rm(i,input,Residuals,Fitted,fixed_model)
  
} # i



BLUEs = bind_rows(BLUEs)

write.csv(
  BLUEs,
  file = "./02_BLUEs_CLUSTERS/BLUEs.csv",
  row.names = FALSE
)

Clusters = spread(
  BLUEs %>% select(-Std_error),
  key = "Trait",
  value = "BLUE"
)

# split white and yellow genotypes ####


input = na.omit(Clusters$YELLOWNESS)


WSS = list()


for(i in seq(1:10)) {
  
  
  WSS[[i]] = data.frame(
    k = i,
    WSS = sum(kmeans(input,i)$withinss)
  )
  
}

WSS = bind_rows(WSS) %>%
  mutate(
    k = factor(k,levels = as.character(sort(k)))
  )

plot(
  seq_along(WSS$k), WSS$WSS, type="b",
  xlab="Number of groups",
  ylab="Sum of squares within a group"
)

# use k = 2 as expected for white yellow
set.seed(1)
kmeans_res = kmeans(input,2)

centers = data.frame(
  cluster = c(1,2),
  center = kmeans_res$centers
)
centers = centers[order(centers$center),] %>%
  mutate(
    white_yellow = c("White","Yellow")
  )

Clusters = data.frame(
  YELLOWNESS = input,
  cluster = kmeans_res$cluster
) %>%
  left_join(
    centers[,c("cluster","white_yellow")]
  ) %>%
  left_join(
    Clusters
  ) %>%
  select(
    names(Clusters),white_yellow
  )


# DM k-means clustering ####



input = na.omit(Clusters$DM)


WSS = list()


for(i in seq(1:10)) {
  
  
  WSS[[i]] = data.frame(
    k = i,
    WSS = sum(kmeans(input,i)$withinss)
  )
  
}

WSS = bind_rows(WSS) %>%
  mutate(
    k = factor(k,levels = as.character(sort(k)))
  )

plot(
  seq_along(WSS$k), WSS$WSS, type="b",
  xlab="Number of groups",
  ylab="Sum of squares within a group"
)

set.seed(1)
kmeans_res = kmeans(
  input,3
)

centers = data.frame(
  cluster = c(1,2,3),
  center = kmeans_res$centers
)
centers = centers[order(centers$center),] %>%
  mutate(
    DM_group = c("Low","Medium", "High")
  )

Clusters = data.frame(
  DM = input,
  cluster = kmeans_res$cluster
) %>%
  left_join(
    centers[,c("cluster","DM_group")]
  ) %>%
  right_join(Clusters) %>%
  select(
    names(Clusters),DM_group
  )

# RFW k-means clustering ####



input = na.omit(Clusters$RTWT)


WSS = list()


for(i in seq(1:10)) {
  
  
  WSS[[i]] = data.frame(
    k = i,
    WSS = sum(kmeans(input,i)$withinss)
  )
  
}

WSS = bind_rows(WSS) %>%
  mutate(
    k = factor(k,levels = as.character(sort(k)))
  )

plot(
  seq_along(WSS$k), WSS$WSS, type="b",
  xlab="Number of groups",
  ylab="Sum of squares within a group"
)

rm(WSS)

set.seed(1)
kmeans_res = kmeans(input,3)


centers = data.frame(
  cluster = c(1,2,3),
  center = kmeans_res$centers
)
centers = centers[order(centers$center),] %>%
  mutate(
    RTWT_group = c("Low","Medium", "High")
  )

Clusters = data.frame(
  RTWT = input,
  cluster = kmeans_res$cluster
) %>%
  left_join(
    centers[,c("cluster","RTWT_group")]
  ) %>%
  right_join(Clusters) %>%
  select(
    names(Clusters),RTWT_group
  )

write.csv(
  Clusters,
  file = "./02_BLUEs_CLUSTERS/CLUSTERS_ALL_GTs.csv",
  row.names = FALSE
)

# selected Genotypes ####

genotypes = c(
  "IITA-TMS-ZAR010116",
  "IITA-TMS-IBA990304",
  "IITA-TMS-IBA071313",
  "IITA-TMS-IBA982101",
  "TMS13F1343P0002",
  "IITA-TMS-IBA980581",
  "TMS13F1088P0007",
  "IITA-TMS-IBA090576"
)


genotypes = Clusters[Clusters$accession_name %in% genotypes,]

genotypes = genotypes[order(genotypes$DM),]

write.csv(
  genotypes,
  file = "./02_BLUEs_CLUSTERS/CLUSTERS_SELECTED_GENOTYPES.csv",
  row.names = FALSE
)


# plots ####

dir.create(
  "./XX_FIGURES/02_TRAIT_CLUSTERS/",
  recursive = TRUE
)

plot_df = BLUEs %>%
  inner_join(
    Clusters %>%
      select(accession_name, white_yellow, RTWT_group, DM_group)
  ) %>%
  mutate(
    white_yellow = factor(white_yellow,levels = c("Yellow", "White")),
    RTWT_group = factor(RTWT_group, levels = c("Low", "Medium", "High"))
  )


ggplot(
  plot_df[plot_df$Trait == "RTWT",],
  aes(
    x = BLUE,
    y = white_yellow
  )
) +
  geom_point(
    shape=21,
    aes(fill = RTWT_group),
    size = 1,
    stroke = 0.25
  ) +
  labs(
    x = "Root fresh weight BLUE [kg]",
    y = "White/Yellow",
    fill = "Group"
  ) +
  theme_bw() +
  scale_fill_manual(values = c("blue","white","red")) +
  theme(
    axis.title = element_text(size = 6, face = "bold"),
    axis.text = element_text(size = 5),
    legend.position = "right",
    legend.justification = "top",
    legend.title = element_text(size = 5, face = "bold"),
    legend.text = element_text(size = 4),
    legend.background = element_rect(fill = "transparent", color = "black", size = 0.25),
    legend.key.size = unit(.25,"cm"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent", color = NA),
    panel.border = element_blank(),
    axis.line = element_line(color = 'black', size = 0.25),
    strip.text = element_blank(),
    strip.background = element_rect(fill = NA, color = NA)
  )
ggsave(
  path = "./XX_FIGURES/02_TRAIT_CLUSTERS/",
  filename = "01_RTWT.png",
  units = "cm",
  width = 10,
  height = 7,
  dpi = 900
)

plot_df = BLUEs[BLUEs$Trait %in% c("RTWT", "DM", "SHTWT","HI"),] %>%
  inner_join(
    genotypes %>%
      select(accession_name, white_yellow, RTWT_group, DM_group)
  ) %>%
  mutate(
    accession_name = factor(accession_name, levels = unique(genotypes$accession_name)),
    Trait = factor(Trait, levels = c("DM","RTWT","SHTWT","HI")),
    white_yellow = factor(white_yellow,levels = c("Yellow", "White")),
    DM_group = factor(DM_group, levels = c("Low", "Medium", "High"))
  )

levels(plot_df$Trait) = c("Dry matter content [%]", "Root fresh weight [kg]", "Shoot fresh weight [kg]", "Harvest index")

ggplot(
  plot_df,
  aes(
    x = accession_name,
    y = BLUE,
    fill = DM_group
  )
) +
  geom_errorbar(
    aes(
      ymin = BLUE - Std_error,
      ymax = BLUE + Std_error
    ),
    size = 0.2,
    width = 0.1
  ) +
  geom_point(
    shape=21,
    aes(fill = DM_group),
    size = 1.5,
    stroke = 0.25
  )  +
  labs(
    x = "Genotype",
    y = paste0("BLUE \u00B1SE"),
    fill = "DM% Group"
  ) +
  scale_y_continuous(
    #limits = c(0,NA)
  ) + 
  facet_wrap(~Trait, scales = "free_y") +
  theme_bw() +
  scale_fill_manual(values = c("blue","white","red")) +
  theme(
    axis.title = element_text(size = 6, face = "bold"),
    axis.text.x = element_text(size = 5, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 5),
    legend.position = "right",
    legend.justification = "top",
    legend.title = element_text(size = 5, face = "bold"),
    legend.text = element_text(size = 4),
    legend.background = element_rect(fill = "transparent", color = "black", size = 0.25),
    legend.key.size = unit(.25,"cm"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent", color = NA),
    panel.border = element_blank(),
    axis.line = element_line(color = 'black', size = 0.25),
    strip.text = element_text(size = 6, face = "bold"),
    strip.background = element_rect(fill = NA, color = NA)
  )
ggsave(
  path = "./XX_FIGURES/02_TRAIT_CLUSTERS/",
  filename = "02_Selected_GTs_all.png",
  units = "cm",
  width = 12,
  height = 10,
  dpi = 900
)

# selected traits

plot_df = BLUEs[BLUEs$Trait %in% c("RTWT", "DM"),] %>%
  inner_join(
    genotypes %>%
      select(accession_name, white_yellow, RTWT_group, DM_group)
  ) %>%
  mutate(
    accession_name = factor(accession_name, levels = unique(genotypes$accession_name)),
    Trait = factor(Trait, levels = c("DM","RTWT")),
    white_yellow = factor(white_yellow,levels = c("Yellow", "White")),
    DM_group = factor(DM_group, levels = c("Low", "Medium", "High"))
  )

levels(plot_df$Trait) = c("Dry matter content [%]", "Root fresh weight [kg]")

ggplot(
  plot_df,
  aes(
    x = accession_name,
    y = BLUE,
    fill = DM_group
  )
) +
  geom_errorbar(
    aes(
      ymin = BLUE - Std_error,
      ymax = BLUE + Std_error
    ),
    size = 0.2,
    width = 0.1
  ) +
  geom_point(
    shape=21,
    aes(fill = DM_group),
    size = 1.5,
    stroke = 0.25
  )  +
  labs(
    x = "Genotype",
    y = paste0("BLUE \u00B1SE"),
    fill = "DM% Group"
  ) +
  scale_y_continuous(
    #limits = c(0,NA)
  ) + 
  facet_wrap(~Trait, scales = "free_y") +
  theme_bw() +
  scale_fill_manual(values = c("blue","white","red")) +
  theme(
    axis.title = element_text(size = 6, face = "bold"),
    axis.text.x = element_text(size = 5, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 5),
    legend.position = "right",
    legend.justification = "top",
    legend.title = element_text(size = 5, face = "bold"),
    legend.text = element_text(size = 4),
    legend.background = element_rect(fill = "transparent", color = "black", size = 0.25),
    legend.key.size = unit(.25,"cm"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent", color = NA),
    panel.border = element_blank(),
    axis.line = element_line(color = 'black', size = 0.25),
    strip.text = element_text(size = 6, face = "bold"),
    strip.background = element_rect(fill = NA, color = NA)
  )
ggsave(
  path = "./XX_FIGURES/02_TRAIT_CLUSTERS/",
  filename = "02_Selected_GTs_DM_RTWT.png",
  units = "cm",
  width = 14,
  height = 8,
  dpi = 900
)




# white vs yellow



plot_df = BLUEs %>%
  inner_join(
    Clusters %>%
      select(accession_name, white_yellow)
  ) %>%
  mutate(
    white_yellow = factor(white_yellow,levels = c("White", "Yellow")),
    accession_name = factor(accession_name, levels = unique(Clusters[order(Clusters$YELLOWNESS),]$accession_name))
  )


ggplot(
  plot_df[plot_df$Trait == "YELLOWNESS",],
  aes(
    x = accession_name,
    y = BLUE
  )
) +
  geom_errorbar(
    aes(
      ymin = BLUE - Std_error,
      ymax = BLUE + Std_error
    ),
    size = 0.2,
    width = 0.2
  ) +
  geom_point(
    shape=21,
    aes(fill = white_yellow),
    size = 1,
    stroke = 0.25
  ) +
  labs(
    y = "Yellownes BLUE \u00B1SE",
    x = "Genotype",
    fill = "Group"
  ) +
  theme_bw() +
  scale_fill_manual(values = c("white", "yellow")) +
  theme(
    axis.title = element_text(size = 6, face = "bold"),
    axis.text.x = element_text(size = 5, angle = 90, hjust = 1, vjust = .5),
    axis.text.y = element_text(size = 5),
    legend.position = "right",
    legend.justification = "top",
    legend.title = element_text(size = 5, face = "bold"),
    legend.text = element_text(size = 4),
    legend.background = element_rect(fill = "transparent", color = "black", size = 0.25),
    legend.key.size = unit(.25,"cm"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent", color = NA),
    panel.border = element_blank(),
    axis.line = element_line(color = 'black', size = 0.25),
    strip.text = element_blank(),
    strip.background = element_rect(fill = NA, color = NA)  )
ggsave(
  path = "./XX_FIGURES/02_TRAIT_CLUSTERS",
  filename = "00_YELLOWNESS.png",
  units = "cm",
  width = 10,
  height = 7,
  dpi = 900
)

