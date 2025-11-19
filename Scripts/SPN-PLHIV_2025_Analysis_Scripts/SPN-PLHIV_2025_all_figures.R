#################################### LOAD REQUIRED PACKAGES ####################
library(tidyverse)
library(ggpubr)
#################################### FIGURE 1 ##################################
###################################################################### Figure 1a
Fig1a <- magick::image_read("HIV-PAPER/Figures/Figure 1/Fig1a.png")
Fig1a <- ggdraw()+ 
  draw_image(Fig1a, scale = 1.8)+
  labs(title = "a.")+
  theme(
    plot.title = element_text(hjust = 0.01,vjust = .01, size = 50, face = "bold"),  # Customize title
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")  # Optional: adjust margins
  )
Fig1a
###################################################################### Figure 1b 
Figure_1_data <- read_csv("HIV-PAPER/Data/SPN-PLHIV_2025_data/Figure_1_data.csv")
Fig1b <- Figure_1_data %>%
  filter(
    `HIV Status`=="HIV-"
  ) %>%
  pivot_longer(cols = c("NER","TER","MER"),
               names_to = "Immune cells",
               values_to = "Immune cell to epithelial cell ratio") %>%
  ggplot(aes(x=factor(`Immune cells`,
                      levels=c("NER","MER","TER")),
             y=log10(`Immune cell to epithelial cell ratio`),
             shape=factor(`Immune cells`,
                          levels=c("NER","MER","TER")),
             color=factor(`Immune cells`,
                          levels=c("NER","MER","TER")),
             fill=factor(`Immune cells`,
                         levels=c("NER","MER","TER"))))+
  geom_boxplot(width=0.5,notch = F, outlier.shape = NA,fill="white")+
  geom_jitter(position = position_jitter(width = 0.2),size=7) +
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=10,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 'wilcox.test',
           label = "{paste0('p = ', signif(p, digits = 4))}", 
           #label = "{ifelse(p < 0.0001, 'p < 0.0001', sprintf('p = %.4f', p))}",
           label.size = 10,
           tip.length = 0.01,
           hide.ns = F,
           p.adjust.method = "holm",
           vjust = -0.2)+
  scale_color_manual(values = c('#000','#000','#000'))+
  scale_fill_manual(values = c('#A9A9A9','#A9A9A9','#A9A9A9'))+
  scale_shape_manual(values = c(21, 24, 22))+
  labs(x='',
       y=expression("Immune cell to epithelial cell ratio"~"("~"Log"[10]~")"),
       fill='Immune cells',
       title = "b.")+
  scale_y_continuous(limits = c(-3,1.3))+
  scale_x_discrete(labels=c("MER"="CD14+\nMonocytes",
                            "TER"="CD3+\n T cells",
                            "NER"="CD66b+\nNeutrophils"))+
  theme_classic()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        plot.title = element_text(size = 50, face = "bold"),
        axis.ticks.length = unit(0.5,"cm"),
        legend.text = element_text(size=25),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 30))

Fig1b

###################################################################### Figure 1c
Fig1c <- Figure_1_data %>%
  ggplot(aes(x=factor(`HIV Status`,levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y')),
             y=log10(as.numeric(`NER`)),
             color=factor(`HIV Status`,
                          levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y'))))+
  geom_boxplot(width=0.5,notch = F, outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.2),size=7) +
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=7,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 'wilcox.test',
           label = "{paste0('p = ', signif(p, digits = 4))}", 
           #label = "{ifelse(p < 0.0001, 'p < 0.0001', sprintf('p = %.4f', p))}",
           label.size = 10,
           tip.length = 0.01,
           hide.ns = F,
           p.adjust.method = "holm",
           vjust = -0.2)+
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  labs(x='',
       y=expression("Neutrophil to epithelial cell ratio"~"("~"Log"[10]~")"),
       title = "c.")+
  scale_x_discrete(labels=c("HIV-"="HIV-",
                            "PLHIV ART <3m"="PLHIV\nART<3m",
                            "PLHIV ART >1y"="PLHIV\nART>1yr"))+
  theme_classic()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        plot.title = element_text(size = 50, face = "bold"),
        axis.ticks.length = unit(0.5,"cm"),
        legend.text = element_text(size=25),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 30))

Fig1c

###################################################################### Figure 1d
Fig1d <- Figure_1_data %>%
  ggplot(aes(x=factor(`HIV Status`,levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y')),
             y=log10(as.numeric(`MER`)),
             color=factor(`HIV Status`,
                          levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y'))))+
  geom_boxplot(width=0.5,notch = F, outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.2),size=7) +
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=7,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 'wilcox.test',
           label = "{ifelse(p < 0.0001, 'p < 0.0001', sprintf('p = %.4f', p))}",
           label.size = 10,
           tip.length = 0.01,
           hide.ns = F)+
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  labs(x='',
       y=expression("Monocyte to epithelial cell ratio"~"("~"Log"[10]~")"),
       title = "d.")+
  scale_x_discrete(labels=c("HIV-"="HIV-",
                            "PLHIV ART <3m"="PLHIV\nART<3m",
                            "PLHIV ART >1y"="PLHIV\nART>1yr"))+
  theme_classic()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        plot.title = element_text(size = 50, face = "bold"),
        axis.ticks.length = unit(0.5,"cm"),
        legend.text = element_text(size=25),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 30))

Fig1d

###################################################################### Figure 1e
Fig1e <- Figure_1_data %>%
  ggplot(aes(x=factor(`HIV Status`,levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y')),
             y=as.numeric(`CD4+`),
             color=factor(`HIV Status`,
                          levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y'))))+
  geom_boxplot(width=0.5,notch = F,outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.2),size=7) +
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=12,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 't.test',
           label = "{ifelse(p < 0.0001, 'p < 0.0001', sprintf('p = %.4f', p))}",
           label.size = 10,
           tip.length = 0.01,
           hide.ns = F)+
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  labs(x='',
       y=expression("Frequency of CD4"^"+"~"T cells (%)"),
       title = "e.")+
  #scale_y_continuous(limits = c(0,60),breaks = c(0,10,20,30,40,50,60))+
  scale_x_discrete(labels=c("HIV-"="HIV-",
                            "PLHIV ART <3m"="PLHIV\nART<3m",
                            "PLHIV ART >1y"="PLHIV\nART>1yr"))+
  theme_classic()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        plot.title = element_text(size = 50, face = "bold"),
        axis.ticks.length = unit(0.5,"cm"),
        legend.text = element_text(size=25),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 30))
Fig1e

###################################################################### Figure 1f
Fig1f <- Figure_1_data %>%
  ggplot(aes(x=factor(`HIV Status`,levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y')),
             y=as.numeric(`CD8+`),
             color=factor(`HIV Status`,
                          levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y'))))+
  geom_boxplot(width=0.5,notch = F,outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.2),size=7) +
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=12,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 't.test',
           label = "{ifelse(p < 0.0001, 'p < 0.0001', sprintf('p = %.4f', p))}",
           label.size = 10,
           tip.length = 0.01,
           hide.ns = F)+
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  labs(x='',
       y=expression("Frequency of CD8"^"+"~"T cells (%)"),
       title = "f.")+
  scale_y_continuous(limits = c(20,120),breaks = c(20,40,60,80,100,120))+
  scale_x_discrete(labels=c("HIV-"="HIV-",
                            "PLHIV ART <3m"="PLHIV\nART<3m",
                            "PLHIV ART >1y"="PLHIV\nART>1yr"))+
  theme_classic()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        plot.title = element_text(size = 50, face = "bold"),
        axis.ticks.length = unit(0.5,"cm"),
        legend.text = element_text(size=25),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 30))

Fig1f

###################################################################### Figure 1g
Fig1g <- Figure_1_data %>%
  ggplot(aes(x=factor(`HIV Status`,levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y')),
             y=as.numeric(`CD3+Mait`),
             color=factor(`HIV Status`,
                          levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y'))))+
  geom_boxplot(width=0.5,notch = F,outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.3),size=7) +
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=12,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 't.test',
           label = "{ifelse(p < 0.0001, 'p < 0.0001', sprintf('p = %.4f', p))}",
           label.size = 10,
           tip.length = 0.01,
           hide.ns = F)+
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  labs(x='',
       y=expression("Frequency of MAIT cells (%)"),
       title = "g.")+
  scale_y_continuous(limits = c(0,40),breaks = c(0,10,20,30,40))+
  scale_x_discrete(labels=c("HIV-"="HIV-",
                            "PLHIV ART <3m"="PLHIV\nART<3m",
                            "PLHIV ART >1y"="PLHIV\nART>1yr"))+
  theme_classic()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        plot.title = element_text(size = 50, face = "bold"),
        axis.ticks.length = unit(0.5,"cm"),
        legend.text = element_text(size=25),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 30))

Fig1g

###################################################################### Figure 1h
Fig1h <- Figure_1_data %>%
  ggplot(aes(x=factor(`HIV Status`,levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y')),
             y=as.numeric(`CD3+TCRgd+`),
             color=factor(`HIV Status`,
                          levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y'))))+
  geom_boxplot(width=0.5,notch = F,outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.3),size=7) +
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=12,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 't.test',
           label = "{ifelse(p < 0.0001, 'p < 0.0001', sprintf('p = %.4f', p))}",
           label.size = 10,
           tip.length = 0.01,
           hide.ns = F)+
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  labs(x='',
       y=expression("Frequency of TCR"~gamma~delta~""^"+"~"T cells"),
       title = "h.")+
  scale_x_discrete(labels=c("HIV-"="HIV-",
                            "PLHIV ART <3m"="PLHIV\nART<3m",
                            "PLHIV ART >1y"="PLHIV\nART>1yr"))+
  theme_classic()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        plot.title = element_text(size = 50, face = "bold"),
        axis.ticks.length = unit(0.5,"cm"),
        legend.text = element_text(size=25),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 30))

Fig1h

###################################################################### Figure 1i
Fig1i <- Figure_1_data %>%
  ggplot(aes(x=factor(`HIV Status`,levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y')),
             y=as.numeric(`CD3+CD56+`),
             color=factor(`HIV Status`,
                          levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y'))))+
  geom_boxplot(width=0.5,notch = F,outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.3),size=7) +
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=12,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 't.test',
           label = "{ifelse(p < 0.0001, 'p < 0.0001', sprintf('p = %.4f', p))}",
           label.size = 10,
           tip.length = 0.01,
           hide.ns = F)+
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  labs(x='',
       y=expression("Frequency of NK T cells (%)"),
       title = "i.")+
  scale_x_discrete(labels=c("HIV-"="HIV-",
                            "PLHIV ART <3m"="PLHIV\nART<3m",
                            "PLHIV ART >1y"="PLHIV\nART>1yr"))+
  theme_classic()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        plot.title = element_text(size = 50, face = "bold"),
        axis.ticks.length = unit(0.5,"cm"),
        legend.text = element_text(size=25),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 30))
Fig1i
#################################### SUPPLEMENTARY FIGURE 1 #########################
############################################################# Supplementary Figure 1a
Supp_Figure_1a_data <- magick::image_read("HIV-PAPER/Data/SPN-PLHIV_2025_data/Supp_Figure_1a_data.png")
Supp_Figure_1a_data <- ggdraw()+ 
  draw_image(Supp_Figure_1a_data, scale = 1.0)+
  labs(title = "a.",
       subtitle = "Myeloid panel gating strategy")+
  theme(
    plot.title = element_text(hjust = 0.01,vjust = .01, size = 50, face = "bold"),  # Customize title
    plot.subtitle = element_text(hjust = 0.5, size = 50),  # Customize title
    plot.margin = unit(c(0,0,0,0), "cm")  # Optional: adjust margins
  )
Supp_Figure_1a_data
############################################################# Supplementary Figure 1b
Supp_Figure_1b_data <- magick::image_read("HIV-PAPER/Data/SPN-PLHIV_2025_data/Supp_Figure_1b_data.png")
Supp_Figure_1b_data <- ggdraw()+ 
  draw_image(Supp_Figure_1b_data, scale = 1.0)+
  labs(title = "b.",
       subtitle = "T cell panel gating strategy")+
  theme(
    plot.title = element_text(hjust = 0.0,vjust = .01, size = 50, face = "bold"),  # Customize title
    plot.subtitle = element_text(hjust = 0.5, size = 50),  # Customize title
    plot.margin = unit(c(0,0,0,0), "cm")  # Optional: adjust margins
  )
Supp_Figure_1b_data

#################################### FIGURE 2 ##################################

Figure_2a_data <- magick::image_read("HIV-PAPER/Data/SPN-PLHIV_2025_data/Figure_2a_data.png")
Figure_2a_data <- ggdraw()+ 
  draw_image(Figure_2a_data, scale = 1.2)+
  labs(title = "a.")+
  theme(
    plot.title = element_text(hjust = -0.1, vjust = 0.1, size = 50, face = "bold"),  # Customize title
    plot.margin = unit(c(0,0,0,0), "cm")  # Optional: adjust margins
  )
Figure_2a_data


Figure_2c_data <- read_csv("HIV-PAPER/Data/SPN-PLHIV_2025_data/Figure_2c_data.csv")
Fig2c <- Figure_2c_data %>%
  ggplot(aes(x=factor(`HIV Status`,levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y')),
             y=log10(as.numeric(`Myeloperoxidase (pg/mL)`)),
             color=factor(`HIV Status`,
                          levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y'))))+
  geom_boxplot(width=0.5,notch = F,outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.1),size=7) +
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=7,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 't.test',
           label = "{ifelse(p < 0.0001, 'p < 0.0001', sprintf('p = %.4f', p))}",
           label.size = 10,
           tip.length = 0.01,
           hide.ns = F)+
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  labs(#x='',
    y=expression("MPO concentration"~"(Log"[10]~"(pg/ml)"),
    title = "c."
  )+
  scale_x_discrete(labels=c("HIV-"="HIV-",
                            "PLHIV ART <3m"="PLHIV\nART<3m",
                            "PLHIV ART >1y"="PLHIV\nART>1y"))+
  theme_classic()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        plot.title = element_text(size = 55, face = "bold", 
                                  hjust = -0.20,
                                  vjust = 0.5),
        legend.text = element_text(size=25),
        axis.ticks.length = unit(0.5,'cm'),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title.y = element_text(size = 30))
Fig2c

# Figure 2d
Figure_2d_data <- read_csv("HIV-PAPER/Data/SPN-PLHIV_2025_data/Figure_2d_data.csv")
Fig2d <- Figure_2d_data %>%
  ggplot(aes(log10(`Myeloperoxidase (pg/mL)`),log10(NER))) +
  geom_point(size=7,aes(color=`HIV Status`))+
  geom_smooth(method = "lm")+
  stat_cor(method = "spearman",size=10)+
  labs(y= expression("Neutrophil to epithelial cell ratio"~"(Log"[10]~")"),
       x=expression("MPO concentration"~"(Log"[10]~"pg/ml)"),title = "d."
  )+
  scale_y_continuous(limits = c(-2,1))+
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  theme_classic()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        plot.title = element_text(size = 55, face = "bold", 
                                  hjust = -0.20,
                                  vjust = 0.5),
        legend.text = element_text(size=25),
        axis.ticks.length = unit(0.5,'cm'),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 30))
Fig2d

# Figure 2e
Figure_2e_data <- read_csv("HIV-PAPER/Data/SPN-PLHIV_2025_data/Figure_2e_data.csv")
Fig2e <- Figure_2e_data %>%
  ggplot(aes(log10(`Myeloperoxidase (pg/mL)`),log10(MER))) +
  geom_point(size=7,aes(color=`HIV Status`))+
  geom_smooth(method = "lm")+
  stat_cor(method = "spearman",size=10)+
  labs(y= expression("Monocyte to epithelial cell ratio"~"(Log"[10]~")"),
       x=expression("MPO concentration"~"(Log"[10]~"pg/ml)"),title = "e."
  )+
  scale_y_continuous(limits = c(-3,0))+
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  theme_classic()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        plot.title = element_text(size = 55, face = "bold", 
                                  hjust = -0.20,
                                  vjust = 0.5),
        legend.text = element_text(size=25),
        axis.ticks.length = unit(0.5,'cm'),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 30))
Fig2e

Figure_2f_to_2h_data <- read_csv("HIV-PAPER/Data/SPN-PLHIV_2025_data/Figure_2f_to_2h_data.csv")
# Cytokine correlations
Cytokine_correlations <- Figure_2f_to_2h_data %>%
  # Select specified columns
  dplyr::select(`TNF-a`, `IL-6`, `IL-8`, `IL-1B`, `IFN-y`, `IL-2`, `IL-12p70`, 
                `IL-13`, `IL-4`, `IL-10`, `Neutrophils`, `Monocytes`, `HIV Status`) %>%
  # Remove rows with NA values
  na.omit() %>%
  # Apply log10 transformation to numeric columns, adding a small constant to avoid log(0)
  dplyr::mutate(across(where(is.numeric), ~ log10(.x + 1e-6)))
# Convert only numeric columns, excludingg the specified colum (HIV Status)  
char_column <- "HIV Status"
Cytokine_correlations <- Cytokine_correlations %>%
  mutate(across(-all_of(char_column), ~ as.numeric(as.character(.))))
# Remove rows with any "0" values across numeric colums
Cytokine_correlations <- Cytokine_correlations %>%
  filter(if_all(-all_of(char_column), ~ . != 0))
# Split the data by HIV Status
Cytokine_correlations <- split(Cytokine_correlations,
                               Cytokine_correlations$`HIV Status`)

# Calculate correlation matrices and p-values
cor_results <- lapply(Cytokine_correlations, function(df) {
  data <- df %>% dplyr::select(-`HIV Status`)
  Hmisc::rcorr(as.matrix(data)) # Calculates correlations and p-values
})

# Set up plot layout: one plot for each HIV status
num_plots <- length(cor_results)
par(mfrow = c(1,3))  # Adjust for the number of statuses
col_grad <- colorRampPalette(c("blue", "white", "red"))

# Plot each correlation matrix with p-values
for (status in names(cor_results)) {
  cor_mat <- cor_results[[status]]$r   # Correlation coefficients
  p_mat <- cor_results[[status]]$P  # P-values
}

png("HIV-PAPER/Data/SPN-PLHIV_2025_data/Fig2f.png", width = 800, height = 800)
# Add correlation plot
col_grad <- colorRampPalette(c("blue", "white", "red"))
corrplot(cor_results$`HIV-`$r,
         method = "circle",
         type = "lower",
         col = col_grad(200),
         title = paste("HIV-"),
         mar = c(0, 0, 2, 0),
         cl.lim = c(-1, 1),
         addgrid.col = "black",
         tl.cex = 1.9,
         tl.col = "black",
         number.cex = 0.7,
         p.mat = cor_results$`HIV-`$P,
         sig.level = 0.05,
         insig = "label_sig")
dev.off()

# Add correlation plot
png("HIV-PAPER/Data/SPN-PLHIV_2025_data/Fig2g.png", width = 800, height = 800)
corrplot(cor_results$`PLHIV ART <3m`$r,
         method = "circle",
         type = "lower",
         col = col_grad(200),
         title = paste("PLHIV ART<3m"),
         mar = c(0, 0, 2, 0),
         cl.lim = c(-1, 1),
         addgrid.col = "black",
         tl.cex = 1.9,
         tl.col = "black",  
         number.cex = 0.7,
         p.mat = cor_results$`PLHIV ART <3m`$P,            # Add p-value matrix
         sig.level = 0.05,        # Highlight significant correlations
         insig = "label_sig")     # Show p-values as labels
dev.off()


# Add correlation plot
png("HIV-PAPER/Data/SPN-PLHIV_2025_data/Fig2h.png", width = 800, height = 800)
corrplot(cor_results$`PLHIV ART >1y`$r,
         method = "circle",
         type = "lower",
         col = col_grad(200),
         title = paste("PLHIV ART>1yr"),
         mar = c(0, 0, 2, 0),
         cl.lim = c(-1, 1),
         addgrid.col = "black",
         tl.cex = 1.9,
         tl.col = "black",  
         #tl.font = 2,
         number.cex = 0.7,
         p.mat = cor_results$`PLHIV ART >1y`$P,            # Add p-value matrix
         sig.level = 0.05,        # Highlight significant correlations
         insig = "label_sig")     # Show p-values as labels
dev.off()

par(mfrow = c(1, 1))

#################################### SUPPLEMENTARY FIGURE 2 ####################
Supp_Figure_2_data <- read_csv("HIV-PAPER/Data/SPN-PLHIV_2025_data/Supp_Figure_2_data.csv")
Supplementary_Fig2b <- Supp_Figure_2_data %>%
  ggplot(aes(x=`HIV Status`, 
             y=log10(`TNF-a`), 
             color=`HIV Status`)) +
  geom_boxplot(width=0.5,notch = F, outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.2),size=7) +
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=7,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 't.test',
           label = "{sprintf('p = %.3f', p)}",
           label.size = 10,
           tip.length = 0.01,
           hide.ns = F)+
  geom_hline(yintercept = log10(0.001823205), linetype='dashed', linewidth=1)+
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  scale_x_discrete(labels=c("HIV-"="HIV-",
                            "PLHIV ART <3m"="PLHIV\nART<3m",
                            "PLHIV ART >1y"="PLHIV\nART>1yr"))+
  labs(y="Cytokine concentration"~"(Log"[10]~"pg/ml)",
       title = "b.",
       subtitle = "TNF-"~alpha,
       x='')+
  theme_classic()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size = 30,face = 'bold'),
        plot.title = element_text(size = 50, face = "bold"),
        plot.subtitle = element_text(size = 30, hjust = .5),
        legend.text = element_text(size = 25),
        axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 23),
        axis.title = element_text(size = 30))
Supplementary_Fig2b

# Figure 2c
Supplementary_Fig2c <- Supp_Figure_2_data %>%
  ggplot(aes(x=`HIV Status`, 
             y=log10(`IL-6`), 
             color=`HIV Status`)) +
  geom_boxplot(width=0.5,notch = F,outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.2),size=7) +
  geom_hline(yintercept = log10(0.005131728), linetype='dashed', linewidth=1)+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=12,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(
    method = 'wilcox.test',
    label = "{sprintf('p = %.3f', p)}",
    label.size = 8,
    tip.length = 0.01,
    hide.ns = F,
    p.adjust.method = "BH"
  )+
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  scale_x_discrete(labels=c("HIV-"="HIV-",
                            "PLHIV ART <3m"="PLHIV\nART<3m",
                            "PLHIV ART >1y"="PLHIV\nART>1yr"))+
  labs(y="Cytokine concentration"~"(Log"[10]~"pg/ml)",
       title = "c.",
       subtitle = "IL-6",
       x='')+
  theme_classic()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        plot.title = element_text(size = 50, face = "bold"),
        plot.subtitle = element_text(size = 30, hjust = .5),
        legend.text = element_text(size=25),
        axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 23),
        axis.title = element_text(size = 30))
Supplementary_Fig2c

# Figure 2d
Supplementary_Fig2d <- Supp_Figure_2_data %>%
  ggplot(aes(x=`HIV Status`, 
             y=log10(`IL-8`), 
             color=`HIV Status`)) +
  geom_boxplot(width=0.5,notch = F,outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.2),size=7) +
  geom_hline(yintercept = log10(7.060539e-02), linetype='dashed', linewidth=1)+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=12,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(
    method = 'wilcox.test',
    label = "{sprintf('p = %.3f', p)}",
    label.size = 8,
    tip.length = 0.01,
    hide.ns = F,
    p.adjust.method = "BH"
  )+
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  scale_x_discrete(labels=c("HIV-"="HIV-",
                            "PLHIV ART <3m"="PLHIV\nART<3m",
                            "PLHIV ART >1y"="PLHIV\nART>1yr"))+
  labs(y="Cytokine concentration"~"(Log"[10]~"pg/ml)",
       title = "d.",
       subtitle = "IL-8",
       x='')+
  theme_classic()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        plot.title = element_text(size = 50, face = "bold"),
        plot.subtitle = element_text(size = 30, hjust = .5),
        legend.text = element_text(size = 25),
        axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 23),
        axis.title = element_text(size = 30))
Supplementary_Fig2d

# Figure 2e
Supplementary_Fig2e <- Supp_Figure_2_data %>%
  ggplot(aes(x=`HIV Status`, 
             y=log10(`IL-1B`), 
             color=`HIV Status`)) +
  geom_boxplot(width=0.5,notch = F,outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.2),size=7) +
  #geom_hline(yintercept = log10(0.000101521), linetype='dashed', linewidth=1)+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=12,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(
    method = 'wilcox.test',
    label = "{sprintf('p = %.3f', p)}",
    label.size = 8,
    tip.length = 0.01,
    hide.ns = F,
    p.adjust.method = "BH"
  )+
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  scale_x_discrete(labels=c("HIV-"="HIV-",
                            "PLHIV ART <3m"="PLHIV\nART<3m",
                            "PLHIV ART >1y"="PLHIV\nART>1yr"))+
  labs(y="Cytokine concentration"~"(Log"[10]~"pg/ml)",
       subtitle = "IL-1"~beta,
       title = "e.",
       x='')+
  theme_classic()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        plot.title = element_text(size = 50, face = "bold"),
        plot.subtitle = element_text(size = 30, hjust = .5),
        legend.text = element_text(size=25),
        axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 23),
        axis.title = element_text(size = 30))
Supplementary_Fig2e

# Figure 2f
Supplementary_Fig2f <- Supp_Figure_2_data %>%
  ggplot(aes(x=`HIV Status`, 
             y=log10(`IFN-y`), 
             color=`HIV Status`)) +
  geom_boxplot(width=0.5,notch = F,outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.2),size=7) +
  geom_hline(yintercept = log10(0.000101521), linetype='dashed', linewidth=1)+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=12,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(
    method = 'wilcox.test',
    label = "{sprintf('p = %.3f', p)}",
    label.size = 8,
    tip.length = 0.01,
    hide.ns = F,
    p.adjust.method = "BH"
  )+
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  scale_x_discrete(labels=c("HIV-"="HIV-",
                            "PLHIV ART <3m"="PLHIV\nART<3m",
                            "PLHIV ART >1y"="PLHIV\nART>1yr"))+
  labs(y="Cytokine concentration"~"(Log"[10]~"pg/ml)",
       subtitle = "IFN-"~gamma,
       title = "f.",
       x='')+
  theme_classic()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        plot.title = element_text(size = 50, face = "bold"),
        plot.subtitle = element_text(size = 30, hjust = .5),
        legend.text = element_text(size=25),
        axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 23),
        axis.title = element_text(size = 30))
Supplementary_Fig2f

# Figure 2g
# Interleukin-2
Supplementary_Fig2g <- Supp_Figure_2_data %>%
  ggplot(aes(x=`HIV Status`, 
             y=log10(`IL-2`), 
             color=`HIV Status`)) +
  geom_boxplot(width=0.5,notch = F,outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.2),size=7) +
  geom_hline(yintercept = log10(0.01794849), linetype='dashed', linewidth=1)+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=12,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(
    method = 'wilcox.test',
    label = "{sprintf('p = %.3f', p)}",
    label.size = 8,
    tip.length = 0.01,
    hide.ns = F,
    p.adjust.method = "BH"
  )+
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  scale_x_discrete(labels=c("HIV-"="HIV-",
                            "PLHIV ART <3m"="PLHIV\nART<3m",
                            "PLHIV ART >1y"="PLHIV\nART>1yr"))+
  labs(y="Cytokine concentration"~"(Log"[10]~"pg/ml)",
       title = "g.",
       subtitle = "IL-2",
       x='')+
  theme_classic()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        plot.title = element_text(size = 50, face = "bold"),
        plot.subtitle = element_text(size = 30, hjust = .5),
        legend.text = element_text(size=25),
        axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 23),
        axis.title = element_text(size = 30))
Supplementary_Fig2g

# Figure 2h

Supplementary_Fig2h <- Supp_Figure_2_data %>%
  ggplot(aes(x=`HIV Status`, 
             y=log10(`IL-12p70`), 
             color=`HIV Status`)) +
  geom_boxplot(width=0.5,notch = F,outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.2),size=7) +
  geom_hline(yintercept = log10(0.000322027), linetype='dashed', linewidth=1)+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=12,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(
    method = 'wilcox.test',
    label = "{sprintf('p = %.3f', p)}",
    label.size = 8,
    tip.length = 0.01,
    hide.ns = F,
    p.adjust.method = "BH"
  )+
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  scale_x_discrete(labels=c("HIV-"="HIV-",
                            "PLHIV ART <3m"="PLHIV\nART<3m",
                            "PLHIV ART >1y"="PLHIV\nART>1yr"))+
  labs(y="Cytokine concentration"~"(Log"[10]~"pg/ml)",
       subtitle = "IL-12p70",
       title = "h.",
       x='')+
  theme_classic()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        plot.title = element_text(size = 50, face = "bold"),
        plot.subtitle = element_text(size = 30, hjust = .5),
        legend.text = element_text(size=25),
        axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 15,face = "bold"),
        axis.title = element_text(size = 30))
Supplementary_Fig2h

# Figure 2i
Supplementary_Fig2i <- Supp_Figure_2_data %>%
  ggplot(aes(x=`HIV Status`, 
             y=log10(`IL-4`), 
             color=`HIV Status`)) +
  geom_boxplot(width=0.5,notch = F,outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.2),size=7) +
  geom_hline(yintercept = log10(0.000143841), linetype='dashed', linewidth=1)+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=12,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(
    method = 'wilcox.test',
    label = "{sprintf('p = %.3f', p)}",
    label.size = 10,
    tip.length = 0.01,
    hide.ns = F,
    p.adjust.method = "BH"
  )+
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  scale_x_discrete(labels=c("HIV-"="HIV-",
                            "PLHIV ART <3m"="PLHIV\nART<3m",
                            "PLHIV ART >1y"="PLHIV\nART>1yr"))+
  labs(y="Cytokine concentration"~"(Log"[10]~"pg/ml)",
       subtitle = "IL-4",
       title="i.",
       x='')+
  theme_classic()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        plot.title = element_text(size = 50, face = "bold"),
        plot.subtitle = element_text(size = 30, hjust = .5),
        legend.text = element_text(size=25),
        axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 15,face = "bold"),
        axis.title = element_text(size = 30))
Supplementary_Fig2i

# Figure 2j
Supplementary_Fig2j <- Supp_Figure_2_data %>%
  ggplot(aes(x=`HIV Status`, 
             y=log10(`IL-13`), 
             color=`HIV Status`)) +
  geom_boxplot(width=0.5,notch = F,outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.2),size=7) +
  geom_hline(yintercept = log10(0.03676782), linetype='dashed', linewidth=1)+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=12,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(
    method = 'wilcox.test',
    label = "{sprintf('p = %.3f', p)}",
    label.size = 8,
    tip.length = 0.01,
    hide.ns = F,
    p.adjust.method = "BH"
  )+
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  scale_x_discrete(labels=c("HIV-"="HIV-",
                            "PLHIV ART <3m"="PLHIV\nART<3m",
                            "PLHIV ART >1y"="PLHIV\nART>1yr"))+
  labs(y="Cytokine concentration"~"(Log"[10]~"pg/ml)",
       subtitle = "IL-13",
       title = "j.",
       x='')+
  theme_classic()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        plot.title = element_text(size = 50, face = "bold"),
        plot.subtitle = element_text(size = 30, hjust = .5),
        legend.text = element_text(size=25),
        axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 15,face = "bold"),
        axis.title = element_text(size = 30))
Supplementary_Fig2j

# Figure 2k
Supplementary_Fig2k <- Supp_Figure_2_data %>%
  ggplot(aes(x=`HIV Status`, 
             y=log10(`IL-10`), 
             color=`HIV Status`)) +
  geom_boxplot(width=0.5,notch = F,outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.2),size=7) +
  geom_hline(yintercept = log10(2.783030e-04), linetype='dashed', linewidth=1)+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=12,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(
    method = 'wilcox.test',
    label = "{sprintf('p = %.3f', p)}",
    label.size = 8,
    tip.length = 0.01,
    hide.ns = F,
    p.adjust.method = "BH"
  )+
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  scale_x_discrete(labels=c("HIV-"="HIV-",
                            "PLHIV ART <3m"="PLHIV\nART<3m",
                            "PLHIV ART >1y"="PLHIV\nART>1yr"))+
  labs(y="Cytokine concentration"~"(Log"[10]~"pg/ml)",
       subtitle = "IL-10",
       title="k.",
       x='')+
  theme_classic()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        plot.title = element_text(size = 50, face = "bold"),
        plot.subtitle = element_text(size = 30, hjust = .5),
        legend.text = element_text(size=25),
        axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 15,face = "bold"),
        axis.title = element_text(size = 30))
Supplementary_Fig2k

#################################### FIGURE 3 ##################################
Figure_3_data <- readRDS("HIV-PAPER/Data/SPN-PLHIV_2025_data/Figure_3def_5abcd_data.rds")
all_merged_subset_labelled_new <- Figure_3_data
# Define a consistent color palette for 17 clusters
all_merged_subset_labelled_new$Clusters <- paste0(all_merged_subset_labelled_new@active.ident)
all_merged_subset_labelled_new <- RenameIdents(all_merged_subset_labelled_new,
                                               "Mono/Mac"="Monocytes/Macrophages")
all_merged_subset_labelled_new$Clusters <- paste0(all_merged_subset_labelled_new@active.ident)
Clusters <- unique(all_merged_subset_labelled_new$Clusters)
num_clusters <- 12
base_palette <- brewer.pal(min(num_clusters,12),'Paired')
names(base_palette) <- Clusters
Fig3a <- DimPlot(
  all_merged_subset_labelled_new,
  reduction = 'umap.harmony',
  label = TRUE,
  repel = TRUE,
  shuffle = TRUE,
  alpha = 2,
  label.size = 10)+
  scale_color_manual(values = base_palette)+
  #NoLegend()+
  labs(x='UMAP-1',
       y='UMAP-2'#,
       #title = "a."
  )+
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = c(0.70,0.2),
        plot.title = element_text(size = 50,face = "bold"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 20),  # Increase cluster names size
        legend.title = element_text(size = 16, face = "bold"),  # Increase legend title size
        legend.key.size = unit(1, "cm"),  # Increase the size of legend keys (points)
        legend.key.height = unit(1, "cm"),  # Height of legend keys
        legend.key.width = unit(1, "cm"))  # Width of legend keys
Fig3a <- as.ggplot(Fig3a)
Fig3a
# Save Figure 3a
ggsave(Fig3a,filename="HIV-PAPER/Figures/Figure 3/Fig3a.png",
       width = 12,height = 9,dpi = 300)
ggsave(Fig3a,filename="HIV-PAPER/Figures/Figure 3/Fig3a.pdf",
       width = 12,height = 9,dpi = 300)

# Figure 3b (Main Bargraph)
# Calculate cluster proportion 
cluster_df <- as.data.frame(table(all_merged_subset_labelled_new$Clusters,
                                  all_merged_subset_labelled_new$HIV_Status))
colnames(cluster_df)<-c('Cluster','HIV Status','Count')
cluster_df$`Wrapped_HIVstatus` <- str_wrap(cluster_df$`HIV Status`, width = 50)
Fig3b <- cluster_df %>%
  ggplot(aes(x=`HIV Status`,
             y=Count, fill=factor(Cluster,levels=c("Goblet cells","Secretory cells",
                                                   "CD3+ T cells","Mono/Mac","B cells",
                                                   "Ciliated cells","Neutrophils","Squamous cells",
                                                   "Deuterosomal cells","Ionocytes","Dendritic cells","Basal cells"))))+
  geom_col(position = "fill")+
  scale_fill_manual(values = base_palette)+
  scale_x_discrete(labels = c("HIV-"="HIV-",
                              "HIV+ ART<3 Months"="PLHIV\nART<3m",
                              "HIV+ ART>1 Year"="PLHIV\nART>1yr"))+
  #scale_x_discrete(labels = cluster_df$Wrapped_HIVstatus)+
  labs(x='',
       y='Frequency of cells(%)',
       fill="Cell cluster",
       title = "b.")+
  theme_bw()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=8,face = 'bold'),
        plot.title = element_text(size = 50,face = "bold",hjust = -0.15),
        legend.text = element_text(size=6),
        axis.text.x = element_text(size = 40),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 40))+
  guides(
    size = guide_legend(
      title = "Cell cluster",
      title.theme = element_text(size = 10),
      keyheight = unit(1, "cm")))
Fig3b

# Save Figure 4b
ggsave(Fig3b,filename="HIV-PAPER/Figures/Figure 3/Fig3b.png",
       width = 10,height = 12,dpi = 300)
ggsave(Fig3b,filename="HIV-PAPER/Figures/Figure 3/Fig3b.pdf",
       width = 10,height = 12,dpi = 300)

# Figure 3c (vlnPlot of celltype markers)
celltype_markers <- c(
  "CYP2F1","SERPINB3","MUC5AC",
  "SPDEF","VMO1","AQP5","BPIFB1",
  "PTPRC","CD3E","CD3D","CD8A","CD4",
  "TYROBP","HLA-DPA1","HLA-DPB1",
  "MS4A1","CD79A","CD19",
  "CAPS","CFAP157","DNAAF1","FOXJ1",
  "G0S2","CXCL8","CSF3R","KRT4",
  "SPRR3","SPRR2A","KRT6A",
  "CDC20B","CCNO","MSMB",
  "CFTR","SCNN1B","RARRES2",
  "CD83","CD86","LYZ","FCER1G",
  "PCP4L1","CSRP2","SERPINF1","ITGA6",
  "KRT5"
)

Fig3c <- Seurat::VlnPlot(all_merged_subset_labelled_new,
                         features = celltype_markers,
                         stack = T,
                         alpha = 2,
                         sort = F,
                         log = F,
                         layer = "data",
                         fill.by = "ident",
                         cols = base_palette)+
  NoLegend()+
  labs(x="SCT normalised Log-transformed expression"#,
       #title = "c."
  )+
  theme(axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_line(),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 30),
        strip.text = element_text(size = 9, face = "plain"),
        strip.text.x = element_text(angle = 80, vjust = 0.5,size = 14))
Fig3c

# Save Figure 3c
ggsave(Fig3c,filename="HIV-PAPER/Figures/Figure 3/Fig3c.png",
       width = 20,height = 6,dpi = 300)
ggsave(Fig3c,filename="HIV-PAPER/Figures/Figure 3/Fig3c.pdf",
       width = 20,height = 6,dpi = 300)

# Figure 3d (Immune cells UMAP)
Figure_3d_data <- readRDS("HIV-PAPER/Data/SPN-PLHIV_2025_data/Figure_5abcd_data.rds")
DimPlot(Figure_3d_data)
Immune_cells <- Figure_3d_data
Immune_cells <- RenameIdents(Immune_cells,
                             "Mono"="Macrophages",
                             "Mac"="Monocytes")
Immune_cells$Clusters <- paste0(Immune_cells@active.ident)
Clusters <- unique(Immune_cells$Clusters)
num_clusters <- 9
base_palette <- brewer.pal(min(num_clusters,9),'Paired')
names(base_palette) <- Clusters
Fig3d <- Seurat::DimPlot(
  Immune_cells, 
  reduction = "umap",
  label = TRUE,
  repel = TRUE,
  shuffle = TRUE,
  alpha = 2,
  label.size = 10)+
  labs(x="UMAP-1",
       y="UMAP-2"
  )+
  scale_color_manual(values = base_palette)+
  theme_bw()+
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 50,face = 'bold'),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 20),  # Increase cluster names size
        legend.title = element_text(size = 16, face = "bold"),  # Increase legend title size
        legend.key.size = unit(1, "cm"),  # Increase the size of legend keys (points)
        legend.key.height = unit(1, "cm"),  # Height of legend keys
        legend.key.width = unit(1, "cm"))  # Width of legend keys)
Fig3d <- as.ggplot(Fig3d)
Fig3d
# Save Figure 3d
ggsave(Fig3d,filename="HIV-PAPER/Figures/Figure 3/Fig3d.png",
       width = 8,height = 8,dpi = 300)
ggsave(Fig3d,filename="HIV-PAPER/Figures/Figure 3/Fig3d.pdf",
       width = 8,height = 8,dpi = 300)

# Figure 3e (Immune cell proportion bargraph)
# Calculate cluster proportion 
Immune_cluster_df <- as.data.frame(table(Immune_cells$Clusters,
                                         Immune_cells$HIV_Status))
colnames(Immune_cluster_df)<-c('Cluster','HIV Status','Count')
Immune_cluster_df$`Wrapped_HIVstatus` <- str_wrap(Immune_cluster_df$`HIV Status`, width = 50)
Fig3e <- Immune_cluster_df %>%
  ggplot(aes(x=`HIV Status`,
             y=Count, fill=factor(Cluster,levels=c("Mast cells","Monocytes","B cells","DCs",
                                                   "Macrophages","Plasma B cells", "Neutrophils",
                                                   "NK T cells","CD8+ T cells"))))+
  geom_col(position = "fill")+
  scale_fill_manual(values = base_palette)+
  scale_x_discrete(labels = c("HIV-"="HIV-",
                              "HIV+ ART<3 Months"="PLHIV\nART<3m",
                              "HIV+ ART>1 Year"="PLHIV\nART>1yr"))+
  #scale_x_discrete(labels = cluster_df$Wrapped_HIVstatus)+
  labs(x='',
       y='Frequency of immune cells (%)',
       fill="Cell cluster",title = "e."
  )+
  theme_bw()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=8,face = 'bold'),
        legend.text = element_text(size=6),
        axis.text.x = element_text(size=40),
        plot.title = element_text(size = 50,face = 'bold', hjust = -0.15),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 40))+
  guides(
    size = guide_legend(
      title = "Cell cluster",
      title.theme = element_text(size = 10),
      keyheight = unit(1, "cm")))
Fig3e

# Save Figure 3e
ggsave(Fig3e,filename="HIV-PAPER/Figures/Figure 3/Fig3e.png",
       width = 10,height = 12,dpi = 300)
ggsave(Fig3e,filename="HIV-PAPER/Figures/Figure 3/Fig3e.pdf",
       width = 10,height = 12,dpi = 300)

# Figure 3g (Neutrophil epithelial ratio from scRNA sequencing data)
cluster_df <- as.data.frame(table(all_merged_subset_labelled_new$Clusters,
                                  all_merged_subset_labelled_new$sample))


colnames(cluster_df) <- c("Clusters","sample","Count")
cluster_df <- cluster_df %>%
  tidyr::pivot_wider(names_from = "Clusters",
                     values_from = "Count") %>%
  dplyr::mutate(Epithelial_cells=`Goblet cells`+`Squamous cells`+`Ciliated cells`+
                  `Ionocytes`+`Basal cells`+`Secretory cells`,
                `Neut:Epithelial`=Neutrophils/Epithelial_cells,
                `Tcell:Goblet` = `CD3+ T cells`/`Goblet cells`) %>%
  dplyr::mutate(`HIV Status`=ifelse(grepl("CUF130",sample),"ART>1y",
                                    ifelse(grepl("CUF131",sample),"HIV-",
                                           ifelse(grepl("CUF134",sample),"ART<3m",
                                                  ifelse(grepl("CUF135",sample),"ART>1y",
                                                         ifelse(grepl("CUF137",sample),"HIV-",
                                                                ifelse(grepl("CUF13J",sample),"ART<3m",
                                                                       ifelse(grepl("CUF13K",sample),"ART<3m",
                                                                              ifelse(grepl("CUF136",sample),"HIV-",
                                                                                     ifelse(grepl("CUF13I",sample),"ART>1y",
                                                                                            ifelse(grepl("CUG11X",sample),"ART>1y","ART>1y")))))))))))
# Neutrophil to goblet cell ratio
Fig3f <- cluster_df %>%
  ggplot(aes(factor(`HIV Status`,
                    levels = c("HIV-",
                               "ART<3m",
                               "ART>1y")),
             log10(`Neut:Epithelial`),
             color=factor(`HIV Status`,
                          levels = c("HIV-",
                                     "ART<3m",
                                     "ART>1y"))))+
  geom_boxplot(width=0.5,notch = F, outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.2),size=7) +
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=7,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 't.test',
           label = "{ifelse(p < 0.0001, 'p < 0.0001', sprintf('p = %.4f', p))}",
           label.size = 10,
           tip.length = 0.01,
           hide.ns = F)+
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  labs(x='',
       y=expression("Neutrophil to epithelial cell ratio"~"(log"[10]~")"),
       title = "f.")+
  scale_x_discrete(labels=c("HIV-"="HIV-",
                            "ART<3m"="PLHIV\nART<3m",
                            "ART>1y"="PLHIV\nART>1yr"))+
  theme_classic()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        axis.ticks.length = unit(0.5,'cm'),
        legend.text = element_text(size=25), 
        plot.title = element_text(size = 50,face = "bold",hjust = -0.15),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 30))
Fig3f
ggsave(Fig3f,filename="HIV-PAPER/Figures/Figure 3/Fig3f.png",
       width = 14,height = 18,dpi = 300)
ggsave(Fig3f,filename="HIV-PAPER/Figures/Figure 3/Fig3f.pdf",
       width = 14,height = 18,dpi = 300)
#################################### SUPPLEMENTARY FIGURE 3 ####################
Supp_Figure_3_data <- readRDS("HIV-PAPER/Data/SPN-PLHIV_2025_data/Figure_3abc_4_supp_fig3_data.rds")
Supp_Figure_3_data <- FindVariableFeatures(Supp_Figure_3_data,
                                           selection.method = 'vst',
                                           verbose = TRUE)
genes <- VariableFeatures(Supp_Figure_3_data)
toplot <- SeuratExtend::CalcStats(Supp_Figure_3_data,
                                  features = genes,
                                  method = 'zscore',
                                  order = 'p',
                                  n = 5)

Supplementary_Figure_3a <- Seurat::VlnPlot(Supp_Figure_3_data,
                                  features = rownames(toplot),
                                  stack = T,
                                  alpha = 2,
                                  sort = F,
                                  log = F,
                                  layer = "data",
                                  fill.by = "ident"
)+
  NoLegend()+
  labs(x="SCT normalised Log-transformed expression"#,
       #title = "c."
  )+
  theme(axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_line(),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 30),
        strip.text = element_text(size = 9, face = "plain"),
        strip.text.x = element_text(angle = 80, vjust = 0.5,size = 14))
Supplementary_Figure_3a
#################################### FIGURE 4a #################################
####### Cell to cell communication between  epithelial cells and to nasal neutrophils
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(multinichenetr)
library(nichenetr)
library(Seurat)
library(gridGraphics)
library(patchwork)
# Load the data
Figure_4_data <- readRDS("HIV-PAPER/Data/SPN-PLHIV_2025_data/Figure_3_data.rds")
all_merged_subset_labelled_new <- Figure_4_data
# Load the organism
organism = "human"
if(organism == "human"){
  
  lr_network_all = 
    readRDS(url(
      "https://zenodo.org/record/10229222/files/lr_network_human_allInfo_30112033.rds"
    )) %>% 
    mutate(
      ligand = convert_alias_to_symbols(ligand, organism = organism), 
      receptor = convert_alias_to_symbols(receptor, organism = organism))
  
  lr_network_all = lr_network_all  %>% 
    mutate(ligand = make.names(ligand), receptor = make.names(receptor)) 
  
  lr_network = lr_network_all %>% 
    distinct(ligand, receptor)
  
  ligand_target_matrix = readRDS(url(
    "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"
  ))
  
  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  
  lr_network = lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
  ligand_target_matrix = ligand_target_matrix[, lr_network$ligand %>% unique()]
  
} else if(organism == "mouse"){
  
  lr_network_all = readRDS(url(
    "https://zenodo.org/record/10229222/files/lr_network_mouse_allInfo_30112033.rds"
  )) %>% 
    mutate(
      ligand = convert_alias_to_symbols(ligand, organism = organism), 
      receptor = convert_alias_to_symbols(receptor, organism = organism))
  
  lr_network_all = lr_network_all  %>% 
    mutate(ligand = make.names(ligand), receptor = make.names(receptor)) 
  lr_network = lr_network_all %>% 
    distinct(ligand, receptor)
  
  ligand_target_matrix = readRDS(url(
    "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"
  ))
  
  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  
  lr_network = lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
  ligand_target_matrix = ligand_target_matrix[, lr_network$ligand %>% unique()]
  
}

# Prepare the seurat object for conversion to single cell experiment object
all_merged_subset_labelled_new$Clusters <- paste0(all_merged_subset_labelled_new@active.ident)
all_merged_subset_labelled_new$HIV_Status <- paste0(all_merged_subset_labelled_new$HIV_Status)
all_merged_subset_labelled_new$sample <- paste0(all_merged_subset_labelled_new$sample)
# Make names syntactically valid
all_merged_subset_labelled_new$HIV_Status <- make.names(all_merged_subset_labelled_new$HIV_Status)
all_merged_subset_labelled_new$sample <- make.names(all_merged_subset_labelled_new$sample)
all_merged_subset_labelled_new$Clusters <- make.names(all_merged_subset_labelled_new$Clusters)
# Convert to a single cell experiment
sce <- Seurat::as.SingleCellExperiment(all_merged_subset_labelled_new)
SingleCellExperiment::colData(sce)
# Prepare cell-cell communication analysis
sample_id="sample"
group_id="HIV_Status"
celltype_id="Clusters"
covariates=NA
batches=NA
# Define contracts
contrasts_oi = c("'HIV..ART.3.Months-(HIV..ART.1.Year+HIV.)/2','HIV..ART.1.Year-(HIV..ART.3.Months+HIV.)/2','HIV.-(HIV..ART.1.Year+HIV..ART.3.Months)/2'")
contrast_tbl = tibble(contrast =
                        c("HIV..ART.3.Months-(HIV..ART.1.Year+HIV.)/2","HIV..ART.1.Year-(HIV..ART.3.Months+HIV.)/2","HIV.-(HIV..ART.1.Year+HIV..ART.3.Months)/2"),
                      group = c("HIV..ART.3.Months","HIV..ART.1.Year","HIV."))

senders_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
receivers_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
sce = sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in%
            c(senders_oi,receivers_oi)]

conditions_keep = c("HIV..ART.3.Months","HIV..ART.1.Year","HIV.")
sce = sce[,SummarizedExperiment::colData(sce)[,group_id] %in% 
            conditions_keep]
# Running Multinichenet core analysis
# Cell-type filtering
min_cells = 3
abundance_info =multinichenetr::get_abundance_info(
  sce = sce, 
  sample_id = sample_id, 
  group_id = group_id, 
  celltype_id = celltype_id, 
  min_cells = min_cells, 
  senders_oi = senders_oi, 
  receivers_oi = receivers_oi, 
  batches = batches)

abundance_info$abund_plot_sample

abundance_df_summarized = abundance_info$abundance_data %>%
  dplyr::mutate(keep = as.logical(keep)) %>%
  dplyr::group_by(group_id,celltype_id) %>%
  dplyr::summarise(samples_present = sum((keep)))

celltypes_absent_one_condition = abundance_df_summarized %>%
  dplyr::filter(samples_present == 0) %>% 
  dplyr::pull(celltype_id) %>%
  unique()

celltypes_present_one_condition = abundance_df_summarized %>%
  dplyr::filter(samples_present >= 2) %>%
  dplyr::pull(celltype_id) %>%
  unique()

total_nr_conditions = SummarizedExperiment::colData(sce)[,group_id] %>%
  unique() %>%
  length()

absent_celltypes <- abundance_df_summarized %>%
  dplyr::filter(samples_present < 2) %>% 
  dplyr::group_by(celltype_id) %>%
  dplyr::summarize(n = n(), .groups = "drop") %>%
  dplyr::filter(n == total_nr_conditions) %>%
  pull(celltype_id)

analyse_condition_specific_celltypes = TRUE
if(analyse_condition_specific_celltypes == TRUE){
  senders_oi = senders_oi %>% dplyr::setdiff(absent_celltypes)
  receivers_oi = receivers_oi %>% dplyr::setdiff(absent_celltypes)
} else {
  senders_oi = senders_oi %>% 
    dplyr::setdiff(union(absent_celltypes, condition_specific_celltypes))
  receivers_oi = receivers_oi %>% 
    dplyr::setdiff(union(absent_celltypes, condition_specific_celltypes))
}

sce = sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% 
            c(senders_oi, receivers_oi)
]

min_sample_prop = 0.50
fraction_cutoff = 0.05

frq_list = get_frac_exprs(
  sce = sce, 
  sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id, 
  batches = batches, 
  min_cells = min_cells, 
  fraction_cutoff = fraction_cutoff, min_sample_prop = min_sample_prop)

genes_oi = frq_list$expressed_df %>% 
  filter(expressed == TRUE) %>% dplyr::pull(gene) %>% unique() 
sce = sce[genes_oi, ]

# Pseudobulk expression calculation
abundance_expression_info = multinichenetr::process_abundance_expression_info(
  sce = sce, 
  sample_id = sample_id, 
  group_id = group_id, 
  celltype_id = celltype_id, 
  min_cells = min_cells, 
  senders_oi = senders_oi, 
  receivers_oi = receivers_oi, 
  lr_network = lr_network, 
  batches = batches, 
  frq_list = frq_list, 
  abundance_info = abundance_info)

abundance_expression_info$celltype_info$pb_df %>% head()
abundance_expression_info$celltype_info$pb_df_group %>% head()
abundance_expression_info$sender_receiver_info$pb_df %>% head()
abundance_expression_info$sender_receiver_info$pb_df_group %>% head()

# Differential expression (DE) analysis
DE_info = multinichenetr::get_DE_info(
  sce = sce, 
  sample_id = sample_id, 
  group_id = group_id, 
  celltype_id = celltype_id, 
  batches = batches, covariates = covariates, 
  contrasts_oi = contrasts_oi, 
  min_cells = min_cells, 
  expressed_df = frq_list$expressed_df)
# Check DE results
DE_info$celltype_de$de_output_tidy %>% head()
DE_info$hist_pvals


empirical_pval = FALSE
if(empirical_pval == TRUE){
  DE_info_emp = get_empirical_pvals(DE_info$celltype_de$de_output_tidy)
  celltype_de = DE_info_emp$de_output_tidy_emp %>% select(-p_val, -p_adj) %>% 
    rename(p_val = p_emp, p_adj = p_adj_emp)
} else {
  celltype_de = DE_info$celltype_de$de_output_tidy
} 

# Combine DE information for ligand-senders and receptors-receivers
sender_receiver_de = multinichenetr::combine_sender_receiver_de(
  sender_de = celltype_de,
  receiver_de = celltype_de,
  senders_oi = senders_oi,
  receivers_oi = receivers_oi,
  lr_network = lr_network)
sender_receiver_de %>% head(20)

# Ligand activity prediction: use DE analysis output to predict the activity
# of ligands in receiver cell types and infer their potential targer genes
logFC_threshold = 0.50
p_val_threshold = 0.05
p_val_adj = FALSE 
geneset_assessment = contrast_tbl$contrast %>% 
  lapply(
    process_geneset_data, 
    celltype_de, logFC_threshold, p_val_adj, p_val_threshold
  ) %>% 
  bind_rows() 
geneset_assessment

# In case i want to use adjusted p values
geneset_assessment_adjustedPval = contrast_tbl$contrast %>% 
  lapply(
    process_geneset_data, 
    celltype_de, logFC_threshold, p_val_adj = TRUE, p_val_threshold
  ) %>% 
  bind_rows() 
geneset_assessment_adjustedPval

# Perform the ligand activity analysis and ligand-target inference
top_n_target = 250

verbose = TRUE
cores_system = 1
n.cores = min(cores_system, celltype_de$cluster_id %>% unique() %>% length()) 

# Running the ligand activity prediction
ligand_activities_targets_DEgenes = suppressMessages(suppressWarnings(
  multinichenetr::get_ligand_activities_targets_DEgenes(
    receiver_de = celltype_de,
    receivers_oi = intersect(receivers_oi, celltype_de$cluster_id %>% unique()),
    ligand_target_matrix = ligand_target_matrix,
    logFC_threshold = logFC_threshold,
    p_val_threshold = p_val_threshold,
    p_val_adj = p_val_adj,
    top_n_target = top_n_target,
    verbose = verbose, 
    n.cores = n.cores)))

ligand_activities_targets_DEgenes$ligand_activities %>% head(20)

# Prioritization: Rank cell-cell communication patters through multi-criteria prioritization

ligand_activity_down = FALSE
sender_receiver_tbl = sender_receiver_de %>% dplyr::distinct(sender, receiver)

metadata_combined = SummarizedExperiment::colData(sce) %>% tibble::as_tibble()

if(!is.na(batches)){
  grouping_tbl = metadata_combined[,c(sample_id, group_id, batches)] %>% 
    tibble::as_tibble() %>% distinct()
  colnames(grouping_tbl) = c("sample","group",batches)
} else {
  grouping_tbl = metadata_combined[,c(sample_id, group_id)] %>% 
    tibble::as_tibble() %>% distinct()
  colnames(grouping_tbl) = c("sample","group")
}

prioritization_tables = suppressMessages(multinichenetr::generate_prioritization_tables(
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  contrast_tbl = contrast_tbl,
  sender_receiver_tbl = sender_receiver_tbl,
  grouping_tbl = grouping_tbl,
  scenario = "regular", # all prioritization criteria will be weighted equally
  fraction_cutoff = fraction_cutoff, 
  abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
  abundance_data_sender = abundance_expression_info$abundance_data_sender,
  ligand_activity_down = ligand_activity_down))
# Check output tables
prioritization_tables$group_prioritization_tbl %>% head(20)

# Calculate cross-samples expression correlation between ligand-receptor pairs and target genes
lr_target_prior_cor = multinichenetr::lr_target_prior_cor_inference(
  receivers_oi = prioritization_tables$group_prioritization_tbl$receiver %>% unique(), 
  abundance_expression_info = abundance_expression_info, 
  celltype_de = celltype_de, 
  grouping_tbl = grouping_tbl, 
  prioritization_tables = prioritization_tables, 
  ligand_target_matrix = ligand_target_matrix, 
  logFC_threshold = logFC_threshold, 
  p_val_threshold = p_val_threshold, 
  p_val_adj = p_val_adj)

# Save the output of multinichenetr
Path = "scRNAseq_Results/"
multinichenet_output = list(
  celltype_info = abundance_expression_info$celltype_info,
  celltype_de = celltype_de,
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de =  sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  prioritization_tables = prioritization_tables,
  grouping_tbl = grouping_tbl,
  lr_target_prior_cor = lr_target_prior_cor) 
multinichenet_output = multinichenetr::make_lite_output(multinichenet_output)

save = TRUE
if(save == TRUE){
  saveRDS(multinichenet_output, paste0(Path, "multinichenet_output.rds"))
  
}

multinichenet_output <- readRDS("scRNAseq_Results/multinichenet_output.rds")

# Visualization of differential cell-cell interactions
prioritized_tbl_oi_all = multinichenetr::get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  top_n = 250, 
  rank_per_group = T)

prioritized_tbl_oi = 
  multinichenet_output$prioritization_tables$group_prioritization_tbl %>%
  filter(id %in% prioritized_tbl_oi_all$id) %>%
  distinct(id, sender, receiver, ligand, receptor, group) %>% 
  left_join(prioritized_tbl_oi_all)
prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] = 0

Senders <- c("Basal.cells","Ciliated.cells","Goblet.cells","Neutrophils","Secretory.cells","Squamous.cells")
Receivers <- ("Neutrophils")

Filtered_prioritization_tbl_oi <- prioritized_tbl_oi %>%
  dplyr::filter(sender %in% Senders & receiver %in% Receivers)

senders_receivers = union(Filtered_prioritization_tbl_oi$sender %>% unique(), Filtered_prioritization_tbl_oi$receiver %>% unique()) %>% sort()

colors_sender = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)
colors_receiver = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)

circos_list = multinichenetr::make_circos_group_comparison(Filtered_prioritization_tbl_oi, colors_sender, colors_receiver)



# Save the output of Circos_list
Fig4a1_circos <- circos_list$HIV.
Fig4a3_circos <- circos_list$HIV..ART.1.Year
Fig4a2_circos <- circos_list$HIV..ART.3.Months


recordedplot_to_grob <- function(recorded_plot) {
  gridGraphics::grid.echo(recorded_plot)
  grid::grid.grab()
}

Fig4a1_circos <- recordedplot_to_grob(Fig4a1_circos)
Fig4a2_circos <- recordedplot_to_grob(Fig4a2_circos)
Fig4a3_circos <- recordedplot_to_grob(Fig4a3_circos)


# Example: Convert grob to ggplot
ggplot_from_grob <- function(grob) {
  ggplot() +
    annotation_custom(grob) +
    theme_void()  # Remove axes, grid lines, etc.
}

# Convert your grob
Fig4a1_circos <- ggplot_from_grob(Fig4a1_circos)
Fig4a2_circos <- ggplot_from_grob(Fig4a2_circos)
Fig4a3_circos <- ggplot_from_grob(Fig4a3_circos)

# Combine Multiple ggplot Objects
Fig4a <- Fig4a1_circos + Fig4a2_circos + Fig4a3_circos + plot_layout(ncol = 3)
print(Fig4a)

#################################### SUPPLEMENTARY FIGURE 4a ####################
####### Cell to cell communication between  epithelial cells and to nasal neutrophils
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(multinichenetr)
library(nichenetr)
library(Seurat)
library(gridGraphics)
library(patchwork)

# Load the data
Supp_figure_4_data <- readRDS("HIV-PAPER/Data/SPN-PLHIV_2025_data/Supp_figure_4_data.rds")
Immune_cells <- Supp_figure_4_data
# Load the organism
organism = "human"
if(organism == "human"){
  
  lr_network_all = 
    readRDS(url(
      "https://zenodo.org/record/10229222/files/lr_network_human_allInfo_30112033.rds"
    )) %>% 
    mutate(
      ligand = convert_alias_to_symbols(ligand, organism = organism), 
      receptor = convert_alias_to_symbols(receptor, organism = organism))
  
  lr_network_all = lr_network_all  %>% 
    mutate(ligand = make.names(ligand), receptor = make.names(receptor)) 
  
  lr_network = lr_network_all %>% 
    distinct(ligand, receptor)
  
  ligand_target_matrix = readRDS(url(
    "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"
  ))
  
  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  
  lr_network = lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
  ligand_target_matrix = ligand_target_matrix[, lr_network$ligand %>% unique()]
  
} else if(organism == "mouse"){
  
  lr_network_all = readRDS(url(
    "https://zenodo.org/record/10229222/files/lr_network_mouse_allInfo_30112033.rds"
  )) %>% 
    mutate(
      ligand = convert_alias_to_symbols(ligand, organism = organism), 
      receptor = convert_alias_to_symbols(receptor, organism = organism))
  
  lr_network_all = lr_network_all  %>% 
    mutate(ligand = make.names(ligand), receptor = make.names(receptor)) 
  lr_network = lr_network_all %>% 
    distinct(ligand, receptor)
  
  ligand_target_matrix = readRDS(url(
    "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"
  ))
  
  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  
  lr_network = lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
  ligand_target_matrix = ligand_target_matrix[, lr_network$ligand %>% unique()]
  
}

# Prepare the seurat object for conversion to single cell experiment object
Immune_cells$Clusters <- paste0(Immune_cells@active.ident)
Immune_cells$HIV_Status <- paste0(Immune_cells$HIV_Status)
Immune_cells$sample <- paste0(Immune_cells$sample)
# Make names syntactically valid
Immune_cells$HIV_Status <- make.names(Immune_cells$HIV_Status)
Immune_cells$sample <- make.names(Immune_cells$sample)
Immune_cells$Clusters <- make.names(Immune_cells$Clusters)
# Convert to a single cell experiment
sce <- Seurat::as.SingleCellExperiment(Immune_cells)
SingleCellExperiment::colData(sce)
# Prepare cell-cell communication analysis
sample_id="sample"
group_id="HIV_Status"
celltype_id="Clusters"
covariates=NA
batches=NA
# Define contracts
contrasts_oi = c("'HIV..ART.3.Months-(HIV..ART.1.Year+HIV.)/2','HIV..ART.1.Year-(HIV..ART.3.Months+HIV.)/2','HIV.-(HIV..ART.1.Year+HIV..ART.3.Months)/2'")
contrast_tbl = tibble(contrast =
                        c("HIV..ART.3.Months-(HIV..ART.1.Year+HIV.)/2","HIV..ART.1.Year-(HIV..ART.3.Months+HIV.)/2","HIV.-(HIV..ART.1.Year+HIV..ART.3.Months)/2"),
                      group = c("HIV..ART.3.Months","HIV..ART.1.Year","HIV."))

senders_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
receivers_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
sce = sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in%
            c(senders_oi,receivers_oi)]

conditions_keep = c("HIV..ART.3.Months","HIV..ART.1.Year","HIV.")
sce = sce[,SummarizedExperiment::colData(sce)[,group_id] %in% 
            conditions_keep]
# Running Multinichenet core analysis
# Cell-type filtering
min_cells = 3
abundance_info =multinichenetr::get_abundance_info(
  sce = sce, 
  sample_id = sample_id, 
  group_id = group_id, 
  celltype_id = celltype_id, 
  min_cells = min_cells, 
  senders_oi = senders_oi, 
  receivers_oi = receivers_oi, 
  batches = batches)

abundance_info$abund_plot_sample

abundance_df_summarized = abundance_info$abundance_data %>%
  dplyr::mutate(keep = as.logical(keep)) %>%
  dplyr::group_by(group_id,celltype_id) %>%
  dplyr::summarise(samples_present = sum((keep)))

celltypes_absent_one_condition = abundance_df_summarized %>%
  dplyr::filter(samples_present == 0) %>% 
  dplyr::pull(celltype_id) %>%
  unique()

celltypes_present_one_condition = abundance_df_summarized %>%
  dplyr::filter(samples_present >= 2) %>%
  dplyr::pull(celltype_id) %>%
  unique()

total_nr_conditions = SummarizedExperiment::colData(sce)[,group_id] %>%
  unique() %>%
  length()

absent_celltypes <- abundance_df_summarized %>%
  dplyr::filter(samples_present < 2) %>% 
  dplyr::group_by(celltype_id) %>%
  dplyr::summarize(n = n(), .groups = "drop") %>%
  dplyr::filter(n == total_nr_conditions) %>%
  pull(celltype_id)

analyse_condition_specific_celltypes = TRUE
if(analyse_condition_specific_celltypes == TRUE){
  senders_oi = senders_oi %>% dplyr::setdiff(absent_celltypes)
  receivers_oi = receivers_oi %>% dplyr::setdiff(absent_celltypes)
} else {
  senders_oi = senders_oi %>% 
    dplyr::setdiff(union(absent_celltypes, condition_specific_celltypes))
  receivers_oi = receivers_oi %>% 
    dplyr::setdiff(union(absent_celltypes, condition_specific_celltypes))
}

sce = sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% 
            c(senders_oi, receivers_oi)
]

min_sample_prop = 0.50
fraction_cutoff = 0.05

frq_list = get_frac_exprs(
  sce = sce, 
  sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id, 
  batches = batches, 
  min_cells = min_cells, 
  fraction_cutoff = fraction_cutoff, min_sample_prop = min_sample_prop)

genes_oi = frq_list$expressed_df %>% 
  filter(expressed == TRUE) %>% dplyr::pull(gene) %>% unique() 
sce = sce[genes_oi, ]

# Pseudobulk expression calculation
abundance_expression_info = multinichenetr::process_abundance_expression_info(
  sce = sce, 
  sample_id = sample_id, 
  group_id = group_id, 
  celltype_id = celltype_id, 
  min_cells = min_cells, 
  senders_oi = senders_oi, 
  receivers_oi = receivers_oi, 
  lr_network = lr_network, 
  batches = batches, 
  frq_list = frq_list, 
  abundance_info = abundance_info)

abundance_expression_info$celltype_info$pb_df %>% head()
abundance_expression_info$celltype_info$pb_df_group %>% head()
abundance_expression_info$sender_receiver_info$pb_df %>% head()
abundance_expression_info$sender_receiver_info$pb_df_group %>% head()

# Differential expression (DE) analysis
DE_info = multinichenetr::get_DE_info(
  sce = sce, 
  sample_id = sample_id, 
  group_id = group_id, 
  celltype_id = celltype_id, 
  batches = batches, covariates = covariates, 
  contrasts_oi = contrasts_oi, 
  min_cells = min_cells, 
  expressed_df = frq_list$expressed_df)
# Check DE results
DE_info$celltype_de$de_output_tidy %>% head()
DE_info$hist_pvals


empirical_pval = FALSE
if(empirical_pval == TRUE){
  DE_info_emp = get_empirical_pvals(DE_info$celltype_de$de_output_tidy)
  celltype_de = DE_info_emp$de_output_tidy_emp %>% select(-p_val, -p_adj) %>% 
    rename(p_val = p_emp, p_adj = p_adj_emp)
} else {
  celltype_de = DE_info$celltype_de$de_output_tidy
} 

# Combine DE information for ligand-senders and receptors-receivers
sender_receiver_de = multinichenetr::combine_sender_receiver_de(
  sender_de = celltype_de,
  receiver_de = celltype_de,
  senders_oi = senders_oi,
  receivers_oi = receivers_oi,
  lr_network = lr_network)
sender_receiver_de %>% head(20)

# Ligand activity prediction: use DE analysis output to predict the activity
# of ligands in receiver cell types and infer their potential targer genes
logFC_threshold = 0.50
p_val_threshold = 0.05
p_val_adj = FALSE 
geneset_assessment = contrast_tbl$contrast %>% 
  lapply(
    process_geneset_data, 
    celltype_de, logFC_threshold, p_val_adj, p_val_threshold
  ) %>% 
  bind_rows() 
geneset_assessment

# In case i want to use adjusted p values
geneset_assessment_adjustedPval = contrast_tbl$contrast %>% 
  lapply(
    process_geneset_data, 
    celltype_de, logFC_threshold, p_val_adj = TRUE, p_val_threshold
  ) %>% 
  bind_rows() 
geneset_assessment_adjustedPval

# Perform the ligand activity analysis and ligand-target inference
top_n_target = 250

verbose = TRUE
cores_system = 1
n.cores = min(cores_system, celltype_de$cluster_id %>% unique() %>% length()) 

# Running the ligand activity prediction
ligand_activities_targets_DEgenes = suppressMessages(suppressWarnings(
  multinichenetr::get_ligand_activities_targets_DEgenes(
    receiver_de = celltype_de,
    receivers_oi = intersect(receivers_oi, celltype_de$cluster_id %>% unique()),
    ligand_target_matrix = ligand_target_matrix,
    logFC_threshold = logFC_threshold,
    p_val_threshold = p_val_threshold,
    p_val_adj = p_val_adj,
    top_n_target = top_n_target,
    verbose = verbose, 
    n.cores = n.cores)))

ligand_activities_targets_DEgenes$ligand_activities %>% head(20)

# Prioritization: Rank cell-cell communication patters through multi-criteria prioritization

ligand_activity_down = FALSE
sender_receiver_tbl = sender_receiver_de %>% dplyr::distinct(sender, receiver)

metadata_combined = SummarizedExperiment::colData(sce) %>% tibble::as_tibble()

if(!is.na(batches)){
  grouping_tbl = metadata_combined[,c(sample_id, group_id, batches)] %>% 
    tibble::as_tibble() %>% distinct()
  colnames(grouping_tbl) = c("sample","group",batches)
} else {
  grouping_tbl = metadata_combined[,c(sample_id, group_id)] %>% 
    tibble::as_tibble() %>% distinct()
  colnames(grouping_tbl) = c("sample","group")
}

prioritization_tables = suppressMessages(multinichenetr::generate_prioritization_tables(
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  contrast_tbl = contrast_tbl,
  sender_receiver_tbl = sender_receiver_tbl,
  grouping_tbl = grouping_tbl,
  scenario = "regular", # all prioritization criteria will be weighted equally
  fraction_cutoff = fraction_cutoff, 
  abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
  abundance_data_sender = abundance_expression_info$abundance_data_sender,
  ligand_activity_down = ligand_activity_down))
# Check output tables
prioritization_tables$group_prioritization_tbl %>% head(20)

# Calculate cross-samples expression correlation between ligand-receptor pairs and target genes
lr_target_prior_cor = multinichenetr::lr_target_prior_cor_inference(
  receivers_oi = prioritization_tables$group_prioritization_tbl$receiver %>% unique(), 
  abundance_expression_info = abundance_expression_info, 
  celltype_de = celltype_de, 
  grouping_tbl = grouping_tbl, 
  prioritization_tables = prioritization_tables, 
  ligand_target_matrix = ligand_target_matrix, 
  logFC_threshold = logFC_threshold, 
  p_val_threshold = p_val_threshold, 
  p_val_adj = p_val_adj)

# Save the output of multinichenetr
Path = "scRNAseq_Results/"
multinichenet_output = list(
  celltype_info = abundance_expression_info$celltype_info,
  celltype_de = celltype_de,
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de =  sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  prioritization_tables = prioritization_tables,
  grouping_tbl = grouping_tbl,
  lr_target_prior_cor = lr_target_prior_cor) 
multinichenet_output = multinichenetr::make_lite_output(multinichenet_output)

save = TRUE
if(save == TRUE){
  saveRDS(multinichenet_output, paste0(Path, "multinichenet_output.rds"))
  
}

multinichenet_output <- readRDS("scRNAseq_Results/multinichenet_output.rds")

# Visualization of differential cell-cell interactions
prioritized_tbl_oi_all = multinichenetr::get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  top_n = 25, 
  rank_per_group = T)

prioritized_tbl_oi = 
  multinichenet_output$prioritization_tables$group_prioritization_tbl %>%
  filter(id %in% prioritized_tbl_oi_all$id) %>%
  distinct(id, sender, receiver, ligand, receptor, group) %>% 
  left_join(prioritized_tbl_oi_all)
prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] = 0


senders_receivers = union(prioritized_tbl_oi$sender %>% unique(), prioritized_tbl_oi$receiver %>% unique()) %>% sort()

colors_sender = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)
colors_receiver = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)

circos_list = multinichenetr::make_circos_group_comparison(prioritized_tbl_oi, colors_sender, colors_receiver)



# Save the output of Circos_list
Supp_Fig4a1_circos <- circos_list$HIV.
Supp_Fig4a3_circos <- circos_list$HIV..ART.1.Year
Supp_Fig4a2_circos <- circos_list$HIV..ART.3.Months


recordedplot_to_grob <- function(recorded_plot) {
  gridGraphics::grid.echo(recorded_plot)
  grid::grid.grab()
}

Supp_Fig4a1_circos <- recordedplot_to_grob(Supp_Fig4a1_circos)
Supp_Fig4a2_circos <- recordedplot_to_grob(Supp_Fig4a2_circos)
Supp_Fig4a3_circos <- recordedplot_to_grob(Supp_Fig4a3_circos)


# Example: Convert grob to ggplot
ggplot_from_grob <- function(grob) {
  ggplot() +
    annotation_custom(grob) +
    theme_void()  # Remove axes, grid lines, etc.
}

# Convert your grob
Supp_Fig4a1_circos <- ggplot_from_grob(Supp_Fig4a1_circos)
Supp_Fig4a2_circos <- ggplot_from_grob(Supp_Fig4a2_circos)
Supp_Fig4a3_circos <- ggplot_from_grob(Supp_Fig4a3_circos)

# Combine Multiple ggplot Objects
Supp_Fig4a <- Supp_Fig4a1_circos + Supp_Fig4a2_circos + Supp_Fig4a3_circos + plot_layout(ncol = 3)
print(Supp_Fig4a)



#################################### FIGURE 5 ##################################
# Differential gene expression using MAST

# Load libraries
library(Seurat)
library(dplyr)



Figure_3def_5abcd_data <- readRDS("HIV-PAPER/Data/SPN-PLHIV_2025_data/Figure_3def_5abcd_data.rds")
Immune_cells <- Figure_5abcd_data
# Immune cells
# Function to run FindMarkers for a single cluster comparing HIV+ ART<3 Months with HIV-
run_seurat_mast_de3mvsneg <- function(seurat_obj, cluster_name) {
  # Subset to the cluster
  cluster_obj <- subset(seurat_obj, subset = Clusters == cluster_name)
  
  # Check cell count
  if (ncol(cluster_obj) < 20) {
    warning(paste("Cluster", cluster_name, "has fewer than 20 cells. Skipping."))
    return(NULL)
  }
  
  # Check HIV_status distribution
  cell_counts <- table(cluster_obj@meta.data$HIV_Status)
  print(paste("Cluster:", cluster_name))
  print(cell_counts)
  
  if (length(cell_counts) < 3 || any(cell_counts == 0)) {
    warning(paste("Cluster", cluster_name, "missing some HIV_Status levels. Skipping."))
    return(NULL)
  }
  
  # Run FindMarkers for all pairwise comparisons within the cluster
  de_results <- Seurat::FindMarkers(cluster_obj,
                                    ident.1 = c("HIV+ ART<3 Months"),
                                    ident.2 = "HIV-",
                                    group.by = "HIV_Status",
                                    test.use = "MAST",
                                    logfc.threshold = 0.1,  # Adjust as needed
                                    min.pct = 0.05)         # Filter genes expressed in 5% of cells
  
  # Add cluster info and adjust p-values
  de_results <- de_results %>%
    mutate(cluster = cluster_name,
           adj_pval = p.adjust(p_val, method = "BH")) %>%
    rownames_to_column("gene") %>%
    dplyr::select(cluster, gene, avg_log2FC, p_val, adj_pval, pct.1, pct.2)
  
  return(de_results)
}

# Get all unique clusters
clusters <- unique(Immune_cells@meta.data$Clusters)
print(paste("Found", length(clusters), "clusters:", paste(clusters, collapse = ", ")))

# Run FindMarkers for all clusters
de_results_list <- lapply(clusters, function(cluster) {
  tryCatch({
    result <- run_seurat_mast_de3mvsneg(Immune_cells, cluster)
    if (!is.null(result)) {
      write.csv(result, paste0("scRNAseq_Results/Seurat_MAST_DE_Analysis/seurat_mast_3mvsneg_de_", gsub(" ", "_", cluster), ".csv"), row.names = FALSE)
      return(result)
    }
  }, error = function(e) {
    warning(paste("Error in cluster", cluster, ":", e$message))
    return(NULL)
  })
})

# Name the list
names(de_results_list) <- clusters


# Function to run FindMarkers for a single cluster comparing HIV+ ART>1 Year with HIV-
run_seurat_mast_de1yvsneg <- function(seurat_obj, cluster_name) {
  # Subset to the cluster
  cluster_obj <- subset(seurat_obj, subset = Clusters == cluster_name)
  
  # Check cell count
  if (ncol(cluster_obj) < 20) {
    warning(paste("Cluster", cluster_name, "has fewer than 20 cells. Skipping."))
    return(NULL)
  }
  
  # Check HIV_status distribution
  cell_counts <- table(cluster_obj@meta.data$HIV_Status)
  print(paste("Cluster:", cluster_name))
  print(cell_counts)
  
  if (length(cell_counts) < 3 || any(cell_counts == 0)) {
    warning(paste("Cluster", cluster_name, "missing some HIV_Status levels. Skipping."))
    return(NULL)
  }
  
  # Run FindMarkers for all pairwise comparisons within the cluster
  de_results <- FindMarkers(cluster_obj,
                            ident.1 = c("HIV+ ART>1 Year"),
                            ident.2 = "HIV-",
                            group.by = "HIV_Status",
                            test.use = "MAST",
                            logfc.threshold = 0.1,  # Adjust as needed
                            min.pct = 0.05)         # Filter genes expressed in 5% of cells
  
  # Add cluster info and adjust p-values
  de_results <- de_results %>%
    mutate(cluster = cluster_name,
           adj_pval = p.adjust(p_val, method = "BH")) %>%
    rownames_to_column("gene") %>%
    dplyr::select(cluster, gene, avg_log2FC, p_val, adj_pval, pct.1, pct.2)
  
  return(de_results)
}

# Get all unique clusters
clusters <- unique(Immune_cells@meta.data$Clusters)
print(paste("Found", length(clusters), "clusters:", paste(clusters, collapse = ", ")))

# Run FindMarkers for all clusters
de_results_list <- lapply(clusters, function(cluster) {
  tryCatch({
    result <- run_seurat_mast_de1yvsneg(Immune_cells, cluster)
    if (!is.null(result)) {
      write.csv(result, paste0("scRNAseq_Results/Seurat_MAST_DE_Analysis/seurat_mast_1yvsneg_de_", gsub(" ", "_", cluster), ".csv"), row.names = FALSE)
      return(result)
    }
  }, error = function(e) {
    warning(paste("Error in cluster", cluster, ":", e$message))
    return(NULL)
  })
})

# Name the list
names(de_results_list) <- clusters


# Function to run FindMarkers for a single cluster comparing HIV+ ART<3 Months with HIV+ ART>1 Year
run_seurat_mast_de3mvs1y <- function(seurat_obj, cluster_name) {
  # Subset to the cluster
  cluster_obj <- subset(seurat_obj, subset = Clusters == cluster_name)
  
  # Check cell count
  if (ncol(cluster_obj) < 20) {
    warning(paste("Cluster", cluster_name, "has fewer than 20 cells. Skipping."))
    return(NULL)
  }
  
  # Check HIV_status distribution
  cell_counts <- table(cluster_obj@meta.data$HIV_Status)
  print(paste("Cluster:", cluster_name))
  print(cell_counts)
  
  if (length(cell_counts) < 3 || any(cell_counts == 0)) {
    warning(paste("Cluster", cluster_name, "missing some HIV_Status levels. Skipping."))
    return(NULL)
  }
  
  # Run FindMarkers for all pairwise comparisons within the cluster
  de_results <- FindMarkers(cluster_obj,
                            ident.1 = c("HIV+ ART<3 Months"),
                            ident.2 = "HIV+ ART>1 Year",
                            group.by = "HIV_Status",
                            test.use = "MAST",
                            logfc.threshold = 0.1,  # Adjust as needed
                            min.pct = 0.05)         # Filter genes expressed in 5% of cells
  
  # Add cluster info and adjust p-values
  de_results <- de_results %>%
    mutate(cluster = cluster_name,
           adj_pval = p.adjust(p_val, method = "BH")) %>%
    rownames_to_column("gene") %>%
    dplyr::select(cluster, gene, avg_log2FC, p_val, adj_pval, pct.1, pct.2)
  
  return(de_results)
}

# Get all unique clusters
clusters <- unique(Immune_cells@meta.data$Clusters)
print(paste("Found", length(clusters), "clusters:", paste(clusters, collapse = ", ")))

# Run FindMarkers for all clusters
de_results_list <- lapply(clusters, function(cluster) {
  tryCatch({
    result <- run_seurat_mast_de3mvs1y(Immune_cells, cluster)
    if (!is.null(result)) {
      write.csv(result, paste0("scRNAseq_Results/Seurat_MAST_DE_Analysis/seurat_mast_3mvs1y_de_", gsub(" ", "_", cluster), ".csv"), row.names = FALSE)
      return(result)
    }
  }, error = function(e) {
    warning(paste("Error in cluster", cluster, ":", e$message))
    return(NULL)
  })
})

# Name the list
names(de_results_list) <- clusters

# Volcanoplot Neutrophils 3mvsneg
seurat_mast_3mvsneg_de_Neutrophils <- read_csv("scRNAseq_Results/Seurat_MAST_DE_Analysis/seurat_mast_3mvsneg_de_Neutrophils.csv") %>%
  #dplyr::filter(pct.1>.05) %>%
  dplyr::mutate(Significance = ifelse(p_val<.05 & avg_log2FC>.25, "Upregulated",
                                      ifelse(p_val<.05 & avg_log2FC< -.25,"Downregulated","Not significant"))) %>%
  dplyr::mutate(label = ifelse(p_val<.05 & abs(avg_log2FC)>.25,gene,NA))

Fig5b <- seurat_mast_3mvsneg_de_Neutrophils %>%
  ggplot(aes(avg_log2FC,-log10(p_val),color = Significance))+
  geom_point(size=2.5,alpha=.5)+
  geom_text_repel(aes(label = label),size=8)+
  geom_hline(yintercept = 1.3,linetype='dashed')+
  geom_vline(xintercept = c(-.25,.25),linetype='dashed')+
  scale_color_manual(values = c("blue","grey","red"))+
  labs(x=expression("Average"~"log"[2]~"FoldChange"),
       y=expression("-"~"log"[10]~"adjusted P-value"),
       title = 'b.',
       subtitle = "PLHIV-ART<3m vs HIV-")+
  scale_x_continuous(limits = c(-5,3))+
  theme_classic()+
  theme(plot.title = element_text(size = 60, face = "bold",hjust = -0.1),
        legend.position = 'none',
        axis.ticks.length = unit(0.5,'cm'),
        plot.subtitle = element_text(size = 30,hjust = .5),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 35))
Fig5b

# Save Figure 5b
ggsave(Fig5b,filename="HIV-PAPER/Figures/Figure 5/Fig5b.png",
       width = 14,height = 18,dpi = 300)
ggsave(Fig5b,filename="HIV-PAPER/Figures/Figure 5/Fig5b.pdf",
       width = 14,height = 18,dpi = 300)

# Volcanoplot Neutrophils 1yvsneg
seurat_mast_1yvsneg_de_Neutrophils <- read_csv("scRNAseq_Results/Seurat_MAST_DE_Analysis/seurat_mast_1yvsneg_de_Neutrophils.csv") %>%
  #dplyr::filter(pct.1>.05) %>%
  dplyr::mutate(Significance = ifelse(p_val<.05 & avg_log2FC>.25, "Upregulated",
                                      ifelse(p_val<.05 & avg_log2FC< -.25,"Downregulated","Not significant"))) %>%
  dplyr::mutate(label = ifelse(p_val<.05 & abs(avg_log2FC)>.25,gene,NA))

Fig5a <-seurat_mast_1yvsneg_de_Neutrophils %>%
  ggplot(aes(avg_log2FC,-log10(p_val),color = Significance))+
  geom_point(size=2.5,alpha=.5)+
  #geom_text_repel(data = seurat_mast_3mvsneg_de_Neutrophils %>%
  #filter(gene %in% c("CXCL8", "ISG15","IL1B",
  #"IFITM1","NFKBIZ","HLA-F","TNFAIP6","HLA-F",
  #"HLA-A","LCN2","ADM","SERPINB3","TGM2","ICAM1")),
  #aes(label = gene),#vjust = -1, 
  #size = 3,color="black",face="bold") + 
  geom_text_repel(aes(label = label),size=8)+
  geom_hline(yintercept = 1.3,linetype='dashed')+
  geom_vline(xintercept = c(-.25,.25),linetype='dashed')+
  scale_color_manual(values = c("blue","grey","red"))+
  labs(x=expression("Average"~"log"[2]~"FoldChange"),
       y=expression("-"~"log"[10]~"adjusted P-value"),
       title = 'a.',
       subtitle = "PLHIV-ART>1yr vs HIV-")+
  scale_x_continuous(limits = c(-5,3))+
  theme_classic()+
  theme(plot.title = element_text(size = 60, face = "bold",hjust = -0.1),
        legend.position = 'none',
        axis.ticks.length = unit(0.5,'cm'),
        plot.subtitle = element_text(size = 30,hjust = .5),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 35))
Fig5a
# Save Figure 6a2
ggsave(Fig5a,filename="HIV-PAPER/Figures/Figure 5/Fig5a.png",
       width = 14,height = 18,dpi = 300)
ggsave(Fig5a,filename="HIV-PAPER/Figures/Figure 5/Fig5a.pdf",
       width = 14,height = 18,dpi = 300)


# Volcanoplot Neutrophils 3mvs1y
seurat_mast_3mvs1y_de_Neutrophils <- read_csv("scRNAseq_Results/Seurat_MAST_DE_Analysis/seurat_mast_3mvs1y_de_Neutrophils.csv") %>%
  #dplyr::filter(pct.1>.05) %>%
  dplyr::mutate(Significance = ifelse(adj_pval<.05 & avg_log2FC>.25, "Upregulated",
                                      ifelse(adj_pval<.05 & avg_log2FC< -.25,"Downregulated","Not significant"))) %>%
  dplyr::mutate(label = ifelse(adj_pval<.05 & abs(avg_log2FC)>.25,gene,NA))

Fig5c <- seurat_mast_3mvs1y_de_Neutrophils %>%
  ggplot(aes(avg_log2FC,-log10(p_val),color = Significance))+
  geom_point(size=2.5,alpha=.5)+
  geom_text_repel(aes(label = label),size=8)+
  geom_hline(yintercept = 1.3,linetype='dashed')+
  geom_vline(xintercept = c(-.25,.25),linetype='dashed')+
  scale_color_manual(values = c("blue","grey","red"))+
  labs(x=expression("Average"~"log"[2]~"FoldChange"),
       y=expression("-"~"log"[10]~"adjusted P-value"),
       title = 'c.',
       subtitle = "PLHIV-ART<3m vs PLHIV-ART>1yr")+
  scale_x_continuous(limits = c(-5,3))+
  theme_classic()+
  theme(plot.title = element_text(size = 60, face = "bold",hjust = -0.1),
        legend.position = 'none',
        axis.ticks.length = unit(0.5,'cm'),
        plot.subtitle = element_text(size = 30,hjust = .5),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 35))
Fig5c
# Save Figure 5c
ggsave(Fig5c,filename="HIV-PAPER/Figures/Figure 5/Fig5c.png",
       width = 14,height = 18,dpi = 300)
ggsave(Fig5c,filename="HIV-PAPER/Figures/Figure 5/Fig5c.pdf",
       width = 14,height = 18,dpi = 300)


# CEMITool
# CEMiTool in Neutrophils
# Aggregate expression
Neutrophils <- subset(Figure_4_data, idents = "Neutrophils")

Neutrophil_cts <- Seurat::AggregateExpression(
  object = Neutrophils,
  group.by = c("HIV_Status","sample"),
  assays = "RNA",
  slot = "counts",
  return.seurat = F)
# Convert to dataframe
Neutrophil_cts <- as.data.frame(Neutrophil_cts$RNA)
# generate sample level metadata
colData <- data.frame(samples = colnames(Neutrophil_cts))
colData <- colData %>%
  dplyr::mutate(HIV_Status = ifelse(grepl("HIV-", samples), "HIV-", 
                                    ifelse(grepl("ART<3",samples),"PLHIV on ART<3m","PLHIV on ART>1y"))) %>%
  dplyr:: mutate(Sample_Name = str_extract(samples, "(?<=_).*")) %>%
  column_to_rownames(var = "samples")
# Create a DESeq2 object
Neutrophil_dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = Neutrophil_cts,
  colData = colData,
  design = ~ HIV_Status)
# Run DESeq2
Neutrophil_dds <- DESeq2::DESeq(Neutrophil_dds)
# Extract data from the dds object
Neutrophils_out <- DESeq2::rlog(Neutrophil_dds,blind = FALSE)
# Create a matrix
Neutrophils_matrix <- assay(Neutrophil_dds)
Neutrophils_matrix <- na.omit(Neutrophils_matrix)
Neutrophils_matrix <- as.data.frame(Neutrophils_matrix)
# Create a sample annotation file
sample_annot <- colData %>%
  dplyr::mutate(SampleName = rownames(.)) %>%
  dplyr::mutate(Class = ifelse(grepl("HIV-", SampleName), "HIV-", 
                               ifelse(grepl("ART<3",SampleName),"PLHIV on ART<3m","PLHIV on ART>1y"))) %>%
  dplyr::select(SampleName,Class)
write.csv(sample_annot,"scRNAseq_Results/Sample_annotation.csv",row.names = F)
# read GMT file
gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
gmt_in <- read_gmt(gmt_fname)
# read interactions
int_fname <- system.file("extdata", "interactions.tsv", package = "CEMiTool")
int_df <- read.delim(int_fname)
head(int_df)
# Run CEMiTool
Neutrophil_cem <- CEMiTool::cemitool(Neutrophils_matrix,
                                     annot = sample_annot,
                                     gmt = gmt_in,
                                     interactions = int_df,
                                     filter = TRUE,
                                     filter_pval = 0.05,
                                     #apply_vst = T,
                                     cor_method = "pearson",
                                     cor_function = "cor",
                                     network_type = "signed",
                                     gsea_scale = T,
                                     force_beta = T,
                                     ora_pval = 0.05,
                                     min_ngen = 5,
                                     gsea_min_size = 3,
                                     gsea_max_size = 800,
                                     plot = T,
                                     center_func = 'median',
                                     plot_diagnostics = T)

Neutrophil_cem <- readRDS("scRNAseq_Results/Neutrophil_cem.rds")
view(Neutrophil_cem@enrichment)

gsea <- show_plot(Neutrophil_cem, "gsea")

es <- Neutrophil_cem@enrichment$es %>%
  pivot_longer(cols = c("HIV-","PLHIV on ART<3m","PLHIV on ART>1y"),
               names_to = "HIV_Status",
               values_to = "ES")

Nes <- Neutrophil_cem@enrichment$nes %>%
  pivot_longer(cols = c("HIV-","PLHIV on ART<3m","PLHIV on ART>1y"),
               names_to = "HIV_Status",
               values_to = "NES")

Adjusted_Pvalue <- Neutrophil_cem@enrichment$padj %>%
  pivot_longer(cols = c("HIV-","PLHIV on ART<3m","PLHIV on ART>1y"),
               names_to = "HIV_Status",
               values_to = "adj_pvalue")

Neutrophil_enrichment <- Nes %>%
  merge(es, by=c('pathway','HIV_Status'), all=T) %>%
  merge(Adjusted_Pvalue, by=c('pathway','HIV_Status'), all=T)

Fig5d <- Neutrophil_enrichment %>%
  ggplot(aes(x=HIV_Status,y=pathway#,fill = NES
  ))+
  geom_point(aes(col = NES,
                 size = -log10(adj_pvalue)))+
  scale_size(range = c(0, 40),guide = "none")+
  scale_color_gradient2(low = "#0000B4",mid="white", high = "#931500")+
  scale_fill_gradient2(low = "#0000B4",mid="white", high = "#931500")+
  labs(subtitle = "Neutrophil enriched modules",
       title = "d.")+
  scale_x_discrete(labels=c("HIV-"="HIV-",
                            "PLHIV on ART<3m"="PLHIV\nART<3m",
                            "PLHIV on ART>1y"="PLHIV\nART>1yr"))+
  theme_classic()+
  theme(plot.title = element_text(size = 60, face = "bold",hjust = -0.1),
        #legend.position = 'none',
        axis.title = element_blank(),
        axis.ticks.length = unit(0.5,'cm'),
        plot.subtitle = element_text(size = 30, hjust = .5),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30))+
  guides(
    color = guide_colorbar(
      title = "NES",
      title.theme = element_text(size = 30),
      label.theme = element_text(size = 30),
      barwidth = 2.7,
      barsize = 4,
      barheight=20),
    size = guide_legend(
      title = "-log"[10]~"adj p-value",
      title.theme = element_text(size = 30),        # Size of the legend title
      label.theme = element_text(size = 30),        # Size of the legend labels
      override.aes = list(size = c(15,25,40))))
Fig5d
ggsave(Fig5d,filename="HIV-PAPER/Figures/Figure 5/Fig5d.png",
       width = 10,height = 12, dpi = 300,units = "in")
ggsave(Fig5d,filename="HIV-PAPER/Figures/Figure 5/Fig5d.pdf",
       width = 10,height = 12, dpi = 300,units = "in")


Neutrophil_M1_pathways <- Neutrophil_cem@ora %>%
  dplyr::filter(Module=="M1",p.adjust<0.05)

Neutrophil_cem@ora$ID <- str_wrap(Neutrophil_cem@ora$ID, width = 40)
id_colors <- c("black","black","black","black","black","black","black","black","black","black")
Fig5e <- Neutrophil_cem@ora %>%
  dplyr::filter(Module=="M1", p.adjust<0.2) %>%
  sort("p.adjust") %>%
  slice_head(n=10) %>%
  ggplot(aes(-log10(p.adjust),fct_reorder(ID,-log10(p.adjust)), #fill = -log10(p.adjust),
             fill=FoldEnrichment))+
  geom_col()+
  scale_fill_gradient(low = "#4A527F", high = "#931500")+
  #scale_fill_gradient(low = "grey", high = "darkblue")+
  labs(x = "-log10(padj)", 
       title = "e.",
       subtitle = "M1 over-represented pathways")+
  theme_classic2()+
  theme(panel.grid = element_blank(),
        legend.position = c(0.70,0.25),
        axis.title.y = element_blank(),
        axis.ticks.length = unit(0.5,'cm'),
        axis.title.x = element_text(size = 30),
        axis.text.y = element_text(size = 23),
        axis.text.x = element_text(size = 23),
        plot.title = element_text(size = 60, face = 'bold',hjust = -0.5),
        plot.subtitle = element_text(hjust = 0.5,size = 30))+
  guides(fill = guide_colorbar(
    title = "FoldEnrichment",
    title.theme = element_text(size = 30),
    label.theme = element_text(size = 23),
    barwidth = 1.5,
    barsize = 1.5,
    barheight=15.5))+
  geom_vline(xintercept = -log10(0.05), linetype='dashed')
Fig5e
ggsave(Fig5e,filename="HIV-PAPER/Figures/Figure 5/Fig5e.png",
       width = 15,height = 14,dpi = 300,units = "in")
ggsave(Fig5e,filename="HIV-PAPER/Figures/Figure 5/Fig5e.pdf",
       width = 15,height = 14,dpi = 300,units = "in")


Neutrophil_M2_pathways <- Neutrophil_cem@ora %>%
  dplyr::filter(p.adjust<.05,
                Module=="M2")
id_colors <- c("black","black","black","black","black","black","black","black","black","black")
Neutrophil_cem@ora$ID <- str_wrap(Neutrophil_cem@ora$ID, width = 40)
Fig5f <- Neutrophil_cem@ora %>%
  dplyr::filter(Module=="M2", p.adjust<0.1) %>%
  sort("p.adjust") %>%
  slice_head(n=10) %>%
  ggplot(aes(-log10(p.adjust),fct_reorder(ID,-log10(p.adjust)), #fill = -log10(p.adjust),
             fill=FoldEnrichment))+
  geom_col()+
  scale_fill_gradient(low = "#4A527F", high = "#931500")+
  #scale_fill_gradient(low = "grey", high = "darkblue")+
  labs(x = "-log10(padj)", 
       title = "f.",
       subtitle = "M2 over-represented pathways")+
  theme_classic2()+
  theme(panel.grid = element_blank(),
        legend.position = c(0.73,0.25),
        axis.title.y = element_blank(),
        axis.ticks.length = unit(0.5,'cm'),
        axis.title.x = element_text(size = 30),
        axis.text.y = element_text(size = 23),
        axis.text.x = element_text(size = 23),
        plot.title = element_text(size = 60, face = 'bold',hjust = -0.5),
        plot.subtitle = element_text(hjust = 0.5,size = 30))+
  guides(fill = guide_colorbar(
    title = "FoldEnrichment",
    title.theme = element_text(size = 30),
    label.theme = element_text(size = 23),
    barwidth = 1.5,
    barsize = 1.5,
    barheight=15.5))+
  geom_vline(xintercept = -log10(0.05), linetype='dashed')
Fig5f

ggsave(Fig5f,filename="HIV-PAPER/Figures/Figure 5/Fig5f.png",
       width = 15,height = 14,dpi = 300,units = "in")
ggsave(Fig5f,filename="HIV-PAPER/Figures/Figure 5/Fig5f.pdf",
       width = 15,height = 14,dpi = 300,units = "in")

# Gene Module scores
# Function to extract gene list from CEMiTool ORA results for a given pathway
extract_gene_list <- function(cem_object, pathway_description) {
  cem_object@ora %>%
    dplyr::filter(Description == pathway_description) %>%
    dplyr::pull(geneID) %>%
    strsplit(split = "/") %>%
    unlist()
}

M1_Neut_pathways <- Neutrophil_cem@ora %>%
  dplyr::filter(Module == "M1", p.adjust<0.1)
M2_Neut_pathways <- Neutrophil_cem@ora %>%
  dplyr::filter(Module == "M2", p.adjust<0.1)

# Define gene lists using the function (M1)
Hemostasis <- extract_gene_list(Neutrophil_cem, "Hemostasis")
SASP <- extract_gene_list(Neutrophil_cem, "Senescence-Associated Secretory Phenotype (SASP)")
# Define gene lists using the function (M2)
Cytokine_signalling <- extract_gene_list(Neutrophil_cem, "Cytokine Signaling in Immune system")
IFN_signalling <- extract_gene_list(Neutrophil_cem, "Interferon Signaling")
IFNy_signalling <- extract_gene_list(Neutrophil_cem, "Interferon gamma signaling")
IFNab_signalling <- extract_gene_list(Neutrophil_cem, "Interferon alpha/beta signaling")
APC_MHC1 <- extract_gene_list(Neutrophil_cem, "Antigen Presentation: Folding, assembly and peptide loading of class I MHC")
Immunoregulatory_interactions <- extract_gene_list(Neutrophil_cem, "Immunoregulatory interactions between a Lymphoid and a non-Lymphoid cell")
APC_Cross_presentation <- extract_gene_list(Neutrophil_cem, "Antigen processing-Cross presentation")
TLR4_signalling <- extract_gene_list(Neutrophil_cem, "Activated TLR4 signalling")
RIG_MDA5 <- extract_gene_list(Neutrophil_cem, "RIG-I/MDA5 mediated induction of IFN-alpha/beta pathways")
MyD88 <- extract_gene_list(Neutrophil_cem, "MyD88:Mal cascade initiated on plasma membrane")
NFKB <- extract_gene_list(Neutrophil_cem, "NF-kB is activated and signals survival")
Nef_downregulation_MHC1 <- extract_gene_list(Neutrophil_cem, "Nef mediated downregulation of MHC class I complex cell surface expression")
NOTCH1 <- extract_gene_list(Neutrophil_cem, "NOTCH1 Intracellular Domain Regulates Transcription")
ER_phagosome <- extract_gene_list(Neutrophil_cem, "ER-Phagosome pathway")
NOTCH2 <- extract_gene_list(Neutrophil_cem, "NOTCH2 intracellular domain regulates transcription")
Nef_HIV_replication <- extract_gene_list(Neutrophil_cem, "The role of Nef in HIV-1 replication and disease pathogenesis")
Downregulation_TGFbeta <- extract_gene_list(Neutrophil_cem, "Downregulation of TGF-beta receptor signaling")
TLR2 <- extract_gene_list(Neutrophil_cem, "Toll Like Receptor 2 (TLR2) Cascade")


# Create a named list of all gene lists
module_features <- list(
  Hemostasis = Hemostasis,
  SASP = SASP,
  Cytokine_signalling = Cytokine_signalling,
  IFN_signalling = IFN_signalling,
  IFNy_signalling = IFNy_signalling,
  IFNab_signalling = IFNab_signalling,
  APC_MHC1 = APC_MHC1,
  Immunoregulatory_interactions = Immunoregulatory_interactions,
  APC_Cross_presentation = APC_Cross_presentation,
  TLR4_signalling = TLR4_signalling,
  RIG_MDA5 = RIG_MDA5,
  MyD88 = MyD88,
  NFKB = NFKB,
  Nef_downregulation_MHC1 = Nef_downregulation_MHC1,
  NOTCH1 = NOTCH1,
  ER_phagosome = ER_phagosome,
  Nef_HIV_replication = Nef_HIV_replication,
  Downregulation_TGFbeta = Downregulation_TGFbeta,
  TLR2 = TLR2,
  NOTCH2 = NOTCH2
)

load("Data/Single_Cell_Data/Immune_cells.RData")
Figure_5g_to_5j_data <- Immune_cells
saveRDS(Figure_5g_to_5j_data, file="HIV-PAPER/Data/SPN-PLHIV_2025_Data/Figure_5g_to_5j_data.rds")
# Add all module scores in one call
Figure_5g_to_5j_data <- AddModuleScore(
  object = Figure_5g_to_5j_data,
  features = module_features,
  name = names(module_features),  # Names will be appended with numbers (e.g., "Glycan_biosynthesis1")
  slot = "data"
)


Fig5g <- Figure_5g_to_5j_data@meta.data %>%
  dplyr::filter(Clusters=="Neutrophils",
                SASP2>0) %>%
  ggplot(aes(HIV_Status,SASP2,fill = HIV_Status))+
  geom_violin(trim = T, scale = "width")+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=43,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 't.test',
           label = "{ifelse(p < 0.0001, 'p < 0.0001', sprintf('p = %.4f', p))}",
           label.size = 10,
           tip.length = 0.01,
           step.increase = 0.2,
           hide.ns = F)+
  labs(x='',
       title = "g.",
       subtitle = "SASP",
       y='Module Expression')+
  scale_fill_manual(values = c('#A9A9A9','#941100','#005493'))+
  scale_x_discrete(labels=c("HIV+ ART<3 Months"="PLHIV\nART<3m",
                            "HIV+ ART>1 Year"="PLHIV\nART>1yr",
                            "HIV-"="HIV-"))+
  scale_y_continuous(limits = c(-0.5,4.5))+
  theme_classic2()+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 30),
        axis.ticks.length = unit(0.5,'cm'),
        axis.title.x = element_text(size = 30),
        axis.text.y = element_text(size = 23),
        axis.text.x = element_text(size = 30),
        plot.title = element_text(size = 60, face = 'bold',hjust = -0.1),
        plot.subtitle = element_text(hjust = 0.5,size = 30))

Fig5g
ggsave(Fig5g,filename="HIV-PAPER/Figures/Figure 5/Fig5g.png",
       width = 8,height = 10,dpi = 300,units = "in")
ggsave(Fig5g,filename="HIV-PAPER/Figures/Figure 5/Fig5g.pdf",
       width = 8,height = 10,dpi = 300,units = "in")


Fig5h <- Figure_5g_to_5j_data@meta.data %>%
  dplyr::filter(Clusters=="Neutrophils",
                IFN_signalling4>0) %>%
  ggplot(aes(HIV_Status,IFN_signalling4,fill = HIV_Status))+
  geom_violin(trim = T, scale = "width")+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=43,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 'wilcox.test',
           label = "{ifelse(p < 0.0001, 'p < 0.0001', sprintf('p = %.4f', p))}",
           label.size = 10,
           tip.length = 0.01,
           step.increase = 0.2,
           p.adjust.method = "fdr",
           p.adjust.by = "group",
           hide.ns = F)+
  labs(x='',
       title = "h.",
       subtitle = "IFN signaling",
       y='Module Expression')+
  scale_fill_manual(values = c('#A9A9A9','#941100','#005493'))+
  scale_x_discrete(labels=c("HIV+ ART<3 Months"="PLHIV\nART<3m",
                            "HIV+ ART>1 Year"="PLHIV\nART>1yr",
                            "HIV-"="HIV-"))+
  scale_y_continuous(limits = c(0,4.5))+
  theme_classic2()+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 30),
        axis.ticks.length = unit(0.5,'cm'),
        axis.title.x = element_text(size = 30),
        axis.text.y = element_text(size = 23),
        axis.text.x = element_text(size = 30),
        plot.title = element_text(size = 60, face = 'bold',hjust = -0.1),
        plot.subtitle = element_text(hjust = 0.5,size = 30))

Fig5h
ggsave(Fig5h,filename="HIV-PAPER/Figures/Figure 5/Fig5h.png",
       width = 8,height = 10,dpi = 300,units = "in")
ggsave(Fig5h,filename="HIV-PAPER/Figures/Figure 5/Fig5h.pdf",
       width = 8,height = 10,dpi = 300,units = "in")


Fig5i <- Figure_5g_to_5j_data@meta.data %>%
  dplyr::filter(Clusters=="Neutrophils",
                Nef_HIV_replication17>0) %>%
  ggplot(aes(HIV_Status,IFN_signalling4,fill = HIV_Status))+
  geom_violin(trim = T, scale = "width")+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=43,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 't.test',
           label = "{ifelse(p < 0.0001, 'p < 0.0001', sprintf('p = %.4f', p))}",
           label.size = 10,
           tip.length = 0.01,
           step.increase = 0.2,
           hide.ns = F)+
  labs(x='',
       title = "i.",
       subtitle = "Nef HIV replication",
       y='Module Expression')+
  scale_fill_manual(values = c('#A9A9A9','#941100','#005493'))+
  scale_x_discrete(labels=c("HIV+ ART<3 Months"="PLHIV\nART<3m",
                            "HIV+ ART>1 Year"="PLHIV\nART>1yr",
                            "HIV-"="HIV-"))+
  scale_y_continuous(limits = c(0,4.5))+
  theme_classic2()+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 30),
        axis.ticks.length = unit(0.5,'cm'),
        axis.title.x = element_text(size = 30),
        axis.text.y = element_text(size = 23),
        axis.text.x = element_text(size = 30),
        plot.title = element_text(size = 60, face = 'bold',hjust = -0.1),
        plot.subtitle = element_text(hjust = 0.5,size = 30))

Fig5i
ggsave(Fig5i,filename="HIV-PAPER/Figures/Figure 5/Fig5i.png",
       width = 8,height = 10,dpi = 300,units = "in")
ggsave(Fig5i,filename="HIV-PAPER/Figures/Figure 5/Fig5i.pdf",
       width = 8,height = 10,dpi = 300,units = "in")

# Figure 6d
Fig5j <- Figure_5g_to_5j_data@meta.data %>%
  dplyr::filter(Clusters=="Neutrophils",
                SASP2>0) %>%
  ggplot(aes(HIV_Status,APC_MHC17,fill = HIV_Status))+
  geom_violin(trim = T, scale = "width")+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=43,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 't.test',
           label = "{ifelse(p < 0.0001, 'p < 0.0001', sprintf('p = %.4f', p))}",
           label.size = 10,
           tip.length = 0.01,
           step.increase = 0.2,
           hide.ns = F)+
  labs(x='',
       title = "j.",
       subtitle = "MHC1 APC",
       y='Module Expression')+
  scale_fill_manual(values = c('#A9A9A9','#941100','#005493'))+
  scale_x_discrete(labels=c("HIV+ ART<3 Months"="PLHIV\nART<3m",
                            "HIV+ ART>1 Year"="PLHIV\nART>1yr",
                            "HIV-"="HIV-"))+
  #scale_y_continuous(limits = c(-0.5,4.5))+
  theme_classic2()+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 30),
        axis.ticks.length = unit(0.5,'cm'),
        axis.title.x = element_text(size = 30),
        axis.text.y = element_text(size = 23),
        axis.text.x = element_text(size = 30),
        plot.title = element_text(size = 60, face = 'bold',hjust = -0.1),
        plot.subtitle = element_text(hjust = 0.5,size = 30))

Fig5j
ggsave(Fig5j,filename="HIV-PAPER/Figures/Figure 5/Fig5j.png",
       width = 8,height = 10,dpi = 300,units = "in")
ggsave(Fig5j,filename="HIV-PAPER/Figures/Figure 5/Fig5j.pdf",
       width = 8,height = 10,dpi = 300,units = "in")


Phagocytosis_data <- read_csv("Data/Main_Files_Thesis/Neutrophil_Phagocytosis.csv")
Figure_5k_5l_data <- read_csv("HIV-PAPER/Data/SPN-PLHIV_2025_data/Figure_5k_5l_data.csv")
# MFI of Phagocytosis
Fig5k_counts <- Figure_5k_5l_data %>%
  dplyr::filter(Time=="45min",
                #Visit=="Week 1",
                Stimulant=="LPS") %>%
  dplyr::mutate(Oxidation_SI = (`MFI+`-`MFI-`)/(2*`rSD-`)) %>%
  dplyr::mutate(Product_Phagocytosis_MF1_Frequency = `MFI phagocytosis`*`Proportion of Phagocytosis`) %>%
  dplyr::select(Product_Phagocytosis_MF1_Frequency,`HIV Status`) %>%
  dplyr::group_by(`HIV Status`) %>%
  dplyr::summarize(n = dplyr::n(), .groups = "drop")
Fig5k_counts


Fig5k <- Figure_5k_5l_data %>%
  dplyr::filter(Time=="45min",
                #Visit=="Week 1",
                Stimulant=="LPS") %>%
  dplyr::mutate(Oxidation_SI = (`MFI+`-`MFI-`)/(2*`rSD-`)) %>%
  dplyr::mutate(Product_Phagocytosis_MF1_Frequency = `MFI phagocytosis`*`Proportion of Phagocytosis`) %>%
  ggplot(aes(x=factor(`HIV Status`,
                      levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y')),
             y=log10(Product_Phagocytosis_MF1_Frequency),
             color=factor(`HIV Status`,
                          levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y'))))+
  geom_boxplot(width=0.5,notch = F, outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.2),size=7) +
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=7,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 'wilcox.test',
           label = "{ifelse(p < 0.0001, 'p < 0.0001', sprintf('p = %.4f', p))}",
           label.size = 10,
           tip.length = 0.01,
           hide.ns = F)+
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  labs(x='',
       y="Phagocytosis Index (Log10)",
       #fill='HIV Status',
       title = "k.")+
  scale_x_discrete(labels=c("HIV-"="HIV-",
                            "PLHIV ART <3m"="PLHIV\nART<3m",
                            "PLHIV ART >1y"="PLHIV\nART>1yr"))+
  #scale_y_continuous(breaks = c(0,5,10,15,20,25),limits = c(0,25))+
  theme_classic()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        plot.title = element_text(size = 60, face = 'bold',hjust = -0.1),
        plot.subtitle = element_text(size = 30, face = "plain",hjust = .5),
        axis.ticks.length = unit(0.5,"cm"),
        legend.text = element_text(size=25),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 30))
Fig5k

ggsave(Fig5k,filename="HIV-PAPER/Figures/Figure 5/Fig5k.png",
       width = 8,height = 10,dpi = 300,units = "in")
ggsave(Fig5k,filename="HIV-PAPER/Figures/Figure 5/Fig5k.pdf",
       width = 8,height = 10,dpi = 300,units = "in")


# MFI of Oxidation
Fig5l_counts <- Figure_5k_5l_data %>%
  dplyr::filter(Time=="45min",
                #Visit=="Week 1",
                Stimulant=="LPS") %>%
  dplyr::mutate(Oxidation_SI = (`MFI+`-`MFI-`)/(2*`rSD-`)) %>%
  dplyr::mutate(Product_Phagocytosis_MF1_Frequency = `MFI phagocytosis`*`Proportion of Phagocytosis`) %>%
  dplyr::select(Oxidation_SI,`HIV Status`) %>%
  dplyr::group_by(`HIV Status`) %>%
  dplyr::summarize(n = dplyr::n(), .groups = "drop")
Fig5l_counts

Fig5l <- Figure_5k_5l_data %>%
  dplyr::filter(Time=="45min",
                #Visit=="Week 1",
                Stimulant=="LPS") %>%
  dplyr::mutate(Oxidation_SI = (`MFI+`-`MFI-`)/(2*`rSD-`)) %>%
  dplyr::mutate(Product_Oxidation_MF1_Frequency = Oxidation_SI*`Reporter out of Phagocytosed`) %>%
  dplyr::filter(log10(Product_Oxidation_MF1_Frequency)<3) %>%
  ggplot(aes(x=factor(`HIV Status`,
                      levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y')),
             y=log10(Product_Oxidation_MF1_Frequency),
             color=factor(`HIV Status`,
                          levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y'))))+
  geom_boxplot(width=0.5,notch = F, outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.2),size=7) +
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=7,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 'wilcox.test',
           label = "{ifelse(p < 0.0001, 'p < 0.0001', sprintf('p = %.4f', p))}",
           label.size = 10,
           tip.length = 0.01,
           hide.ns = F)+
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  labs(x='',
       y="Oxidation Index [Log10]",
       #fill='HIV Status',
       title = "l."#,
       #subtitle = "Frequency of Oxidation * Oxidation MFI"
  )+
  scale_x_discrete(labels=c("HIV-"="HIV-",
                            "PLHIV ART <3m"="PLHIV\nART<3m",
                            "PLHIV ART >1y"="PLHIV\nART>1yr"))+
  #scale_y_continuous(breaks = c(0,5,10,15,20,25),limits = c(0,25))+
  theme_classic()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        plot.title = element_text(size = 60, face = 'bold',hjust = -0.1),
        plot.subtitle = element_text(size = 30, face = "plain",hjust = .5),
        axis.ticks.length = unit(0.5,"cm"),
        legend.text = element_text(size=25),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 30))
Fig5l

ggsave(Fig5l,filename="HIV-PAPER/Figures/Figure 5/Fig5l.png",
       width = 8,height = 10,dpi = 300,units = "in")
ggsave(Fig5l,filename="HIV-PAPER/Figures/Figure 5/Fig5l.pdf",
       width = 8,height = 10,dpi = 300,units = "in")


# Save figure 5 ggsave
ggsave(filename = "HIV-PAPER/Figures/Figure 5/Fig5.png",
       plot = ((Fig5a|Fig5b|Fig5c|plot_layout(ncol = 3,nrow = 1, width = c(1,1,1))))/ 
         (Fig5d|Fig5e|Fig5f|plot_layout(width = c(0.8,1,1)))/
         (Fig5g|Fig5h|Fig5i|Fig5j|plot_layout(width = c(1,1,1,1)))/
         (Fig5k|Fig5l|plot_spacer()|plot_layout(width = c(1,1,1,1))),
       width = 42, height = 50, units = "in", dpi = 300,limitsize = F)

ggsave(filename = "HIV-PAPER/Figures/Figure 5/Fig5.pdf",
       plot = ((Fig5a|Fig5b|Fig5c|plot_layout(ncol = 3,nrow = 1, width = c(1,1,1))))/ 
         (Fig5d|Fig5e|Fig5f|plot_layout(width = c(0.8,1,1)))/
         (Fig5g|Fig5h|Fig5i|Fig5j|plot_layout(width = c(1,1,1,1)))/
         (Fig5k|Fig5l|plot_spacer()|plot_layout(width = c(1,1,1,1))),
       width = 42, height = 50, units = "in", dpi = 300,limitsize = F)

#################################### SUPPLEMENTARY FIGURE 5 #########################
# Save figure 5 ggsave
Extended_Fig5a <- Neutrophil_cem@ora %>%
  dplyr::filter(Module=="M1") %>%
  dplyr::mutate(Significance=ifelse(p.adjust<.05 & FoldEnrichment>0,"Upregulated",
                                    ifelse(p.adjust<.05 & FoldEnrichment<0,"Downregulated","Not significant")),
                label = ifelse(Significance=="Upregulated",Description,NA)) %>%
  ggplot(aes(FoldEnrichment,-log10(p.adjust), color=Significance))+
  geom_point(size=2.5,alpha=.5)+
  geom_hline(yintercept = -log10(0.05), linetype="dashed")+
  geom_text_repel(aes(label = label),size=8,color='black')+
  scale_color_manual(values = c('grey','red'))+
  labs(x="Fold-Enrichment",
       y="-Log"[10]~"(padj)",
       title = 'a.',
       subtitle = "M1 Neutrophil enriched modules\nPLHIV-ART>1yr")+
  theme_classic()+
  theme(plot.title = element_text(size = 50, face = "bold"),
        legend.position = 'none',
        axis.ticks.length = unit(0.5,'cm'),
        plot.subtitle = element_text(size = 30, hjust = .5),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 35))
Extended_Fig5a

ggsave(Extended_Fig5a,filename="HIV-PAPER/Figures/Extended_Figure 5/Extended_Fig5a.png",
       width = 12,height = 15,dpi = 300,units = "in")
ggsave(Extended_Fig5a,filename="HIV-PAPER/Figures/Extended_Figure 5/Extended_Fig5a.pdf",
       width = 12,height = 15,dpi = 300,units = "in")

Extended_Fig5b <- Neutrophil_cem@ora %>%
  dplyr::filter(Module=="M2") %>%
  dplyr::mutate(Significance=ifelse(p.adjust<.05 & FoldEnrichment>0,"Upregulated",
                                    ifelse(p.adjust<.05 & FoldEnrichment<0,"Downregulated","Not significant")),
                label = ifelse(Significance=="Upregulated",Description,NA)) %>%
  ggplot(aes(FoldEnrichment,-log10(p.adjust), color=Significance))+
  geom_point(size=2.5,alpha=.5)+
  geom_hline(yintercept = 1.3, linetype="dashed")+
  geom_text_repel(aes(label = label),size=8,color='black')+
  scale_color_manual(values = c('grey','red'))+
  labs(x="Fold-Enrichment",
       y="-Log"[10]~"(padj)",
       title = 'b.',
       subtitle = "M2 Neutrophil enriched modules\n(PLHIV-ART<3m)")+
  theme_classic()+
  theme(plot.title = element_text(size = 50, face = "bold"),
        legend.position = 'none',
        axis.ticks.length = unit(0.5,'cm'),
        plot.subtitle = element_text(size = 30,hjust = .5),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 35))
Extended_Fig5b
# Save Figure 5b
ggsave(Extended_Fig5b,filename="HIV-PAPER/Figures/Extended_Figure 5/Extended_Fig5b.png",
       width = 12,height = 15,dpi = 300,units = "in")
ggsave(Extended_Fig5b,filename="HIV-PAPER/Figures/Extended_Figure 5/Extended_Fig5b.pdf",
       width = 12,height = 15,dpi = 300,units = "in")

Extended_Fig5c <- Seurat::DotPlot(
  Immune_cells,
  idents = "Neutrophils",
  dot.min = 0,
  dot.scale = 10,
  features = c(
    "SASP2",
    "Hemostasis1",
    "Cytokine_signalling3",
    "IFN_signalling4",
    "IFNy_signalling5",
    "IFNab_signalling6",
    "APC_MHC17",
    "Immunoregulatory_interactions8",
    "APC_Cross_presentation9",
    "TLR4_signalling10",
    "RIG_MDA511",
    "MyD8812",
    "NFKB13",
    "Nef_downregulation_MHC114",
    "NOTCH115",
    "ER_phagosome16",
    "Nef_HIV_replication17",
    "Downregulation_TGFbeta18",
    "TLR219"
  ),
  group.by = "HIV_Status")+
  coord_flip()+
  scale_y_discrete(labels=c("HIV-"="HIV-",
                            "HIV+ ART<3 Months"="PLHIV\nART<3m",
                            "HIV+ ART>1 Year"="PLHIV\nART>1yr"))+
  scale_x_discrete(labels=c("SASP2"="SASP",
                            "Hemostasis1"="Hemostasis",
                            "Cytokine_signalling3"="Cytokine signalling",
                            "IFN_signalling4"="IFN signalling",
                            "IFNy_signalling5"="IFNy signalling",
                            "IFNab_signalling6"="IFN-"~alpha~beta~ "signalling",
                            "APC_MHC17"="APC MHC-1",
                            "Immunoregulatory_interactions8"="Immunoregulatory interactions",
                            "APC_Cross_presentation9"="APC & Cross presentation",
                            "TLR4_signalling10"="TLR4 signalling",
                            "RIG_MDA511"="RIG-MDA5",
                            "MyD8812"="MyD88 signalling",
                            "NFKB13"="NF"~kappa~"B signalling",
                            "Nef_downregulation_MHC114"="Nef downregulation of MHC-1",
                            "NOTCH115"="NOTCH-1 signalling",
                            "ER_phagosome16"="ER-phagosome",
                            "Nef_HIV_replication17"="Nef in HIV replication",
                            "Downregulation_TGFbeta18"="Downregulation of"~"TGF-"~beta,
                            "TLR219"="TLR2 signalling"))+
  labs(title = "c.",
       subtitle = "Gene Module Scores")+
  theme_classic2()+
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.ticks.length = unit(0.5,'cm'),
        axis.text.y = element_text(size = 23),
        axis.text.x = element_text(size = 23),
        plot.title = element_text(size = 50, hjust = -0.5, face = 'bold'),
        plot.subtitle = element_text(hjust = 0.5,size = 30))+
  guides(
    size = guide_legend(
      title = "Percent Expressed",
      title.theme = element_text(size = 30),
      label.theme = element_text(size = 18),
      keyheight = unit(1, "cm")),
    color = guide_colorbar(
      title = "Average Expression",
      title.theme = element_text(size = 30),
      label.theme = element_text(size = 18),
      barwidth = 1))
Extended_Fig5c

ggsave(Extended_Fig5c,filename="HIV-PAPER/Figures/Extended_Figure 5/Extended_Fig5c.png",
       width = 14,height = 14,dpi = 300,units = "in")
ggsave(Extended_Fig5c,filename="HIV-PAPER/Figures/Extended_Figure 5/Extended_Fig5c.pdf",
       width = 14,height = 14,dpi = 300,units = "in")


# Save extended figure 5 ggsave
ggsave(filename = "HIV-PAPER/Figures/Extended_Figure 5/Extended_Fig5.png",
       plot = ((Extended_Fig5a|Extended_Fig5b|Extended_Fig5c|plot_layout(ncol = 3,nrow = 1, width = c(1,1,1))))/
         (plot_spacer()|plot_layout(width = c(3)))/
         (plot_spacer()|plot_layout(width = c(3))),
       width = 32, height = 29, units = "in", dpi = 300,limitsize = F)

ggsave(filename = "HIV-PAPER/Figures/Extended_Figure 5/Extended_Fig5.pdf",
       plot = ((Extended_Fig5a|Extended_Fig5b|Extended_Fig5c|plot_layout(ncol = 3,nrow = 1, width = c(1,1,1))))/
         (plot_spacer()|plot_layout(width = c(3)))/
         (plot_spacer()|plot_layout(width = c(3))),
       width = 32, height = 29, units = "in", dpi = 300,limitsize = F)


#################################### FIGURE 6 ##################################
# Differential expression in T cells
Figure_6abc_data <- readRDS("HIV-PAPER/Data/SPN-PLHIV_2025_data/Figure_6abc_data.rds")
Tcells <- Figure_6abc_data
Idents(Tcells) <- Tcells$HIV_Status
Tcell_3m_vs_Neg <- SeuratExtend::VolcanoPlot(
  Tcells,
  ident.1 = "HIV+ ART<3 Months",
  ident.2 = "HIV-"
)
Sig_Tcell_3m_vs_Neg <- Tcell_3m_vs_Neg$data

Fig6b <- Sig_Tcell_3m_vs_Neg %>%
  as.data.frame() %>%
  dplyr::mutate(Significance = ifelse(p>-log10(0.05) & logFC>.25, "Upregulated",
                                      ifelse(p>-log10(0.05) & logFC< -.25,"Downregulated","Not significant"))) %>%
  dplyr::mutate(label = ifelse(p>-log10(0.05) & abs(logFC)>.25,rownames(.),NA)) %>%
  ggplot(aes(logFC,p,color=Significance))+
  geom_point(size=2.5,alpha=.5)+
  geom_text_repel(aes(label = label),size=8)+
  geom_hline(yintercept = 1.3,linetype='dashed')+
  geom_vline(xintercept = c(-.25,.25),linetype='dashed')+
  scale_color_manual(values = c("blue","grey","red"))+
  labs(x=expression("Average"~"log"[2]~"FoldChange"),
       y=expression("-"~"log"[10]~"adjusted P-value"),
       title = 'b.',
       subtitle = "PLHIV-ART<3m vs HIV-")+
  theme_classic()+
  theme(plot.title = element_text(size = 60, face = "bold",hjust = -0.1),
        legend.position = 'none',
        axis.ticks.length = unit(0.5,'cm'),
        plot.subtitle = element_text(size = 30,hjust = .5),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 35))

Fig6b

# Save Figure 6b
ggsave(Fig6b,filename="HIV-PAPER/Figures/Figure 6/Fig6b.png",
       width = 14,height = 18,dpi = 300)
ggsave(Fig6b,filename="HIV-PAPER/Figures/Figure 6/Fig6b.pdf",
       width = 14,height = 18,dpi = 300)



Tcell_1y_vs_Neg <- SeuratExtend::VolcanoPlot(
  Tcells,
  ident.1 = "HIV+ ART>1 Year",
  ident.2 = "HIV-"
)
Tcell_1y_vs_Neg
Sig_Tcell_1y_vs_Neg <- Tcell_1y_vs_Neg$data

Fig6a <- Sig_Tcell_1y_vs_Neg %>%
  as.data.frame() %>%
  dplyr::mutate(Significance = ifelse(p>-log10(0.05) & logFC>.25, "Upregulated",
                                      ifelse(p>-log10(0.05) & logFC< -.25,"Downregulated","Not significant"))) %>%
  dplyr::mutate(label = ifelse(p>-log10(0.05) & abs(logFC)>.25,rownames(.),NA)) %>%
  ggplot(aes(logFC,p,color=Significance))+
  geom_point(size=2.5,alpha=.5)+
  geom_text_repel(aes(label = label),size=8)+
  geom_hline(yintercept = 1.3,linetype='dashed')+
  geom_vline(xintercept = c(-.25,.25),linetype='dashed')+
  scale_color_manual(values = c("blue","grey","red"))+
  labs(x=expression("Average"~"log"[2]~"FoldChange"),
       y=expression("-"~"log"[10]~"adjusted P-value"),
       title = 'a.',
       subtitle = "PLHIV-ART>1yr vs HIV-")+
  theme_classic()+
  theme(plot.title = element_text(size = 60, face = "bold",hjust = -0.1),
        legend.position = 'none',
        axis.ticks.length = unit(0.5,'cm'),
        plot.subtitle = element_text(size = 30,hjust = .5),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 35))
Fig6a

# Save Figure 6a
ggsave(Fig6a,filename="HIV-PAPER/Figures/Figure 6/Fig6a.png",
       width = 14,height = 18,dpi = 300)
ggsave(Fig6a,filename="HIV-PAPER/Figures/Figure 6/Fig6a.pdf",
       width = 14,height = 18,dpi = 300)


Tcell_3m_vs_1y <- SeuratExtend::VolcanoPlot(
  Tcells,
  ident.1 = "HIV+ ART<3 Months",
  ident.2 = "HIV+ ART>1 Year"
)
Tcell_3m_vs_1y
Sig_Tcell_3m_vs_1y <- Tcell_3m_vs_1y$dat


Fig6c <- Sig_Tcell_3m_vs_1y %>%
  as.data.frame() %>%
  dplyr::mutate(Significance = ifelse(p>-log10(0.05) & logFC>.25, "Upregulated",
                                      ifelse(p>-log10(0.05) & logFC< -.25,"Downregulated","Not significant"))) %>%
  dplyr::mutate(label = ifelse(p>-log10(0.05) & abs(logFC)>.25,rownames(.),NA)) %>%
  ggplot(aes(logFC,p,color=Significance))+
  geom_point(size=2.5,alpha=.5)+
  geom_text_repel(aes(label = label),size=8)+
  geom_hline(yintercept = 1.3,linetype='dashed')+
  geom_vline(xintercept = c(-.25,.25),linetype='dashed')+
  scale_color_manual(values = c("blue","grey","red"))+
  labs(x=expression("Average"~"log"[2]~"FoldChange"),
       y=expression("-"~"log"[10]~"adjusted P-value"),
       title = 'c.',
       subtitle = "PLHIV-ART<3m vs PLHIV-ART>1yr")+
  theme_classic()+
  theme(plot.title = element_text(size = 60, face = "bold",hjust = -0.1),
        legend.position = 'none',
        axis.ticks.length = unit(0.5,'cm'),
        plot.subtitle = element_text(size = 30,hjust = .5),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 35))
Fig6c

# Save Figure 6c
ggsave(Fig6c,filename="HIV-PAPER/Figures/Figure 6/Fig6c.png",
       width = 14,height = 18,dpi = 300)
ggsave(Fig6c,filename="HIV-PAPER/Figures/Figure 6/Fig6c.pdf",
       width = 14,height = 18,dpi = 300)

# CEMITool
Figure_6_data <- readRDS("HIV-PAPER/Data/SPN-PLHIV_2025_data/Figure_3_data.rds")
all_merged_subset_labelled_new <- Figure_6_data
all_merged_subset_labelled_new$Clusters <- paste0(all_merged_subset_labelled_new@active.ident)
Idents(all_merged_subset_labelled_new) <- all_merged_subset_labelled_new$Clusters
Tcells <- subset(all_merged_subset_labelled_new, 
                 idents = c("CD3+ T cells"),
                 invert = F)
Tcell_cts <- Seurat::AggregateExpression(
  object = Tcells,
  group.by = c("HIV_Status","sample"),
  assays = "RNA",
  slot = "counts",
  return.seurat = F)$RNA
Tcell_cts <- as.data.frame(Tcell_cts)
Keep <- rowSums(Tcell_cts) > 3
Tcell_cts <- Tcell_cts[Keep,]

colData <- data.frame(samples = colnames(Tcell_cts))
colData <- colData %>%
  dplyr::mutate(HIV_Status = ifelse(grepl("HIV-", samples), "HIV-", 
                                    ifelse(grepl("ART<3",samples),"PLHIV on ART<3m","PLHIV on ART>1y"))) %>%
  dplyr:: mutate(Sample_Name = str_extract(samples, "(?<=_).*")) %>%
  column_to_rownames(var = "samples")

sample_annot <- colData %>%
  dplyr::mutate(SampleName = rownames(.)) %>%
  dplyr::mutate(Class = ifelse(grepl("HIV-", SampleName), "HIV-", 
                               ifelse(grepl("ART<3",SampleName),"PLHIV on ART<3m","PLHIV on ART>1y"))) %>%
  dplyr::select(SampleName,Class)
write.csv(sample_annot,"scRNAseq_Results/Sample_annotation.csv",row.names = F)
# read GMT file
gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
gmt_in <- read_gmt(gmt_fname)
# read interactions
int_fname <- system.file("extdata", "interactions.tsv", package = "CEMiTool")
int_df <- read.delim(int_fname)
head(int_df)

Tcell_cem <- CEMiTool::cemitool(
  expr = Tcell_cts,
  annot = sample_annot,
  gmt = gmt_in,
  interactions = int_df,
  filter = TRUE,
  filter_pval = 0.05,
  apply_vst = TRUE,
  cor_method = "spearman",  # Spearman is often better for sparse scRNA-seq data
  cor_function = "cor",
  network_type = "signed",
  gsea_scale = TRUE,
  force_beta = TRUE,
  ora_pval = 0.05,
  min_ngen = 5,
  gsea_min_size = 3,
  gsea_max_size = 800,
  plot = TRUE,
  center_func = "median",
  plot_diagnostics = TRUE,
  verbose = TRUE
)

saveRDS(Tcell_cem, file = "scRNAseq_Results/Tcell_cem.rds")
Tcell_cem <- readRDS("scRNAseq_Results/Tcell_cem.rds")
gsea <- CEMiTool::show_plot(Tcell_cem, "gsea")
gsea
ora <- CEMiTool::show_plot(Tcell_cem, "ora")
ora$M1
ora$M2
ora$M3
ora$M4
ora$M5

M1_Tcell_pathways <- Tcell_cem@ora %>%
  dplyr::filter(p.adjust<.05,
                Module=="M1")
M2_Tcell_pathways <- Tcell_cem@ora %>%
  dplyr::filter(p.adjust<.05,
                Module=="M2")

# GSEA
es <- Tcell_cem@enrichment$es %>%
  pivot_longer(cols = c("HIV-","PLHIV on ART<3m","PLHIV on ART>1y"),
               names_to = "HIV_Status",
               values_to = "ES")

Nes <- Tcell_cem@enrichment$nes %>%
  pivot_longer(cols = c("HIV-","PLHIV on ART<3m","PLHIV on ART>1y"),
               names_to = "HIV_Status",
               values_to = "NES")

Adjusted_Pvalue <- Tcell_cem@enrichment$padj %>%
  pivot_longer(cols = c("HIV-","PLHIV on ART<3m","PLHIV on ART>1y"),
               names_to = "HIV_Status",
               values_to = "adj_pvalue")

Tcell_enrichment <- Nes %>%
  merge(es, by=c('pathway','HIV_Status'), all=T) %>%
  merge(Adjusted_Pvalue, by=c('pathway','HIV_Status'), all=T) %>%
  na.omit() %>%
  dplyr::mutate(adj_pvalue=-log10(adj_pvalue))

Module_order <- c("M1","M2","M3","M4","M5","Not.Correlated")
Fig6d <- Tcell_enrichment %>%
  dplyr::filter(pathway!="Not.Correlated") %>%
  ggplot(aes(x=HIV_Status,y=factor(pathway, levels=c("M1","M2","M3","M4","M5"))))+
  geom_point(aes(col = NES,
                 size = adj_pvalue))+
  scale_size(range = c(0, 40),guide = "none")+
  scale_color_gradient2(low = "#0000B4",mid="white", high = "#931500")+
  scale_fill_gradient2(low = "#0000B4",mid="white", high = "#931500")+
  labs(subtitle = "T cell enriched modules",
       title = "d.")+
  scale_x_discrete(labels=c("HIV-"="HIV-",
                            "PLHIV on ART<3m"="PLHIV\nART<3m",
                            "PLHIV on ART>1y"="PLHIV\nART>1yr"))+
  theme_classic()+
  theme(plot.title = element_text(size = 60, face = "bold",hjust = -0.1),
        #legend.position = 'none',
        axis.title = element_blank(),
        axis.ticks.length = unit(0.5,'cm'),
        plot.subtitle = element_text(size = 30, hjust = .5),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30))+
  guides(
    color = guide_colorbar(
      title = "NES",
      title.theme = element_text(size = 30),
      label.theme = element_text(size = 30),
      barwidth = 2.7,
      barsize = 4,
      barheight=20),
    size = guide_legend(
      title = "-log"[10]~"adj p-value",
      title.theme = element_text(size = 30),        # Size of the legend title
      label.theme = element_text(size = 30),        # Size of the legend labels
      override.aes = list(size = c(15,25,40))))
Fig6d
ggsave(Fig6d,filename="HIV-PAPER/Figures/Figure 6/Fig6d.png",
       width = 10,height = 12, dpi = 300,units = "in")
ggsave(Fig6d,filename="HIV-PAPER/Figures/Figure 6/Fig6d.pdf",
       width = 10,height = 12, dpi = 300,units = "in")


Tcell_cem@ora$ID <- str_wrap(Tcell_cem@ora$ID, width = 50)
Fig6e <- Tcell_cem@ora %>%
  dplyr::filter(Module=="M1",
                p.adjust<0.05,
                !grepl("Defective|Glutathione|Arachidonic|Ethanol|Oxidations|Prostaglandins|
           Retino|Visual|DCC|bile|O-|RA|Retinoic",ID,ignore.case=TRUE)) %>%
  dplyr::distinct(ID, .keep_all = TRUE) %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::slice_head(n = 10) %>%
  ggplot(aes(-log10(p.adjust),fct_reorder(ID,-log10(p.adjust)),
             fill=FoldEnrichment))+
  geom_col()+
  scale_fill_gradient(low = "#4A527F", high = "#931500")+
  #scale_fill_gradient(low = "grey", high = "darkblue")+
  labs(x = "-log10(padj)", 
       title = "e.",
       subtitle = "M1 over-represented pathways")+
  theme_classic2()+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.ticks.length = unit(0.5,'cm'),
        axis.title.x = element_text(size = 30),
        axis.text.y = element_text(size = 23),
        axis.text.x = element_text(size = 23),
        plot.title = element_text(size = 60, face = 'bold',hjust = -0.5),
        plot.subtitle = element_text(hjust = 0.5,size = 30))+
  guides(fill = guide_colorbar(
    title = "FoldEnrichment",
    title.theme = element_text(size = 20),
    label.theme = element_text(size = 23),
    barwidth = 1,
    barsize = 1))+
  geom_vline(xintercept = -log10(0.05), linetype='dashed')
Fig6e
ggsave(Fig6e,filename="HIV-PAPER/Figures/Figure 6/Fig6e.png",
       width = 15,height = 14, dpi = 300,units = "in")
ggsave(Fig6e,filename="HIV-PAPER/Figures/Figure 6/Fig6e.pdf",
       width = 15,height = 14, dpi = 300,units = "in")


Tcell_cem@ora$ID <- str_wrap(Tcell_cem@ora$ID, width = 50)
Fig6f <- Tcell_cem@ora %>%
  dplyr::filter(Module=="M2",
                p.adjust<0.05,
                !grepl("Rhodopsin|GPCR|ACTH",ID)) %>%
  dplyr::distinct(ID, .keep_all = TRUE) %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::slice_head(n = 10) %>%
  ggplot(aes(-log10(p.adjust),fct_reorder(ID,-log10(p.adjust)),
             fill=FoldEnrichment))+
  geom_col()+
  scale_fill_gradient(low = "#4A527F", high = "#931500")+
  #scale_fill_gradient(low = "grey", high = "darkblue")+
  labs(x = "-log10(padj)", 
       title = "f.",
       subtitle = "M2 over-represented pathways")+
  theme_classic2()+
  theme(panel.grid = element_blank(),
        legend.position = c(0.80,0.25),
        axis.title.y = element_blank(),
        axis.ticks.length = unit(0.5,'cm'),
        axis.title.x = element_text(size = 30),
        axis.text.y = element_text(size = 23),
        axis.text.x = element_text(size = 23),
        plot.title = element_text(size = 60, face = 'bold',hjust = -0.5),
        plot.subtitle = element_text(hjust = 0.5,size = 30))+
  guides(fill = guide_colorbar(
    title = "FoldEnrichment",
    title.theme = element_text(size = 30),
    label.theme = element_text(size = 23),
    barwidth = 1.5,
    barsize = 1.5,
    barheight=15.5))+
  geom_vline(xintercept = -log10(0.05), linetype='dashed')
Fig6f
ggsave(Fig6f,filename="HIV-PAPER/Figures/Figure 6/Fig6f.png",
       width = 15,height = 14, dpi = 300,units = "in")
ggsave(Fig6f,filename="HIV-PAPER/Figures/Figure 6/Fig6f.pdf",
       width = 15,height = 14, dpi = 300,units = "in")

# Figure 6g
sig_1yvsNeg_genes <- Sig_Tcell_1y_vs_Neg %>%
  as.data.frame() %>%
  dplyr::mutate(Significance = ifelse(p>-log10(0.05) & logFC>.25, "Upregulated",
                                      ifelse(p>-log10(0.05) & logFC< -.25,"Downregulated","Not significant"))) %>%
  dplyr::mutate(label = ifelse(p>-log10(0.05) & abs(logFC)>.25,rownames(.),NA)) %>%
  dplyr::filter(Significance=="Upregulated") %>%
  dplyr::select(p,logFC,rank,Significance)

sig_3mvsNeg_genes <- Sig_Tcell_3m_vs_Neg %>%
  as.data.frame() %>%
  dplyr::mutate(Significance = ifelse(p>-log10(0.05) & logFC>.25, "Upregulated",
                                      ifelse(p>-log10(0.05) & logFC< -.25,"Downregulated","Not significant"))) %>%
  dplyr::mutate(label = ifelse(p>-log10(0.05) & abs(logFC)>.25,rownames(.),NA)) %>%
  dplyr::filter(Significance=="Upregulated") %>%
  dplyr::select(p,logFC,rank,Significance)

sig_3mvs1y_genes <- Sig_Tcell_3m_vs_1y %>%
  as.data.frame() %>%
  dplyr::mutate(Significance = ifelse(p>-log10(0.05) & logFC>.25, "Upregulated",
                                      ifelse(p>-log10(0.05) & logFC< -.25,"Downregulated","Not significant"))) %>%
  dplyr::mutate(label = ifelse(p>-log10(0.05) & abs(logFC)>.25,rownames(.),NA)) %>%
  dplyr::filter(Significance=="Upregulated") %>%
  dplyr::select(p,logFC,rank,Significance)

# Exhaustion and senescence genes
list_of_exhaustion_and_senescence_associated_genes <- c("KLRG1","PDCD5","PDCD7","PDCL3","CTLA4","LAG3","EOMES","PDCD1","TIGIT","TOX","NR4A1","BATF",
                                                        "HAVCR2","CDKN1A", "CDKN2A","TP53","CD160", "ENTPD1","NR4A2","NR4A3","PRDM1","TCF7","S100A4",
                                                        "S100A10","CDKN1A", "CDKN2A", "TP53","RB1","B3GAT1")

sig_Tcell_genes_exhaustion <- rbind(sig_1yvsNeg_genes,sig_3mvsNeg_genes,sig_3mvs1y_genes) %>%
  dplyr::filter(rank %in% list_of_exhaustion_and_senescence_associated_genes) %>% unique
sig_Tcell_genes_exhaustion$rank %>%
  unique()

# T cell activation genes
Tcell_activation <- c("IL2RA","IL2","TNF","IFNG","CD69","GZMB","CD38","HLA-DRA","HLA-DRB1","CD40LG","PRF1",
                      "ICOS","FOS","JUN","NFATC1","NFATC2","BATF","IRF4","STAT1", "STAT3","STAT5A","HIF1A")
sig_Tcell_genes_activation <- rbind(sig_1yvsNeg_genes,sig_3mvsNeg_genes,sig_3mvs1y_genes) %>%
  dplyr::filter(rank %in% Tcell_activation)
sig_Tcell_genes_activation$rank %>%
  unique()

exhaustion <- c("S100A10","PDCL3","CTLA4","PDCD7","PDCD5","LAG3","KLRG1","TIGIT",
                "CDKN1A","NR4A1")
Tcell_activation <- c("STAT1","IRF4","HLA-DRA","HLA-DRB1","JUN","GZMB","FOS")
Intrinsic_apoptosis <- c("BCL2","BAK","CASP3","CASP8","TNFRSF1A","TP53", "BCL2L11","APAF1","NFKB1")
Extrinsic_apoptosis <- c("FAS","TNFRSF1A","TNFRSF10A","TNFRSF10B")

# Apoptosis associated genes
apoptosis_genes <- c("FAS","TNFRSF1A","TNFRSF10A","TNFRSF10B","BCL2","BAK",
                     "CASP3","CASP8","TNFRSF1A","TP53","BCL2L11", "APAF1")
Extrinsic_apoptosis_genes <- c("FAS","TNFRSF1A","TNFRSF10A","TNFRSF10B")
Intrinsic_apoptosis_genes <- c("BCL2","BAK","CASP3","CASP8","TNFRSF1A","TP53","BCL2L11", "APAF1")

sig_Tcell_genes_apoptosis <- rbind(sig_1yvsNeg_genes,sig_3mvsNeg_genes,sig_3mvs1y_genes) %>%
  dplyr::filter(rank %in% apoptosis_genes)
sig_Tcell_genes_apoptosis$rank %>% 
  unique()

Tcells <- AddModuleScore(Tcells,
                         features = list(exhaustion),
                         name = "Exhaustion",
                         slot = "data")

Tcells <- AddModuleScore(Tcells,
                         features = list(Tcell_activation),
                         name = "Activation",
                         slot = "data")
Tcells <- AddModuleScore(Tcells,
                         features = list(Intrinsic_apoptosis),
                         name = "Intrinsic_apoptosis",
                         slot = "data")

Tcells <- AddModuleScore(Tcells,
                         features = list(Extrinsic_apoptosis),
                         name = "Extrinsic_apoptosis",
                         slot = "data")

Fig6g <- Tcells@meta.data %>%
  ggplot(aes(HIV_Status,Exhaustion1,fill = HIV_Status))+
  geom_violin(trim = T, scale = "width")+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=43,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 't.test',
           label = "{ifelse(p < 0.0001, 'p < 0.0001', sprintf('p = %.4f', p))}",
           label.size = 10,
           tip.length = 0.01,
           step.increase = 0.2,
           hide.ns = F)+
  labs(x='',
       title = "g.",
       subtitle = "Exhaustion & Senescence genes",
       y='Module Expression')+
  scale_fill_manual(values = c('#A9A9A9','#941100','#005493'))+
  scale_x_discrete(labels=c("HIV+ ART<3 Months"="PLHIV\nART<3m",
                            "HIV+ ART>1 Year"="PLHIV\nART>1yr",
                            "HIV-"="HIV-"))+
  #scale_y_continuous(limits = c(-0.5,4.5))+
  theme_classic2()+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 30),
        axis.ticks.length = unit(0.5,'cm'),
        axis.title.x = element_text(size = 30),
        axis.text.y = element_text(size = 23),
        axis.text.x = element_text(size = 30),
        plot.title = element_text(size = 50, face = 'bold',hjust = -0.1),
        plot.subtitle = element_text(hjust = 0.5,size = 30))

Fig6g
ggsave(Fig6g,filename="HIV-PAPER/Figures/Figure 6/Fig6g.png",
       width = 9,height = 12, dpi = 300,units = "in")
ggsave(Fig6g,filename="HIV-PAPER/Figures/Figure 6/Fig6g.pdf",
       width = 9,height = 12, dpi = 300,units = "in")

Fig6h <- Tcells@meta.data %>%
  ggplot(aes(HIV_Status,Activation1,fill = HIV_Status))+
  geom_violin(trim = T, scale = "width")+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=43,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 't.test',
           label = "{ifelse(p < 0.0001, 'p < 0.0001', sprintf('p = %.4f', p))}",
           label.size = 10,
           tip.length = 0.01,
           step.increase = 0.2,
           hide.ns = F)+
  labs(x='',
       title = "h.",
       subtitle = "T cell activation genes",
       y='Module Expression')+
  scale_fill_manual(values = c('#A9A9A9','#941100','#005493'))+
  scale_x_discrete(labels=c("HIV+ ART<3 Months"="PLHIV\nART<3m",
                            "HIV+ ART>1 Year"="PLHIV\nART>1yr",
                            "HIV-"="HIV-"))+
  #scale_y_continuous(limits = c(-0.5,4.5))+
  theme_classic2()+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 30),
        axis.ticks.length = unit(0.5,'cm'),
        axis.title.x = element_text(size = 30),
        axis.text.y = element_text(size = 23),
        axis.text.x = element_text(size = 30),
        plot.title = element_text(size = 50, face = 'bold',hjust = -0.1),
        plot.subtitle = element_text(hjust = 0.5,size = 30))

Fig6h
ggsave(Fig6h,filename="HIV-PAPER/Figures/Figure 6/Fig6h.png",
       width = 9,height = 12, dpi = 300,units = "in")
ggsave(Fig6h,filename="HIV-PAPER/Figures/Figure 6/Fig6h.pdf",
       width = 9,height = 12, dpi = 300,units = "in")


Fig6i <- Tcells@meta.data %>%
  ggplot(aes(HIV_Status,Intrinsic_apoptosis1,fill = HIV_Status))+
  geom_violin(trim = T, scale = "width")+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=43,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 't.test',
           label = "{ifelse(p < 0.0001, 'p < 0.0001', sprintf('p = %.4f', p))}",
           label.size = 10,
           tip.length = 0.01,
           step.increase = 0.2,
           hide.ns = F)+
  labs(x='',
       title = "i.",
       subtitle = "Intrinsic apoptosis genes",
       y='Module Expression')+
  scale_fill_manual(values = c('#A9A9A9','#941100','#005493'))+
  scale_x_discrete(labels=c("HIV+ ART<3 Months"="PLHIV\nART<3m",
                            "HIV+ ART>1 Year"="PLHIV\nART>1yr",
                            "HIV-"="HIV-"))+
  #scale_y_continuous(limits = c(-0.5,4.5))+
  theme_classic2()+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 30),
        axis.ticks.length = unit(0.5,'cm'),
        axis.title.x = element_text(size = 30),
        axis.text.y = element_text(size = 23),
        axis.text.x = element_text(size = 30),
        plot.title = element_text(size = 50, face = 'bold',hjust = -0.1),
        plot.subtitle = element_text(hjust = 0.5,size = 30))

Fig6i
ggsave(Fig6i,filename="HIV-PAPER/Figures/Figure 6/Fig6i.png",
       width = 9,height = 12, dpi = 300,units = "in")
ggsave(Fig6i,filename="HIV-PAPER/Figures/Figure 6/Fig6i.pdf",
       width = 9,height = 12, dpi = 300,units = "in")


Fig6j <- Tcells@meta.data %>%
  ggplot(aes(HIV_Status,Extrinsic_apoptosis1,fill = HIV_Status))+
  geom_violin(trim = T, scale = "width")+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=43,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 't.test',
           label = "{ifelse(p < 0.0001, 'p < 0.0001', sprintf('p = %.4f', p))}",
           label.size = 10,
           tip.length = 0.01,
           step.increase = 0.2,
           hide.ns = F)+
  labs(x='',
       title = "j.",
       subtitle = "Extrinsic apoptosis genes",
       y='Module Expression')+
  scale_fill_manual(values = c('#A9A9A9','#941100','#005493'))+
  scale_x_discrete(labels=c("HIV+ ART<3 Months"="PLHIV\nART<3m",
                            "HIV+ ART>1 Year"="PLHIV\nART>1yr",
                            "HIV-"="HIV-"))+
  #scale_y_continuous(limits = c(-0.5,4.5))+
  theme_classic2()+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 30),
        axis.ticks.length = unit(0.5,'cm'),
        axis.title.x = element_text(size = 30),
        axis.text.y = element_text(size = 23),
        axis.text.x = element_text(size = 30),
        plot.title = element_text(size = 50, face = 'bold',hjust = -0.1),
        plot.subtitle = element_text(hjust = 0.5,size = 30))

Fig6j
ggsave(Fig6j,filename="HIV-PAPER/Figures/Figure 6/Fig6j.png",
       width = 9,height = 12, dpi = 300,units = "in")
ggsave(Fig6j,filename="HIV-PAPER/Figures/Figure 6/Fig6j.pdf",
       width = 9,height = 12, dpi = 300,units = "in")

# Save figure 6 ggsave
ggsave(filename = "HIV-PAPER/Figures/Figure 6/Fig6.png",
       plot =((Fig6a|Fig6b|Fig6c|plot_layout(ncol = 3,nrow = 1, width = c(1,1,1))))/ 
         (Fig6d|Fig6e|Fig6f|plot_layout(width = c(0.8,1,1)))/
         (Fig6g|Fig6h|Fig6i|Fig6j|plot_layout(width = c(1,1,1,1)))/
         (plot_spacer()|plot_layout(width = c(4))),
       width = 42, height = 50, units = "in", dpi = 300,limitsize = F)

ggsave(filename = "HIV-PAPER/Figures/Figure 6/Fig6.pdf",
       plot =((Fig6a|Fig6b|Fig6c|plot_layout(ncol = 3,nrow = 1, width = c(1,1,1))))/ 
         (Fig6d|Fig6e|Fig6f|plot_layout(width = c(0.8,1,1)))/
         (Fig6g|Fig6h|Fig6i|Fig6j|plot_layout(width = c(1,1,1,1)))/
         (plot_spacer()|plot_layout(width = c(4))),
       width = 42, height = 50, units = "in", dpi = 300,limitsize = F)




#################################### FIGURE 7 ##################################
# Relationship between Nasopharyngeal pneumococcal carriage and Neutrophil dynamics
Fig7a <- Figure_7a_data %>%
  ggplot(aes(x=factor(`Carriage Status`,levels=c('Carriage Negative','Carriage Positive')),
             y=as.numeric(log10(`NER`)),
             fill=factor(`Carriage Status`,
                         levels=c('Carriage Negative','Carriage Positive'))))+
  geom_violin(scale = 'width',trim = T)+
  geom_jitter(
    width = 0.1, 
    size=5
  )+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=60,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 't.test',
           label = "{ifelse(p < 0.0001, 'p < 0.0001', sprintf('p = %.4f', p))}",
           label.size = 10,
           tip.length = 0.01,
           hide.ns = F)+
  scale_fill_manual(values = c('#D6B48D','#779CCE'))+
  labs(x='',
       title = 'a.',
       y='Neutrophil to epithelial Ratio (log10)',
       fill='HIV Status')+
  scale_x_discrete(labels=c("Carriage Negative"="Carriage\nNegative",
                            "Carriage Positive"="Carriage\nPositive"))+
  facet_wrap(~ `HIV Status`,
             labeller = labeller(`HIV Status` = c("HIV-" = "HIV-",
                                                  "PLHIV ART <3m" = "PLHIV-ART<3m",
                                                  "PLHIV ART >1y" = "PLHIV-ART>yr")))+
  theme_classic2()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        legend.text = element_text(size=25), 
        strip.background = element_blank(),
        axis.ticks.length = unit(0.5,'cm'),
        plot.title = element_text(size = 50,face = 'bold',hjust = -0.07),
        strip.text.x = element_text(size = 30),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 30))
Fig7a
ggsave(Fig7a,filename="HIV-PAPER/Figures/Figure 7/Fig7a.png",
       width = 15,height = 10,dpi = 300,units = "in")  
ggsave(Fig7a,filename="HIV-PAPER/Figures/Figure 7/Fig7a.pdf",
       width = 15,height = 10,dpi = 300,units = "in")  


Fig7b <- Figure_7b_data %>%
  ggplot(aes(x=factor(`Carriage Status`,levels=c('Carriage Negative','Carriage Positive')),
             y=as.numeric(log10(`MER`)),
             fill=factor(`Carriage Status`,
                         levels=c('Carriage Negative','Carriage Positive'))))+
  geom_violin(scale = 'width',trim = T)+
  geom_jitter(
    width = 0.1, 
    size=5
  )+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=60,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 't.test',
           label = "{ifelse(p < 0.0001, 'p < 0.0001', sprintf('p = %.4f', p))}",
           label.size = 10,
           tip.length = 0.01,
           hide.ns = F)+
  scale_fill_manual(values = c('#D6B48D','#779CCE'))+
  labs(x='',
       title = 'b.',
       y='Monocyte to epithelial Ratio (log10)',
       fill='HIV Status')+
  scale_x_discrete(labels=c("Carriage Negative"="Carriage\nNegative",
                            "Carriage Positive"="Carriage\nPositive"))+
  facet_wrap(~ `HIV Status`,
             labeller = labeller(`HIV Status` = c("HIV-" = "HIV-",
                                                  "PLHIV ART <3m" = "PLHIV-ART<3m",
                                                  "PLHIV ART >1y" = "PLHIV-ART>yr")))+
  theme_classic2()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        legend.text = element_text(size=25), 
        strip.background = element_blank(),
        axis.ticks.length = unit(0.5,'cm'),
        plot.title = element_text(size = 50,face = 'bold',hjust = -0.07),
        strip.text.x = element_text(size = 30),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 30))
Fig7b
ggsave(Fig7b,filename="HIV-PAPER/Figures/Figure 7/Fig7b.png",
       width = 15,height = 10,dpi = 300,units = "in")  
ggsave(Fig7b,filename="HIV-PAPER/Figures/Figure 7/Fig7b.pdf",
       width = 15,height = 10,dpi = 300,units = "in")  

Fig7c <- Figure_7c_data %>%
  ggplot(aes(y=`Neutrophil CD11b++`,
             x=log10(Density), color=`HIV Status`))+
  geom_point(size=7)+
  geom_smooth(method = "lm")+
  stat_cor(method = "spearman", size=12, label.y = 90)+
  facet_wrap(~ `HIV Status`,
             labeller = labeller(`HIV Status` = c("HIV-" = "HIV-",
                                                  "PLHIV ART <3m" = "PLHIV-ART<3m",
                                                  "PLHIV ART >1y" = "PLHIV-ART>yr")))+
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  labs(x="Carriage Density (CFU/ML)",
       title = "c.",
       y=expression("Proportion of CD11b"^"++"~"Neutrophils"))+
  theme_classic2()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        legend.text = element_text(size=25), 
        strip.background = element_blank(),
        axis.ticks.length = unit(0.5,'cm'),
        plot.title = element_text(size = 50,face = 'bold',hjust = -0.07),
        strip.text.x = element_text(size = 30),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 30))
Fig7c

ggsave(Fig7c,filename="HIV-PAPER/Figures/Figure 7/Fig7c.png",
       width = 15,height = 10,dpi = 300,units = "in")  
ggsave(Fig7c,filename="HIV-PAPER/Figures/Figure 7/Fig7c.pdf",
       width = 15,height = 10,dpi = 300,units = "in")  


Fig7d <- Figure_7d_data %>%
  ggplot(aes(y=log10(MER),
             x=log10(Density), color=`HIV Status`))+
  geom_point(size=7)+
  geom_smooth(method = "lm")+
  stat_cor(method = "spearman", size=12, label.y =-0.2 )+
  facet_wrap(~ `HIV Status`,
             labeller = labeller(`HIV Status` = c("HIV-" = "HIV-",
                                                  "PLHIV ART <3m" = "PLHIV-ART<3m",
                                                  "PLHIV ART >1y" = "PLHIV-ART>yr")))+
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  labs(x="Carriage Density (CFU/ML)",
       title = "d.",
       y=expression("Monocyte to epithelial cell ratio (log10)"))+
  theme_classic2()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        legend.text = element_text(size=25), 
        strip.background = element_blank(),
        axis.ticks.length = unit(0.5,'cm'),
        plot.title = element_text(size = 50,face = 'bold',hjust = -0.07),
        strip.text.x = element_text(size = 30),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 30))
Fig7d
ggsave(Fig7d,filename="HIV-PAPER/Figures/Figure 7/Fig7d.png",
       width = 15,height = 10,dpi = 300,units = "in")  
ggsave(Fig7d,filename="HIV-PAPER/Figures/Figure 7/Fig7d.pdf",
       width = 15,height = 10,dpi = 300,units = "in")  


Phagocytosis_data <- read_csv("HIV-PAPER/Data/SPN-PLHIV_2025_data/Figure_5k_5l_7e_7f_data.csv")

# MFI of Oxidation
Fig7e <- Phagocytosis_data %>%
  dplyr::mutate(Oxidation_SI = (`MFI+`-`MFI-`)/(2*`rSD-`)) %>%
  dplyr::mutate(Product_Phagocytosis_MF1_Frequency = `MFI phagocytosis`*`Proportion of Phagocytosis`) %>%
  ggplot(aes(x=factor(`Carriage Status`,
                      levels=c('Carriage Negative','Carriage Positive')),
             y=log10(Product_Phagocytosis_MF1_Frequency),
             fill=factor(`Carriage Status`,
                         levels=c('Carriage Negative','Carriage Positive'))))+
  geom_violin(scale = 'width',trim = T)+
  geom_jitter(
    width = 0.1, 
    size=5
  )+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=60,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 'wilcox.test',
           label = "{ifelse(p < 0.0001, 'p < 0.0001', sprintf('p = %.4f', p))}",
           label.size = 10,
           tip.length = 0.01,
           hide.ns = F)+
  scale_fill_manual(values = c('#D6B48D','#779CCE'))+
  labs(x='',
       y="Phagocytosis Index [Log10]" ,
       #fill='HIV Status',
       title = "e.")+
  facet_wrap(~ `HIV Status`,
             labeller = labeller(`HIV Status` = c("HIV-" = "HIV-",
                                                  "PLHIV ART <3m" = "PLHIV-ART<3m",
                                                  "PLHIV ART >1y" = "PLHIV-ART>yr")))+
  scale_x_discrete(labels=c("Carriage Negative"="Carriage\nNegative",
                            "Carriage Positive"="Carriage\nPositive"))+
  theme_classic2()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        legend.text = element_text(size=25), 
        strip.background = element_blank(),
        axis.ticks.length = unit(0.5,'cm'),
        plot.title = element_text(size = 50,face = 'bold',hjust = -0.07),
        plot.subtitle = element_text(size = 30,face = 'plain',hjust = .5),
        strip.text.x = element_text(size = 30),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 30))
Fig7e
ggsave(Fig7e,filename="HIV-PAPER/Figures/Figure 7/Fig7e.png",
       width = 15,height = 10,dpi = 300,units = "in")  
ggsave(Fig7e,filename="HIV-PAPER/Figures/Figure 7/Fig7e.pdf",
       width = 15,height = 10,dpi = 300,units = "in")  

# MFI of Oxidation
Fig7f <- Phagocytosis_data %>%
  dplyr::mutate(Oxidation_SI = (`MFI+`-`MFI-`)/(2*`rSD-`)) %>%
  dplyr::mutate(Product_Oxidation_MF1_Frequency = Oxidation_SI*`Reporter out of Phagocytosed`) %>%
  ggplot(aes(x=factor(`Carriage Status`,
                      levels=c('Carriage Negative','Carriage Positive')),
             y=log10(Product_Oxidation_MF1_Frequency),
             fill=factor(`Carriage Status`,
                         levels=c('Carriage Negative','Carriage Positive'))))+
  geom_violin(scale = 'width',trim = T)+
  geom_jitter(
    #aes(color=factor(`Carriage Status`,
    #levels=c('Carriage Negative','Carriage Positive'))),
    width = 0.1, 
    size=5#, 
    #alpha=0.3
  )+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=60,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 'wilcox.test',
           label = "{ifelse(p < 0.0001, 'p < 0.0001', sprintf('p = %.4f', p))}",
           label.size = 10,
           tip.length = 0.01,
           hide.ns = F)+
  scale_fill_manual(values = c('#D6B48D','#779CCE'))+
  labs(x='',
       y="Oxidation Index [Log10]",
       #fill='HIV Status',
       title = "f.")+
  facet_wrap(~ `HIV Status`,
             labeller = labeller(`HIV Status` = c("HIV-" = "HIV-",
                                                  "PLHIV ART <3m" = "PLHIV-ART<3m",
                                                  "PLHIV ART >1y" = "PLHIV-ART>yr")))+
  scale_x_discrete(labels=c("Carriage Negative"="Carriage\nNegative",
                            "Carriage Positive"="Carriage\nPositive"))+
  theme_classic2()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        legend.text = element_text(size=25), 
        strip.background = element_blank(),
        axis.ticks.length = unit(0.5,'cm'),
        plot.title = element_text(size = 50,face = 'bold',hjust = -0.07),
        plot.subtitle = element_text(size = 30,face = 'plain',hjust = .5),
        strip.text.x = element_text(size = 30),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 30))
Fig7f
ggsave(Fig7f,filename="HIV-PAPER/Figures/Figure 7/Fig7f.png",
       width = 15,height = 10,dpi = 300,units = "in")  
ggsave(Fig7f,filename="HIV-PAPER/Figures/Figure 7/Fig7f.pdf",
       width = 15,height = 10,dpi = 300,units = "in")  


# ggsave 
ggsave(filename = "HIV-PAPER/Figures/Figure 7/Fig7.png",
       plot = ((Fig7a|Fig7b|plot_layout(ncol = 2,nrow = 1, width = c(1,1))))/ 
         (Fig7c|Fig7d|plot_layout(width = c(1,1)))/
         (Fig7e|Fig7f|plot_layout(width = c(1,1))),
       width = 32, height = 31, units = "in", dpi = 300)

ggsave(filename = "HIV-PAPER/Figures/Figure 7/Fig7.pdf",
       plot = ((Fig7a|Fig7b|plot_layout(ncol = 2,nrow = 1, width = c(1,1))))/ 
         (Fig7c|Fig7d|plot_layout(width = c(1,1)))/
         (Fig7e|Fig7f|plot_layout(width = c(1,1))),
       width = 32, height = 31, units = "in", dpi = 300)
