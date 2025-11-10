#################################### LOADING REQUIRED PACKAGES #################
pacman::p_load(char = c("lubridate","gtsummary", "tidyverse", "dplyr", "here", "rio", "scales", "boot", "Matrix","ggpubr",
                        "magrittr",  "mvtnorm", "zoo", "patchwork", "mgcv", "PropCIs", "writexl","DropletUtils","SeuratWrappers",
                        "Seurat","rjson","R2HTML","DT","cowplot","RCurl","glmGamPoi","DESeq2","ggrepel","rPanglaoDB","ComplexHeatmap",
                        "EnhancedVolcano","RColorBrewer","circlize","rmarkdown","biomaRt","biomartr","clusterProfiler","multinichenetr",
                        "AnnotationDbi","org.Hs.eg.db","CEMiTool","enrichplot","pathview","scmap","SingleR","S4Vectors","TMB","muscat",
                        "SingleCellExperiment","apeglm","edgeR","purrr","tibble","png","RColorBrewer","scran","microViz","CommPath",
                        "reshape2","scater","Azimuth","scCATCH","CellChat","SoupX","knitr","DoubletFinder","ggpubr","viridis","GSVA",
                        "ggmin","cluster","foreach","doParallel","BPCells","ggimage","ggbeeswarm","grid","data.table","scriabin",
                        "clusterExperiment","destiny","gam","corrplot","ggthemes","base64enc","Biobase","CATALYST","dittoSeq","viridis",
                        "DelayedArray","DelayedMatrixStats","limma","lme4","batchelor","HDF5Array","terra","ggrastr","nloptr","ggsignif",
                        'lubridate',"gtsummary", 'tidyverse', "dplyr", "here", "rio", "scales", "boot", 
                        "magrittr",  "mvtnorm", "zoo", "patchwork", "mgcv", "PropCIs", "writexl", "presto",
                        "ggsignif", "ggpubr", "ggeasy", "cowplot","ggExtra", "PupillometryR","hrbrthemes", "ggstance",
                        "survival","survminer","sysfonts","showtext","nlme",'glue'))

#################################### SET WORKING DIRECTORY #####################
setwd("/Users/471288/Library/CloudStorage/OneDrive-Malawi-LiverpoolWellcomeTrust/Projects/PhD Work/Nasomune/Thesis Paper Publications/HIV Paper")

#################################### PREPARING INPUT SAMPLES ###################
# read in Micro.csv and MPO NETs .csv data file
Clinical_Data <- read_csv("Data/Main_Files_Thesis/Clinical_Data.csv")
Lab_Data <- read_csv("Data/Main_Files_Thesis/Lab_Data.csv")
Micro <- read_csv("Data/Main_Files_Thesis/Micro.csv")
NETSMPO <- read_csv("Data/Main_Files_Thesis/MPO_NETs.csv")
Cytokines_averages <- read_csv("Data/Main_Files_Thesis/Cytokines_averages.csv") 
Innate_Like_2024 <- read_csv("Data/Main_Files_Thesis/Innate-Like_2024.csv") %>%
  dplyr::select(-Sample) %>%
  dplyr::filter(`Sample Type`=="Nasal") %>%
  dplyr::select(-`Sample Type`,-`Sample Quality`)
# List all .csv files in the folder
CD3_files <- list.files(path = "Data/Main_Files_thesis", pattern = "_HIVPaper_CD3.csv", full.names = TRUE)
CD3_files
# read the CSV files
CD3_files <- lapply(CD3_files, read_csv)
# Combine all the main data
CD3_files <- do.call(rbind, CD3_files) %>%
  dplyr::select(-Sample)
# List all .csv files in the folder
Maincsv_files <- list.files(path = "Data/Main_Files_thesis", pattern = "_Main.csv", full.names = TRUE)
Maincsv_files
# read the CSV files
Maindata_list <- lapply(Maincsv_files, read_csv)
# Combine all the main data
Maindata <- do.call(rbind, Maindata_list) %>%
  dplyr::select(-Sample) %>%
  dplyr::mutate(NER = `Neutrophil count`/`Epithelial count`,
                MER = `Monocyte count`/`Epithelial count`,
                LER = `CD45+ count`/`Epithelial count`) %>%
  dplyr::mutate(NCD45ER = `Neutrophil count`/`LER`,
                MCD45ER = `Monocyte count`/`LER`,
                CD11bER = `Neutrophil CD11b++`/`Epithelial count`) %>%
  dplyr::filter(`Sample Type`=="Nasal",
                `Sample Quality`=="Good") %>%
  dplyr::select(-`Sample Quality`,-`Sample Type`,-Panel)
Maindata$`LAB ID` <- gsub("CHU106","CUH106", Maindata$`LAB ID`)
Maindata$`LAB ID` <- gsub("CUF12BN1","CUF12B", Maindata$`LAB ID`)
# List all .csv files in the folder containing NeutrophilCD10CD63
NeutCD10CD63csv_files <- list.files(path = "Data/Main_Files_thesis", pattern = "_NeutrophilCD10CD63.csv", full.names = TRUE)
NeutCD10CD63csv_files
# read the CSV files
NeutCD10CD63data_list <- lapply(NeutCD10CD63csv_files, read_csv)
# Combine all the NeutCD10CD63 data
NeutCD10CD63 <- do.call(rbind, NeutCD10CD63data_list) %>%
  dplyr::select(-Sample,-`Sample Type`,-`Sample Quality`,-`CD3 Staining`,-`CD11b Staining`,-`Panel`)
# List all .csv files in the folder containing MonocyteCD10CD63
MonoCD10CD63csv_files <- list.files(path = "Data/Main_Files_thesis", pattern = "_MonocyteCD10CD63.csv", full.names = TRUE)
MonoCD10CD63csv_files
# read the CSV files
MonoCD10CD63data_list <- lapply(MonoCD10CD63csv_files, read_csv)
# Combine all the MonoCD10CD63 data
MonoCD10CD63 <- do.call(rbind, MonoCD10CD63data_list) %>%
  dplyr::select(-Sample,-`Sample Type`,-`Sample Quality`,-`CD3 Staining`,-`CD11b Staining`,-`Panel`)
# Merge all the files
Masterfile <- Maindata %>%
  merge(NeutCD10CD63, by=c("LAB ID"),all=T) %>%
  merge(MonoCD10CD63, by=c("LAB ID"),all=T) %>%
  merge(NETSMPO, by=c("LAB ID"),all=T) %>%
  merge(Cytokines_averages, by=c("LAB ID"),all=T) %>%
  merge(Innate_Like_2024, by=c("LAB ID"),all = T) %>%
  merge(Lab_Data, by=c("LAB ID"),all = T) %>%
  merge(Clinical_Data, by=c("LAB ID"),all = T) %>%
  merge(Micro,by=c("LAB ID"),all=T) %>%
  merge(CD3_files, by=c("LAB ID"),all = T)
Masterfile <- Masterfile[!duplicated(Masterfile$`LAB ID`), ]
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
Fig1b_counts <- Masterfile %>%
  dplyr::filter(`HIV Status` == "HIV-",
                `Epithelial count` > 100,
                `CD45+ count` > 200,
                `Visit` == 'Week 1',
                `CD3 Staining` == 'Stained') %>%
  dplyr::mutate(TER = `T cell count` / `Epithelial count`) %>%
  dplyr::select(NER,MER,TER) %>%
  pivot_longer(cols = c("NER", "TER", "MER"),
               names_to = "Immune cells",
               values_to = "Immune cell to epithelial cell ratio") %>%
  dplyr::group_by(`Immune cells`) %>%
  dplyr::summarize(n = dplyr::n(), .groups = "drop")
Fig1b_counts

Fig1b <- Masterfile %>%
  filter(`HIV Status`=="HIV-",
         `Epithelial count`>100,
         `CD45+ count`>200,
         `Visit`=='Week 1',
         `CD3 Staining`=='Stained'
  ) %>%
  mutate(TER=`T cell count`/`Epithelial count`) %>%
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

ggsave(Fig1b,filename="HIV-PAPER/Figures/Figure 1/Fig1b.png",
       width = 9,height = 10,dpi = 300)
ggsave(Fig1b,filename="HIV-PAPER/Figures/Figure 1/Fig1b.pdf",
       width = 9,height = 10,dpi = 300)

###################################################################### Figure 1c
Fig1c_counts <- Masterfile %>%
  filter(`Epithelial count`>100,
         `CD45+ count`>200,
         `Visit`=="Week 1") %>%
  dplyr::select(NER,`HIV Status`) %>%
  dplyr::group_by(`HIV Status`) %>%
  dplyr::summarize(n = dplyr::n(), .groups = "drop")
Fig1c_counts

Fig1c <- Masterfile %>%
  filter(`Epithelial count`>100,
         `CD45+ count`>200,
         `Visit`=="Week 1") %>%
  group_by(`HIV Status`,`Carriage Status`) %>%
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

ggsave(Fig1c,filename="HIV-PAPER/Figures/Figure 1/Fig1c.png",
       width = 8,height = 10,dpi = 300)
ggsave(Fig1c,filename="HIV-PAPER/Figures/Figure 1/Fig1c.pdf",
       width = 8,height = 10,dpi = 300)

###################################################################### Figure 1d
HIV_status_colors <- c('#A9A9A9','#941100','#005493')
Fig1d_counts <- Masterfile %>%
  filter(`Epithelial count`>100,
         `CD45+ count`>200,
         `Visit`=="Week 1") %>%
  dplyr::select(NER,`HIV Status`) %>%
  dplyr::group_by(`HIV Status`) %>%
  dplyr::summarize(n = dplyr::n(), .groups = "drop")
Fig1d_counts

Fig1d <- Masterfile %>%
  filter(`Epithelial count`>100,
         `CD45+ count`>200,
         `Visit`=="Week 1") %>%
  group_by(`HIV Status`,`Carriage Status`) %>%
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
  geom_pwc(method = 't.test',
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

ggsave(Fig1d,filename="HIV-PAPER/Figures/Figure 1/Fig1d.png",
       width = 8,height = 10,dpi = 300)
ggsave(Fig1d,filename="HIV-PAPER/Figures/Figure 1/Fig1d.pdf",
       width = 8,height = 10,dpi = 300)


###################################################################### Figure 1e
Fig1e_counts <- Masterfile %>%
  dplyr::filter(`CD4+`<80,
                `HIV Status`!='NA',
                Visit=="Week 1") %>%
  dplyr::select(`CD4+`,`HIV Status`) %>%
  dplyr::group_by(`HIV Status`) %>%
  dplyr::summarize(n = dplyr::n(), .groups = "drop")
Fig1e_counts

Fig1e <- Masterfile %>%
  filter(`CD4+`<80,
         `HIV Status`!='NA',
         Visit=="Week 1") %>%
  group_by(`HIV Status`,`Carriage Status`) %>%
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
  scale_y_continuous(limits = c(0,60),breaks = c(0,10,20,30,40,50,60))+
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

ggsave(Fig1e,filename="HIV-PAPER/Figures/Figure 1/Fig1e.png",
       width = 8,height = 10,dpi = 300)
ggsave(Fig1e,filename="HIV-PAPER/Figures/Figure 1/Fig1e.pdf",
       width = 8,height = 10,dpi = 300)

###################################################################### Figure 1f
Fig1f_counts <- Masterfile %>%
  dplyr::filter(`CD8+`>20,
                Visit=="Week 1",
                `HIV Status`!="NA") %>%
  dplyr::select(`CD8+`,`HIV Status`) %>%
  dplyr::group_by(`HIV Status`) %>%
  dplyr::summarize(n = dplyr::n(), .groups = "drop")
Fig1f_counts

Fig1f <- Masterfile %>%
  filter(`CD8+`>20,
         Visit=="Week 1",
         `HIV Status`!="NA") %>%
  group_by(`HIV Status`,`Carriage Status`) %>%
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

ggsave(Fig1f,filename="HIV-PAPER/Figures/Figure 1/Fig1f.png",
       width = 8,height = 10,dpi = 300)
ggsave(Fig1f,filename="HIV-PAPER/Figures/Figure 1/Fig1f.pdf",
       width = 8,height = 10,dpi = 300)


###################################################################### Figure 1g
Fig1g_counts <- Masterfile %>%
  dplyr::filter(Visit=="Week 1",
                `HIV Status`!="NA") %>%
  dplyr::select(`CD3+Mait`,`HIV Status`) %>%
  dplyr::group_by(`HIV Status`) %>%
  na.omit() %>%
  dplyr::summarize(n = dplyr::n(), .groups = "drop")
Fig1g_counts

Fig1g <- Masterfile %>%
  filter(`HIV Status`!='NA',
         Visit=="Week 1") %>%
  group_by(`HIV Status`,`Carriage Status`) %>%
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
ggsave(Fig1g,filename="HIV-PAPER/Figures/Figure 1/Fig1g.png",
       width = 8,height = 10,dpi = 300)
ggsave(Fig1g,filename="HIV-PAPER/Figures/Figure 1/Fig1g.pdf",
       width = 8,height = 10,dpi = 300)

###################################################################### Figure 1h
Fig1h_counts <- Masterfile %>%
  dplyr::filter(Visit=="Week 1",
                `HIV Status`!="NA") %>%
  dplyr::select(`CD3+TCRgd+`,`HIV Status`) %>%
  dplyr::group_by(`HIV Status`) %>%
  na.omit() %>%
  dplyr::summarize(n = dplyr::n(), .groups = "drop")
Fig1h_counts

Fig1h <- Masterfile %>%
  filter(`HIV Status`!='NA',
         Visit=="Week 1") %>%
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
ggsave(Fig1h,filename="HIV-PAPER/Figures/Figure 1/Fig1h.png",
       width = 8,height = 10,dpi = 300)
ggsave(Fig1h,filename="HIV-PAPER/Figures/Figure 1/Fig1h.pdf",
       width = 8,height = 10,dpi = 300)

###################################################################### Figure 1i
Fig1i_counts <- Masterfile %>%
  dplyr::filter(Visit=="Week 1",
                `HIV Status`!="NA") %>%
  dplyr::select(`CD3+CD56+`,`HIV Status`) %>%
  dplyr::group_by(`HIV Status`) %>%
  na.omit() %>%
  dplyr::summarize(n = dplyr::n(), .groups = "drop")
Fig1i_counts

Fig1i <- Masterfile %>%
  filter(`HIV Status`!='NA',
         Visit=="Week 1") %>%
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

ggsave(Fig1i,filename="HIV-PAPER/Figures/Figure 1/Fig1i.png",
       width = 8,height = 10,dpi = 300)
ggsave(Fig1i,filename="HIV-PAPER/Figures/Figure 1/Fig1i.pdf",
       width = 8,height = 10,dpi = 300)

###################################################################### Save Fig1
ggsave(filename = "HIV-PAPER/Figures/Figure 1/Fig1.pdf",
       plot = ((Fig1a|Fig1b|plot_layout(ncol = 2,nrow = 1, width = c(3,1))))/ 
         (Fig1c|Fig1d|Fig1e|Fig1f|plot_layout(width = c(1,1,1,1)))/
         (Fig1g|Fig1h|Fig1i|plot_layout(width = c(1,1,1,1))),
       width = 32, height = 29, units = "in", dpi = 300)

ggsave(filename = "HIV-PAPER/Figures/Figure 1/Fig1.png",
       plot = ((Fig1a|Fig1b|plot_layout(ncol = 2,nrow = 1, width = c(3,1))))/ 
         (Fig1c|Fig1d|Fig1e|Fig1f|plot_layout(width = c(1,1,1,1)))/
         (Fig1g|Fig1h|Fig1i|plot_layout(width = c(1,1,1,1))),
       width = 32, height = 29, units = "in", dpi = 300)

#################################### EXTENDED FIGURE 1 #########################
############################################################# Extended Figure 1a
Extended_Fig1a <- magick::image_read("New_Figures/Extended_Figure 1/Extended_Figure1a.png")
Extended_Fig1a <- ggdraw()+ 
  draw_image(Extended_Fig1a, scale = 1.0)+
  labs(title = "a.",
       subtitle = "Myeloid panel gating strategy")+
  theme(
    plot.title = element_text(hjust = 0.01,vjust = .01, size = 50, face = "bold"),  # Customize title
    plot.subtitle = element_text(hjust = 0.5, size = 50),  # Customize title
    plot.margin = unit(c(0,0,0,0), "cm")  # Optional: adjust margins
  )
Extended_Fig1a
############################################################# Extended Figure 1b
Extended_Fig1b <- magick::image_read("New_Figures/Extended_Figure 1/Extended_Fig1b_Tcell_gating.png")
Extended_Fig1b <- ggdraw()+ 
  draw_image(Extended_Fig1b, scale = 1.0)+
  labs(title = "b.",
       subtitle = "T cell panel gating strategy")+
  theme(
    plot.title = element_text(hjust = 0.0,vjust = .01, size = 50, face = "bold"),  # Customize title
    plot.subtitle = element_text(hjust = 0.5, size = 50),  # Customize title
    plot.margin = unit(c(0,0,0,0), "cm")  # Optional: adjust margins
  )
Extended_Fig1b

############################################################# Save Extended_Fig1
ggsave(filename = "HIV-PAPER/Figures/Extended_Figure 1/Extended_Fig1.pdf",
       plot = ((Extended_Fig1a|plot_spacer()|plot_layout(ncol = 2,nrow = 1, width = c(2,0.001))))/
         (Extended_Fig1b|plot_layout(width = c(2,0.001)))/
         (plot_spacer()|plot_layout(width = c(1,1,1,1)))/
         (plot_spacer()|plot_layout(width = c(1))),  
       width = 32, height = 29, units = "in", dpi = 300)

ggsave(filename = "HIV-PAPER/Figures/Extended_Figure 1/Extended_Fig1.png",
       plot = ((Extended_Fig1a|plot_spacer()|plot_layout(ncol = 2,nrow = 1, width = c(2,0.001))))/
         (Extended_Fig1b|plot_layout(width = c(2,0.001)))/
         (plot_spacer()|plot_layout(width = c(1,1,1,1)))/
         (plot_spacer()|plot_layout(width = c(1))),  
       width = 32, height = 29, units = "in", dpi = 300)


#################################### FIGURE 2 ##################################
# Figure 2a
# Load data
Neg <- magick::image_read("New_Figures/Figure 2/Neg.png")
Neg <- ggdraw()+ 
  draw_image(Neg, scale = 1.0)+
  labs(title = "f.",
       subtitle = "HIV-")+
  theme(
    plot.title = element_text(hjust = 0.10,vjust = -12, size = 40, face = "bold"),
    plot.subtitle = element_text(hjust = 0.55, size = 40, vjust = -72,face = "plain"),# Customize title
    plot.margin = unit(c(0.0, 0.0, 0.0, 0.0), "cm")  # Optional: adjust margins
  )
Neg

M3 <- magick::image_read("New_Figures/Figure 2/3M.png")
M3 <- ggdraw()+ 
  draw_image(M3, scale = 1.0)+
  labs(title = "g.",
       subtitle = "PLHIV-ART<3m")+
  theme(
    plot.title = element_text(hjust = 0.10,vjust = -12, size = 40, face = "bold"),
    plot.subtitle = element_text(hjust = 0.7, size = 40, vjust = -72,face = "plain"),# Customize title
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")  # Optional: adjust margins
  )
M3


Y1 <- magick::image_read("New_Figures/Figure 2/1y.png")
Y1 <- ggdraw()+ 
  draw_image(Y1, scale = 1.0)+
  labs(title = "h.",
       subtitle = "PLHIV-ART>1yr")+
  theme(
    plot.title = element_text(hjust = 0.10,vjust = -12, size = 40, face = "bold"),
    plot.subtitle = element_text(hjust = 0.7, size = 40, vjust = -72,face = "plain"),# Customize title
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")  # Optional: adjust margins
  )
Y1

Cytokines <- read_csv("Data/Main_Files_Thesis/Cytokines_averages_zeros_replaced_with_minimum_values.csv")
Micro <- read_csv("Data/Main_Files_Thesis/Micro.csv")

Neutrophil_Monocytes_data <- Masterfile %>%
  dplyr::select(`LAB ID`,NER, MER, `CD45+ count`, `Epithelial count`)

colnames(Neutrophil_Monocytes_data) <- c("LAB ID", "Neutrophils", "Monocytes", "CD45+ count", "Epithelial count")

Cytokines <- Cytokines %>%
  merge(Micro, by=c("LAB ID"),all=F) %>%
  merge(Neutrophil_Monocytes_data, by=c("LAB ID"),all=F)


# Cytokine correlations
Cytokine_correlations <- Cytokines %>%
  # Filter for Week 1
  dplyr::filter(Visit == "Week 1") %>%
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

png("Final_Paper_Figures/Figure 2/HIV-correlogram.png", width = 800, height = 800)
# Add correlation plot
col_grad <- colorRampPalette(c("blue", "white", "red"))
corrplot(cor_results$`HIV-`$r,
         method = "circle",
         type = "lower",
         col = col_grad(200),
         #title = paste("HIV-"),
         mar = c(0, 0, 2, 0),
         cl.lim = c(-1, 1),
         addgrid.col = "black",
         tl.cex = 1.9,
         tl.col = "black",  
         number.cex = 0.7,
         p.mat = cor_results$`HIV-`$P,            # Add p-value matrix
         sig.level = 0.05,        # Highlight significant correlations
         insig = "label_sig")     # Show p-values as label
dev.off()

# Add correlation plot
png("Final_Paper_Figures/Figure 2/PLHIV-ART<3mcorrelogram.png", width = 800, height = 800)
corrplot(cor_results$`PLHIV ART <3m`$r,
         method = "circle",
         type = "lower",
         col = col_grad(200),
         #title = paste("ART<3m"),
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
png("Final_Paper_Figures/Figure 2/PLHIV-ART>1yrcorrelogram.png", width = 800, height = 800)
corrplot(cor_results$`PLHIV ART >1y`$r,
         method = "circle",
         type = "lower",
         col = col_grad(200),
         #title = paste("ART>1y"),
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


`````
Neg <- magick::image_read("New_Figures/Figure 2/Neg.png")
Neg <- ggdraw()+ 
  draw_image(Neg, scale = 1.0)+
  labs(title = "f.",
       subtitle = "HIV-")+
  theme(
    plot.title = element_text(hjust = 0.10,vjust = -12, size = 40, face = "bold"),
    plot.subtitle = element_text(hjust = 0.55, size = 40, vjust = -72,face = "plain"),# Customize title
    plot.margin = unit(c(0.0, 0.0, 0.0, 0.0), "cm")  # Optional: adjust margins
  )
Neg
````

M3 <- magick::image_read("New_Figures/Figure 2/3M.png")
M3 <- ggdraw()+ 
  draw_image(M3, scale = 1.0)+
  labs(title = "g.",
       subtitle = "PLHIV-ART<3m")+
  theme(
    plot.title = element_text(hjust = 0.10,vjust = -12, size = 40, face = "bold"),
    plot.subtitle = element_text(hjust = 0.7, size = 40, vjust = -72,face = "plain"),# Customize title
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")  # Optional: adjust margins
  )
M3
````


Y1 <- magick::image_read("New_Figures/Figure 2/1y.png")
Y1 <- ggdraw()+ 
  draw_image(Y1, scale = 1.0)+
  labs(title = "h.",
       subtitle = "PLHIV-ART>1yr")+
  theme(
    plot.title = element_text(hjust = 0.10,vjust = -12, size = 40, face = "bold"),
    plot.subtitle = element_text(hjust = 0.7, size = 40, vjust = -72,face = "plain"),# Customize title
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")  # Optional: adjust margins
  )
Y1
````


ggsave(filename = "New_Figures/Figure 2/Fig2f-h.png",
       plot = ((Neg|M3|Y1|plot_layout(ncol = 3,nrow = 1, width = c(1,1,1))))/
         (plot_spacer()|plot_layout(width = c(1,1,1)))/
         (plot_spacer()|plot_layout(width = c(1,1,1))),
       width = 32, height = 29, units = "in", dpi = 300)

ggsave(filename = "New_Figures/Figure 2/Fig2f-h.png",
       plot = ((Neg|M3|Y1|plot_layout(ncol = 3,nrow = 1, width = c(1,1,1))))/
         (plot_spacer()|plot_layout(width = c(1,1,1)))/
         (plot_spacer()|plot_layout(width = c(1,1,1))),
       width = 25, height = 29, units = "in", dpi = 300)

ggsave(filename = "New_Figures/Figure 2/Fig2f-h.pdf",
       plot = ((Neg|M3|Y1|plot_layout(ncol = 3,nrow = 1, width = c(1,1,1))))/
         (plot_spacer()|plot_layout(width = c(1,1,1)))/
         (plot_spacer()|plot_layout(width = c(1,1,1))),
       width = 25, height = 29, units = "in", dpi = 300)
````
#################################### EXTENDED FIGURE 2 ########################
# Load data
Cytokines <- read_csv("Data/Main_Files_Thesis/Cytokines_averages_zeros_replaced_with_minimum_values.csv")
Micro <- read_csv("Data/Main_Files_Thesis/Micro.csv")

Neutrophil_Monocytes_data <- Masterfile %>%
  dplyr::select(`LAB ID`,NER, MER, `CD45+ count`, `Epithelial count`)

Cytokines <- Cytokines %>%
  merge(Micro, by=c("LAB ID"),all=F) %>%
  merge(Neutrophil_Monocytes_data, by=c("LAB ID"),all=F)

write.csv(Cytokines,"HIV-PAPER/Data/Cytokines_used_in_paper.csv")

group_counts <- Masterfile %>%
  dplyr::filter(Visit=="Week 1",
                `HIV Status`!="NA") %>%
  dplyr::select(`CD3+CD56+`,`HIV Status`) %>%
  dplyr::group_by(`HIV Status`) %>%
  na.omit() %>%
  dplyr::summarize(n = dplyr::n(), .groups = "drop")
group_counts

# Extended figure 2a
Extended_Fig2a <- Masterfile %>%
  filter(`CD11b Staining`=="Stained",
         #`Visit`=="Week 1",
         `HIV Status`=="HIV-",
         `Neutrophil CD10+CD11b+CD62L-CD63-`>20,
         `Neutrophil CD11b++`>20 &
           `Neutrophil CD11b++`<90
  ) %>%
  ggplot(aes(`Neutrophil CD10+CD11b+CD62L-CD63-`,`Neutrophil CD11b++`))+
  geom_point(size=7)+
  geom_smooth(method = "lm")+
  stat_cor(method = "spearman",size=12)+
  labs(x= expression(atop("Proportion of CD10"^"+"~"CD11b"^"++"~"CD63"^"-"~"CD62L"^"-", 
                          "Neutrophils"
  )
  ),
  y=expression("Proportion of CD11b"^"++"~" Neutrophils"),
  title = "a.")+
  scale_y_continuous(limits = c(50,75),breaks = c(50,60,70))+
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
Extended_Fig2a
ggsave(Extended_Fig2a,filename="HIV-PAPER/Figures/Extended_Figure 2/Extended_Fig2a.png",
       width = 8,height = 12,dpi = 300)
ggsave(Extended_Fig2a,filename="HIV-PAPER/Figures/Extended_Figure 2/Extended_Fig2a.pdf",
       width = 8,height = 12,dpi = 300)


# Each cytokine 
# Figure 2b
TNFa_counts <- Cytokines %>%
  filter(`Visit`=="Week 3") %>%
  dplyr::select(`TNF-a`,`HIV Status`) %>%
  dplyr::group_by(`HIV Status`) %>%
  dplyr::summarize(n = dplyr::n(), .groups = "drop")
TNFa_counts

Extended_Fig2b <- Cytokines %>%
  dplyr::filter(`Visit` %in% c("Week 3")) %>%
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
Extended_Fig2b
ggsave(Extended_Fig2b,filename="HIV-PAPER/Figures/Extended_Figure 2/Extended_Fig2b.png",
       width = 8,height = 12,dpi = 300)
ggsave(Extended_Fig2b,filename="HIV-PAPER/Figures/Extended_Figure 2/Extended_Fig2b.pdf",
       width = 8,height = 12,dpi = 300)

# Figure 2c
IL6_counts <- Cytokines %>%
  filter(`Visit`=="Week 3") %>%
  dplyr::select(`IL-6`,`HIV Status`) %>%
  dplyr::group_by(`HIV Status`) %>%
  dplyr::summarize(n = dplyr::n(), .groups = "drop")
IL6_counts

Extended_Fig2c <- Cytokines %>%
  dplyr::filter(Visit=="Week 3") %>%
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
Extended_Fig2c
ggsave(Extended_Fig2c,filename="HIV-PAPER/Figures/Extended_Figure 2/Extended_Fig2c.png",
       width = 8,height = 12,dpi = 300)
ggsave(Extended_Fig2c,filename="HIV-PAPER/Figures/Extended_Figure 2/Extended_Fig2c.pdf",
       width = 8,height = 12,dpi = 300)
# Figure 2d
IL8_counts <- Cytokines %>%
  filter(`Visit`=="Week 3") %>%
  dplyr::select(`IL-8`,`HIV Status`) %>%
  dplyr::group_by(`HIV Status`) %>%
  dplyr::summarize(n = dplyr::n(), .groups = "drop")
IL8_counts

Extended_Fig2d <- Cytokines %>%
  dplyr::filter(`Visit`=="Week 3") %>%
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
Extended_Fig2d
ggsave(Extended_Fig2d,filename="HIV-PAPER/Figures/Extended_Figure 2/Extended_Fig2d.png",
       width = 8,height = 12,dpi = 300)
ggsave(Extended_Fig2d,filename="HIV-PAPER/Figures/Extended_Figure 2/Extended_Fig2d.pdf",
       width = 8,height = 12,dpi = 300)

# Figure 2e
IL1B_counts <- Cytokines %>%
  filter(`Visit`=="Week 3") %>%
  dplyr::select(`IL-1B`,`HIV Status`) %>%
  dplyr::group_by(`HIV Status`) %>%
  dplyr::summarize(n = dplyr::n(), .groups = "drop")
IL1B_counts

Extended_Fig2e <- Cytokines %>%
  dplyr::filter(Visit=="Week 3") %>%
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
Extended_Fig2e
ggsave(Extended_Fig2e,filename="HIV-PAPER/Figures/Extended_Figure 2/Extended_Fig2e.png",
       width = 8,height = 12,dpi = 300)
ggsave(Extended_Fig2e,filename="HIV-PAPER/Figures/Extended_Figure 2/Extended_Fig2e.pdf",
       width = 8,height = 12,dpi = 300)

# Figure 2f
IFNy_counts <- Cytokines %>%
  filter(`Visit`=="Week 3") %>%
  dplyr::select(`IFN-y`,`HIV Status`) %>%
  dplyr::group_by(`HIV Status`) %>%
  dplyr::summarize(n = dplyr::n(), .groups = "drop")
IFNy_counts

Extended_Fig2f <- Cytokines %>%
  dplyr::filter(Visit=="Week 3") %>%
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
Extended_Fig2f
ggsave(Extended_Fig2f,filename="HIV-PAPER/Figures/Extended_Figure 2/Extended_Fig2f.png",
       width = 8,height = 12,dpi = 300)
ggsave(Extended_Fig2f,filename="HIV-PAPER/Figures/Extended_Figure 2/Extended_Fig2f.pdf",
       width = 8,height = 12,dpi = 300)

# Figure 2g
IL2_counts <- Cytokines %>%
  filter(`Visit`=="Week 3") %>%
  dplyr::select(`IL-2`,`HIV Status`) %>%
  dplyr::group_by(`HIV Status`) %>%
  dplyr::summarize(n = dplyr::n(), .groups = "drop")
IL2_counts

# Interleukin-2
Extended_Fig2g <- Cytokines %>%
  dplyr::filter(Visit=="Week 3") %>%
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
Extended_Fig2g
ggsave(Extended_Fig2g,filename="HIV-PAPER/Figures/Extended_Figure 2/Extended_Fig2g.png",
       width = 8,height = 12,dpi = 300)
ggsave(Extended_Fig2g,filename="HIV-PAPER/Figures/Extended_Figure 2/Extended_Fig2g.pdf",
       width = 8,height = 12,dpi = 300)


# Figure 2h
IL12_counts <- Cytokines %>%
  filter(`Visit`=="Week 3") %>%
  dplyr::select(`IL-12p70`,`HIV Status`) %>%
  dplyr::group_by(`HIV Status`) %>%
  dplyr::summarize(n = dplyr::n(), .groups = "drop")
IL12_counts

Extended_Fig2h <- Cytokines %>%
  dplyr::filter(Visit=="Week 3") %>%
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
Extended_Fig2h
ggsave(Extended_Fig2h,filename="HIV-PAPER/Figures/Extended_Figure 2/Extended_Fig2h.png",
       width = 8,height = 12,dpi = 300)
ggsave(Extended_Fig2h,filename="HIV-PAPER/Figures/Extended_Figure 2/Extended_Fig2h.pdf",
       width = 8,height = 12,dpi = 300)


# Figure 2i
IL4_counts <- Cytokines %>%
  filter(`Visit`=="Week 3") %>%
  dplyr::select(`IL-4`,`HIV Status`) %>%
  dplyr::group_by(`HIV Status`) %>%
  dplyr::summarize(n = dplyr::n(), .groups = "drop")
IL4_counts

Extended_Fig2i <- Cytokines %>%
  dplyr::filter(Visit=="Week 3") %>%
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
Extended_Fig2i
ggsave(Extended_Fig2i,filename="HIV-PAPER/Figures/Extended_Figure 2/Extended_Fig2i.png",
       width = 8,height = 12,dpi = 300)
ggsave(Extended_Fig2i,filename="HIV-PAPER/Figures/Extended_Figure 2/Extended_Fig2i.pdf",
       width = 8,height = 12,dpi = 300)


# Figure 2j
IL13_counts <- Cytokines %>%
  filter(`Visit`=="Week 3") %>%
  dplyr::select(`IL-13`,`HIV Status`) %>%
  dplyr::group_by(`HIV Status`) %>%
  dplyr::summarize(n = dplyr::n(), .groups = "drop")
IL13_counts

Extended_Fig2j <- Cytokines %>%
  dplyr::filter(Visit=="Week 3") %>%
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
Extended_Fig2j
ggsave(Extended_Fig2j,filename="HIV-PAPER/Figures/Extended_Figure 2/Extended_Fig2j.png",
       width = 8,height = 12,dpi = 300)
ggsave(Extended_Fig2j,filename="HIV-PAPER/Figures/Extended_Figure 2/Extended_Fig2j.pdf",
       width = 8,height = 12,dpi = 300)

# Figure 2k
IL10_counts <- Cytokines %>%
  filter(`Visit`=="Week 3") %>%
  dplyr::select(`IL-10`,`HIV Status`) %>%
  dplyr::group_by(`HIV Status`) %>%
  dplyr::summarize(n = dplyr::n(), .groups = "drop")
IL10_counts

Extended_Fig2k <- Cytokines %>%
  dplyr::filter(Visit=="Week 3") %>%
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
Extended_Fig2k
ggsave(Extended_Fig2k,filename="HIV-PAPER/Figures/Extended_Figure 2/Extended_Fig2k.png",
       width = 8,height = 12,dpi = 300)
ggsave(Extended_Fig2k,filename="HIV-PAPER/Figures/Extended_Figure 2/Extended_Fig2k.pdf",
       width = 8,height = 12,dpi = 300)


# Save figure 2 ggsave
ggsave(filename = "HIV-PAPER/Figures/Extended_Figure 2/Extended_Fig2.pdf",
       plot = ((Extended_Fig2a|Extended_Fig2b|Extended_Fig2c|Extended_Fig2d|plot_layout(ncol = 4,nrow = 1, width = c(1,1,1,1))))/ 
         (Extended_Fig2e|Extended_Fig2f|Extended_Fig2g|Extended_Fig2h|plot_layout(width = c(1,1,1,1)))/
         (Extended_Fig2i|Extended_Fig2j|Extended_Fig2k|plot_layout(width = c(1,1,1,1))),
       width = 32, height = 29, units = "in", dpi = 300)

ggsave(filename = "HIV-PAPER/Figures/Extended_Figure 2/Extended_Fig2.png",
       plot = ((Extended_Fig2a|Extended_Fig2b|Extended_Fig2c|Extended_Fig2d|plot_layout(ncol = 4,nrow = 1, width = c(1,1,1,1))))/ 
         (Extended_Fig2e|Extended_Fig2f|Extended_Fig2g|Extended_Fig2h|plot_layout(width = c(1,1,1,1)))/
         (Extended_Fig2i|Extended_Fig2j|Extended_Fig2k|plot_layout(width = c(1,1,1,1))),
       width = 32, height = 29, units = "in", dpi = 300)


#################################### FIGURE 2 ##################################
# Figure 2c
Neut_CD11b_data <- read_csv("Data/Main_Files_Thesis/Neut_CD11b_data.csv")
Fig2b <- Neut_CD11b_data %>%
  ggplot(aes(x=factor(`HIV Status`,levels=c('PLHIV ART <3m', 'PLHIV ART >1y','HIV-')),
             y=as.numeric(`Neut CD11b+CD63-`),
             color=factor(`HIV Status`,
                          levels=c('PLHIV ART <3m', 'PLHIV ART >1y','HIV-'))))+
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
  labs(x='',
       y=expression("Proportion of CD11b"^"++"~"Neutrophils"),
       title = "b.")+
  scale_x_discrete(labels=c("HIV-"="PLHIV\nART<3m",
                            "PLHIV ART <3m"="HIV-",
                            "PLHIV ART >1y"="PLHIV\nART>1yr"))+
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

Fig2b
ggsave(Fig2b,filename="HIV-PAPER/Figures/Figure 2/Fig2b.png",
       width = 8,height = 10,dpi = 300)
ggsave(Fig2b,filename="HIV-PAPER/Figures/Figure 2/Fig2b.pdf",
       width = 8,height = 10,dpi = 300)

Fig2c_counts <- Masterfile %>%
  dplyr::filter(`HIV Status`!="NA") %>%
  dplyr::select(`Myeloperoxidase (pg/mL)`,`HIV Status`) %>%
  dplyr::group_by(`HIV Status`) %>%
  na.omit() %>%
  dplyr::summarize(n = dplyr::n(), .groups = "drop")
Fig2c_counts

Fig2c <- Masterfile %>%
  filter(`HIV Status`!="NA") %>%
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
ggsave(Fig2c,filename="HIV-PAPER/Figures/Figure 2/Fig2c.png",
       width = 8,height = 10,dpi = 300)
ggsave(Fig2c,filename="HIV-PAPER/Figures/Figure 2/Fig2c.pdf",
       width = 8,height = 10,dpi = 300)

# Figure 2d
Fig2d <- Masterfile %>%
  filter(`Epithelial count`>100,
         `CD45+ count`>200,
         `HIV Status`!="NA") %>%
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
ggsave(Fig2d,filename="HIV-PAPER/Figures/Figure 2/Fig2d.png",
       width = 8,height = 10,dpi = 300)
ggsave(Fig2d,filename="HIV-PAPER/Figures/Figure 2/Fig2d.pdf",
       width = 8,height = 10,dpi = 300)

# Figure 2e
Fig2e <- Masterfile %>%
  filter(`Epithelial count`>100,
         `CD45+ count`>200,
         `HIV Status`!="NA") %>%
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
ggsave(Fig2e,filename="HIV-PAPER/Figures/Figure 2/Fig2e.png",
       width = 8,height = 10,dpi = 300)
ggsave(Fig2e,filename="HIV-PAPER/Figures/Figure 2/Fig2e.pdf",
       width = 8,height = 10,dpi = 300)

# Saving Figure 3 using ggsave
Fig2a <- magick::image_read("New_Figures/Figure 3/Fig3a_new.png")
Fig2a <- ggdraw()+ 
  draw_image(Fig2a, scale = 1.2)+
  labs(title = "a.")+
  theme(
    plot.title = element_text(hjust = -0.1, vjust = 0.1, size = 50, face = "bold"),  # Customize title
    plot.margin = unit(c(0,0,0,0), "cm")  # Optional: adjust margins
  )
Fig2a


Fig2b <- magick::image_read("New_Figures/Figure 3/Fig3b.png")
Fig2b <- ggdraw()+ 
  draw_image(Fig3b, scale = 1.0)+
  #labs(title = "b.")+
  theme(
    plot.title = element_text(hjust = 0.05,vjust = 1, size = 50, face = "bold"),  # Customize title
    plot.margin = unit(c(0,0,0,0), "cm")  # Optional: adjust margins
  )
Fig3b

Fig3c <- magick::image_read("New_Figures/Figure 3/Fig3c.png")
Fig3c <- ggdraw()+ 
  draw_image(Fig3c, scale = 1.0)+
  #labs(title = "c.")+
  theme(
    plot.title = element_text(hjust = 0.05,vjust = 1, size = 50, face = "bold"),  # Customize title
    plot.margin = unit(c(0,0,0,0), "cm")  # Optional: adjust margins
  )
Fig3c

Fig3d <- magick::image_read("New_Figures/Figure 3/Fig3d.png")
Fig3d <- ggdraw()+ 
  draw_image(Fig3d, scale = 1.0)+
  #labs(title = "e.")+
  theme(
    plot.title = element_text(hjust = 0.05,vjust = 0, size = 50, face = "bold"),  # Customize title
    plot.margin = unit(c(0,0,0,0), "cm")  # Optional: adjust margins
  )
Fig3d

Fig3e <- magick::image_read("New_Figures/Figure 3/Fig3e.png")
Fig3e <- ggdraw()+ 
  draw_image(Fig3e, scale = 1.0)+
  #labs(title = "e.")+
  theme(
    plot.title = element_text(hjust = 0.05,vjust = 0, size = 50, face = "bold"),  # Customize title
    plot.margin = unit(c(0,0,0,0), "cm")  # Optional: adjust margins
  )
Fig3e

# Save figure 3 ggsave
ggsave(filename = "New_Figures/Figure 3/Fig3.png",
       plot = (plot_spacer()/
                 (plot_spacer()|Fig3a|plot_spacer()|plot_layout(widths = c(0.5,5,0.5)))/
                 (Fig3b|Fig3c|Fig3d|Fig3e|plot_layout(widths = c(1,1,1,1))) +
                 plot_layout(nrow = 3,heights = c(0.001,2,1))),
       width = 31, height = 29, units = "in", dpi = 300)

ggsave(filename = "New_Figures/Figure 3/Fig3.pdf",
       plot = (plot_spacer()/
                 (plot_spacer()|Fig3a|plot_spacer()|plot_layout(widths = c(0.5,5,0.5)))/
                 (Fig3b|Fig3c|Fig3d|Fig3e|plot_layout(widths = c(1,1,1,1))) +
                 plot_layout(nrow = 3,heights = c(0.001,2,1))),
       width = 31, height = 29, units = "in", dpi = 300)


#################################### FIGURE 3 ##################################
all_merged_subset_labelled_new <- readRDS("Data/Single_Cell_Data/all_merged_subset_labelled_new.rds")
samples <- unique(all_merged_subset_labelled_new@meta.data$sample)
clinical_metadata <- Masterfile %>%
  dplyr::filter(`LAB ID` %in% samples) %>%
  dplyr::mutate(`Possession index`=watch+radio+bank+iron+sew+mobile+
                  cd+fanelec+netyn+mattress+bed+bike+motor+car+tv) %>%
  dplyr::mutate(U5 = ifelse(grepl("1",Under_5),"One","2+")) %>%
  dplyr::select(PID,`LAB ID`,`Age`,`Sex`,`HIV Status`,
                `Carriage Status`, U5, `HIV duration`,
                `ART duration`, Location,`Carriage Status`,
                Density,Abs_CD4_Count,HIV_VL,`Possession index`)

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
HIV_status_colors <- c('#A9A9A9','#941100','#005493')
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
load("Data/Single_Cell_Data/Immune_cells.RData")
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
  #NoLegend()+
  labs(x="UMAP-1",y="UMAP-2"#,title = "d."
  )+
  scale_color_manual(values = base_palette)+
  theme_bw()+
  theme(axis.line = element_blank(),
        #legend.position = c(),
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
HIV_status_colors <- c('#A9A9A9','#941100','#005493')
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

# T cell to goblet cell ratio
Fig3g <- cluster_df %>%
  ggplot(aes(factor(`HIV Status`,
                    levels = c("HIV-",
                               "ART<3m",
                               "ART>1y")),
             log10(`Tcell:Goblet`),
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
       y=expression("T cell to epithelial cell ratio"~"(log"[10]~")"),
       fill='HIV Status'#,title = "g."
  )+
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
Fig3g
ggsave(Fig3g,filename="HIV-PAPER/Figures/Figure 3/Fig3g.png",
       width = 14,height = 18,dpi = 300)
ggsave(Fig3g,filename="HIV-PAPER/Figures/Figure 3/Fig3g.png",
       width = 14,height = 18,dpi = 300)

# Save figure 3 in ggsave
Fig3a <- magick::image_read("HIV-PAPER/Figures/Figure 3/Fig3a.png")
Fig3a <- ggdraw()+ 
  draw_image(Fig3a, scale = 1.3)+
  labs(title = "a.")+
  theme(
    plot.title = element_text(hjust = 0.01,vjust = .01, size = 60, face = "bold"),  # Customize title
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")  # Optional: adjust margins
  )
Fig3a

Fig3c <- magick::image_read("HIV-PAPER/Figures/Figure 3/Fig3c.png")
Fig3c <- ggdraw()+ 
  draw_image(Fig3c, scale = 1)+
  labs(title = "c.")+
  theme(
    plot.title = element_text(hjust = 0.01,vjust = .01, size = 60, face = "bold"),  # Customize title
    plot.margin = unit(c(0,0,0,0), "cm")  # Optional: adjust margins
  )
Fig3c

Fig3d <- magick::image_read("HIV-PAPER/Figures/Figure 3/Fig3d.png")
Fig3d <- ggdraw()+ 
  draw_image(Fig3d, scale = 1)+
  labs(title = "d.")+
  theme(
    plot.title = element_text(hjust = 0.01,vjust = .01, size = 60, face = "bold"),  # Customize title
    plot.margin = unit(c(0,0,0,0), "cm")  # Optional: adjust margins
  )
Fig3d

Fig3g <- magick::image_read("HIV-PAPER/Figures/Figure 3/Fig3g.png")
Fig3g <- ggdraw()+ 
  draw_image(Fig3g, scale = 1)+
  labs(title = "g.")+
  theme(
    plot.title = element_text(hjust = 0.01,vjust = .01, size = 60, face = "bold"),  # Customize title
    plot.margin = unit(c(0,0,0,0), "cm")  # Optional: adjust margins
  )
Fig3g
# Save figure 3 ggsave
ggsave(filename = "HIV-PAPER/Figures/Figure 3/Fig3.pdf",
       plot = ((Fig3a|Fig3b|plot_layout(ncol = 2,nrow = 1, width = c(1.5,0.9), heights = c(1.3,1))))/ 
         (Fig3c|plot_layout(width = c(1)))/
         (Fig3d|Fig3e|Fig3f|plot_layout(width = c(1,1,1))),
       width = 35, height = 37, units = "in", dpi = 300)

ggsave(filename = "HIV-PAPER/Figures/Figure 3/Fig3.png",
       plot = ((Fig3a|Fig3b|plot_layout(ncol = 2,nrow = 1, width = c(1.5,0.9), heights = c(1.3,1))))/ 
         (Fig3c|plot_layout(width = c(1)))/
         (Fig3d|Fig3e|Fig3f|plot_layout(width = c(1,1,1))),
       width = 35, height = 37, units = "in", dpi = 300)


#################################### EXTENDED FIGURE 3 #########################
# Extended Figure 3a (vlnPlot of celltype viable features)
all_merged_subset_labelled_new <- readRDS("Data/Single_Cell_Data/all_merged_subset_labelled_new.rds")
all_merged_subset_labelled_new <- FindVariableFeatures(all_merged_subset_labelled_new,
                                                       selection.method = 'vst',
                                                       verbose = TRUE)
genes <- VariableFeatures(all_merged_subset_labelled_new)
toplot <- SeuratExtend::CalcStats(all_merged_subset_labelled_new,
                                  features = genes,
                                  method = 'zscore',
                                  order = 'p',
                                  n = 5)

Extended_Fig3a <- Seurat::VlnPlot(all_merged_subset_labelled_new,
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
Extended_Fig3a

# Save Extended_Figure 3a
ggsave(Extended_Fig3a,filename="HIV-PAPER/Figures/Extended_Figure 3/Extended_Fig3.png",
       width = 20,height = 6,dpi = 300,units = "in")
ggsave(Extended_Fig3a,filename="HIV-PAPER/Figures/Extended_Figure 3/Extended_Fig3.pdf",
       width = 20,height = 6,dpi = 300,units = "in")

#################################### FIGURE 4 ##################################
# Figure 4a (Epithelial cells communicating with Neutrophils)
filt_prioritized_tbl_oi <- readRDS("scRNAseq_Results/filt_prioritized_tbl_oi.rds")
senders_receivers = union(filt_prioritized_tbl_oi$sender %>% unique(), filt_prioritized_tbl_oi$receiver %>% unique()) %>% sort()

colors_sender = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)
colors_receiver = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)
Neut_circos_list = multinichenetr::make_circos_group_comparison(filt_prioritized_tbl_oi, colors_sender, colors_receiver)
Neg <- Neut_circos_list$HIV.
`3M` <- Neut_circos_list$HIV..ART.3.Months
`1Y` <- Neut_circos_list$HIV..ART.1.Year

library(gridGraphics)
recordedplot_to_grob <- function(recorded_plot) {
  gridGraphics::grid.echo(recorded_plot)
  grid::grid.grab()
}

Neg <- recordedplot_to_grob(Neg)
`3M` <- recordedplot_to_grob(`3M`)
`1Y` <- recordedplot_to_grob(`1Y`)

library(ggplot2)

# Example: Convert grob to ggplot
ggplot_from_grob <- function(grob) {
  ggplot() +
    annotation_custom(grob) +
    theme_void()  # Remove axes, grid lines, etc.
}

# Convert your grob
Neg <- ggplot_from_grob(Neg)
`3M` <- ggplot_from_grob(`3M`)
`1Y` <- ggplot_from_grob(`1Y`)
legend <- Neut_circos_list$legend

# Save Fig4a
ggsave(Neg,filename="HIV-PAPER/Figures/Figure 4/Fig4a1.png",
       width = 6,height = 10,dpi = 300)
ggsave(`3M`,filename="HIV-PAPER/Figures/Figure 4/Fig4a2.png",
       width = 6,height = 10,dpi = 300)
ggsave(`1Y`,filename="HIV-PAPER/Figures/Figure 4/Fig4a3.png",
       width = 6,height = 10,dpi = 300)
png("HIV-PAPER/Figures/Figure 4/Fig4a4.png", width = 800, height = 800)
replayPlot(legend)
dev.off()



# Figure 4g
all_merged_subset_labelled_new <- readRDS("Data/Single_Cell_Data/all_merged_subset_labelled_new.rds")
all_merged_subset_labelled_new$Clusters <- paste0(all_merged_subset_labelled_new@active.ident)
Neutrophil_recruitment <- list(
  CHEMOTAXIS_GENES=c("ANXA1","BST2","CD55","CDH1","CXCL1","MFAP2","MIF","TSPAN3","C3",
                     "ADGRE5","ICAM1","IL1B","S100A8","TFF1","ADAM10","CD44","DYSF","FPR1","IGF1R","IL1R1",
                     "IL2RG","ITGAX","ITGB2","LILRA5","MSN","MUC5AC","NOTCH1"
  ))

all_merged_subset_labelled_new <- Seurat::AddModuleScore(
  all_merged_subset_labelled_new, features = list(c(Neutrophil_recruitment$CHEMOTAXIS_GENES)),
  name = "Chemotaxis_Score",assay = "RNA"
)

Fig4b <- all_merged_subset_labelled_new@meta.data %>%
  dplyr::filter(Clusters  %in% c('Basal cells')
  ) %>%
  ggplot(aes(HIV_Status,Chemotaxis_Score1, fill=HIV_Status))+
  geom_violin(trim = T, scale = "width")+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=43,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 'wilcox.test',
           label = "{paste0('p = ', signif(p, digits = 4))}", 
           #label = "{ifelse(p < 0.0001, 'p < 0.0001', sprintf('p = %.4f', p))}",
           label.size = 10,
           tip.length = 0.01,
           step.increase = 0.2,
           p.adjust.method = "bonferroni",
           hide.ns = F)+
  labs(x='',
       y='Gene Module Expression',
       #title = "b.",
       subtitle = "Basal cells")+
  scale_fill_manual(values = c('#A9A9A9','#941100','#005493'))+
  scale_x_discrete(labels=c("HIV+ ART<3 Months"="PLHIV\nART<3m",
                            "HIV+ ART>1 Year"="PLHIV\nART>1yr",
                            "HIV-"="HIV-"))+
  scale_y_continuous(limits = c(-.4,1.25),
                     breaks = c(-0.4,0,0.4,0.8,1.2))+
  theme_classic()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        legend.text = element_text(size=25), 
        axis.ticks.length = unit(0.5,'cm'),
        plot.title = element_text(size = 50, face = 'bold'),
        plot.subtitle = element_text(size = 30,hjust = 0.5),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 30))
Fig4b
ggsave(Fig4b,filename="HIV-PAPER/Figures/Figure 4/Fig4b.png",
       width = 8,height = 10,dpi = 300)
ggsave(Fig4b,filename="HIV-PAPER/Figures/Figure 4/Fig4b.pdf",
       width = 8,height = 10,dpi = 300)

Fig4c <- all_merged_subset_labelled_new@meta.data %>%
  dplyr::filter(Clusters  %in% c('Ciliated cells')
  ) %>%
  ggplot(aes(HIV_Status,Chemotaxis_Score1, fill=HIV_Status))+
  geom_violin(trim = T, scale = "width")+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=43,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 'wilcox.test',
           label = "{paste0('p = ', signif(p, digits = 4))}", 
           #label = "{ifelse(sprintf('p = %.4f', p))}",
           label.size = 10,
           tip.length = 0.01,
           step.increase = 0.2,
           p.adjust.method = "holm",
           hide.ns = F)+
  labs(x='',
       y='Gene Module Expression',
       #title = "b.",
       subtitle = "Ciliated cells")+
  scale_fill_manual(values = c('#A9A9A9','#941100','#005493'))+
  scale_x_discrete(labels=c("HIV+ ART<3 Months"="PLHIV\nART<3m",
                            "HIV+ ART>1 Year"="PLHIV\nART>1yr",
                            "HIV-"="HIV-"))+
  scale_y_continuous(limits = c(-.4,1.25),
                     breaks = c(-0.4,0,0.4,0.8,1.2))+
  theme_classic()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        legend.text = element_text(size=25), 
        axis.ticks.length = unit(0.5,'cm'),
        plot.title = element_text(size = 50, face = 'bold'),
        plot.subtitle = element_text(size = 30,hjust = 0.5),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 30))
Fig4c
ggsave(Fig4c,filename="HIV-PAPER/Figures/Figure 4/Fig4c.png",
       width = 8,height = 10,dpi = 300)
ggsave(Fig4c,filename="HIV-PAPER/Figures/Figure 4/Fig4c.pdf",
       width = 8,height = 10,dpi = 300)

Fig4d <- all_merged_subset_labelled_new@meta.data %>%
  dplyr::filter(Clusters  %in% c('Goblet cells')
  ) %>%
  ggplot(aes(HIV_Status,Chemotaxis_Score1, fill=HIV_Status))+
  geom_violin(trim = T, scale = "width")+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=43,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 'wilcox.test',
           label = "{paste0('p = ', signif(p, digits = 4))}", 
           #label = "{ifelse(p < 0.0001, 'p < 0.0001', sprintf('p = %.4f', p))}",
           label.size = 10,
           tip.length = 0.01,
           step.increase = 0.2,
           p.adjust.method = "holm",
           hide.ns = F)+
  labs(x='',
       y='Gene Module Expression',
       #title = "b.",
       subtitle = "Goblet cells")+
  scale_fill_manual(values = c('#A9A9A9','#941100','#005493'))+
  scale_x_discrete(labels=c("HIV+ ART<3 Months"="PLHIV\nART<3m",
                            "HIV+ ART>1 Year"="PLHIV\nART>1yr",
                            "HIV-"="HIV-"))+
  scale_y_continuous(limits = c(-.4,1.25),
                     breaks = c(-0.4,0,0.4,0.8,1.2))+
  theme_classic()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        legend.text = element_text(size=25), 
        axis.ticks.length = unit(0.5,'cm'),
        plot.title = element_text(size = 50, face = 'bold'),
        plot.subtitle = element_text(size = 30,hjust = 0.5),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 30))
Fig4d
ggsave(Fig4d,filename="HIV-PAPER/Figures/Figure 4/Fig4d.png",
       width = 8,height = 10,dpi = 300)
ggsave(Fig4d,filename="HIV-PAPER/Figures/Figure 4/Fig4d.pdf",
       width = 8,height = 10,dpi = 300)


Fig4e <- all_merged_subset_labelled_new@meta.data %>%
  dplyr::filter(Clusters  %in% c('Squamous cells')
  ) %>%
  ggplot(aes(HIV_Status,Chemotaxis_Score1, fill=HIV_Status))+
  geom_violin(trim = T, scale = "width")+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=43,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 'wilcox.test',
           label = "{paste0('p = ', signif(p, digits = 4))}", 
           #label = "{ifelse(p < 0.0001, 'p < 0.0001', sprintf('p = %.4f', p))}",
           label.size = 10,
           tip.length = 0.01,
           step.increase = 0.2,
           p.adjust.method = "holm",
           hide.ns = F)+
  labs(x='',
       y='Gene Module Expression',
       #title = "b.",
       subtitle = "Squamous cells")+
  scale_fill_manual(values = c('#A9A9A9','#941100','#005493'))+
  scale_x_discrete(labels=c("HIV+ ART<3 Months"="PLHIV\nART<3m",
                            "HIV+ ART>1 Year"="PLHIV\nART>1yr",
                            "HIV-"="HIV-"))+
  scale_y_continuous(limits = c(-.2,1.25),
                     breaks = c(0,0.4,0.8,1.2))+
  theme_classic()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        legend.text = element_text(size=25), 
        axis.ticks.length = unit(0.5,'cm'),
        plot.title = element_text(size = 50, face = 'bold'),
        plot.subtitle = element_text(size = 30,hjust = 0.5),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 30))
Fig4e
ggsave(Fig4e,filename="HIV-PAPER/Figures/Figure 4/Fig4e.png",
       width = 8,height = 10,dpi = 300)
ggsave(Fig4e,filename="HIV-PAPER/Figures/Figure 4/Fig4e.pdf",
       width = 8,height = 10,dpi = 300)


Fig4f <- all_merged_subset_labelled_new@meta.data %>%
  dplyr::filter(Clusters  %in% c('Secretory cells')
  ) %>%
  ggplot(aes(HIV_Status,Chemotaxis_Score1, fill=HIV_Status))+
  geom_violin(trim = T, scale = "width")+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=43,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 'wilcox.test',
           #label = "{paste0('p = ', signif(p, digits = 4))}", 
           label = "{ifelse(p < 0.0001, 'p = 0.0001', sprintf('p = %.4f', p))}",
           label.size = 10,
           tip.length = 0.01,
           step.increase = 0.2,
           p.adjust.method = "holm",
           hide.ns = F)+
  labs(x='',
       y='Gene Module Expression',
       #title = "b.",
       subtitle = "Secretory cells")+
  scale_fill_manual(values = c('#A9A9A9','#941100','#005493'))+
  scale_x_discrete(labels=c("HIV+ ART<3 Months"="PLHIV\nART<3m",
                            "HIV+ ART>1 Year"="PLHIV\nART>1yr",
                            "HIV-"="HIV-"))+
  scale_y_continuous(limits = c(-.3,1.25),
                     breaks = c(0,0.4,0.8,1.2))+
  theme_classic()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        legend.text = element_text(size=25), 
        axis.ticks.length = unit(0.5,'cm'),
        plot.title = element_text(size = 50, face = 'bold'),
        plot.subtitle = element_text(size = 30,hjust = 0.5),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 30))
Fig4f
ggsave(Fig4f,filename="HIV-PAPER/Figures/Figure 4/Fig4f.png",
       width = 8,height = 10,dpi = 300)
ggsave(Fig4f,filename="HIV-PAPER/Figures/Figure 4/Fig4f.pdf",
       width = 8,height = 10,dpi = 300)


# Save figure 4 ggsave
Fig4a <- magick::image_read("HIV-PAPER/Figures/Figure 4/Fig4a.png")
Fig4a <- ggdraw()+ 
  draw_image(Fig4a, scale = 1.2)+
  labs(title = "a.")+
  theme(
    plot.title = element_text(hjust = 0.0, vjust = .01, size = 50, face = "bold"),  # Customize title
    plot.margin = unit(c(0,0,0,0), "cm")  # Optional: adjust margins
  )
Fig4a

Fig4b <- magick::image_read("HIV-PAPER/Figures/Figure 4/Fig4b.png")
Fig4b <- ggdraw()+ 
  draw_image(Fig4b, scale = 1)+
  labs(title = "b.")+
  theme(
    plot.title = element_text(hjust = 0.05, vjust = .01, size = 50, face = "bold"),  # Customize title
    plot.margin = unit(c(0,0,0,0), "cm")  # Optional: adjust margins
  )
Fig4b

Fig4c <- magick::image_read("HIV-PAPER/Figures/Figure 4/Fig4c.png")
Fig4c <- ggdraw()+ 
  draw_image(Fig4c, scale = 1)+
  labs(title = "c.")+
  theme(
    plot.title = element_text(hjust = 0.05, vjust = .01, size = 50, face = "bold"),  # Customize title
    plot.margin = unit(c(0,0,0,0), "cm")  # Optional: adjust margins
  )
Fig4c

Fig4d <- magick::image_read("HIV-PAPER/Figures/Figure 4/Fig4d.png")
Fig4d <- ggdraw()+ 
  draw_image(Fig4d, scale = 1)+
  labs(title = "d.")+
  theme(
    plot.title = element_text(hjust = 0.05, vjust = .01, size = 50, face = "bold"),  # Customize title
    plot.margin = unit(c(0,0,0,0), "cm")  # Optional: adjust margins
  )
Fig4d

Fig4e <- magick::image_read("HIV-PAPER/Figures/Figure 4/Fig4e.png")
Fig4e <- ggdraw()+ 
  draw_image(Fig4e, scale = 1)+
  labs(title = "e.")+
  theme(
    plot.title = element_text(hjust = 0.05, vjust = .01, size = 50, face = "bold"),  # Customize title
    plot.margin = unit(c(0,0,0,0), "cm")  # Optional: adjust margins
  )
Fig4e

Fig4f <- magick::image_read("HIV-PAPER/Figures/Figure 4/Fig4f.png")
Fig4f <- ggdraw()+ 
  draw_image(Fig4f, scale = 1)+
  labs(title = "f.")+
  theme(
    plot.title = element_text(hjust = 0.05, vjust = .01, size = 50, face = "bold"),  # Customize title
    plot.margin = unit(c(0,0,0,0), "cm")  # Optional: adjust margins
  )
Fig4f

ggsave(filename = "HIV-PAPER/Figures/Figure 4/Fig4.png",
       plot = ((Fig4a|plot_spacer()|plot_layout(ncol = 2,nrow = 1, width = c(1,0.0001))))/ 
         (Fig4b|Fig4c|Fig4d|Fig4e|plot_layout(width = c(1,1,1,1)))/
         (Fig4f|plot_layout(width = c(1,1,1,1))),
       width = 32, height = 29, units = "in", dpi = 300)

ggsave(filename = "HIV-PAPER/Figures/Figure 4/Fig4.pdf",
       plot = ((Fig4a|plot_spacer()|plot_layout(ncol = 2,nrow = 1, width = c(1,0.0001))))/ 
         (Fig4b|Fig4c|Fig4d|Fig4e|plot_layout(width = c(1,1,1,1)))/
         (Fig4f|plot_layout(width = c(1,1,1,1))),
       width = 32, height = 29, units = "in", dpi = 300)


#################################### EXTENDED_FIGURE 4 #########################
# Extended Figure 4a (Immune cells communicating with Neutrophils)
Immune_multinichenet_output <- readRDS("scRNAseq_Results/Immune_multinichenet_output.rds")
# Visualization of differential cell-cell interactions
Immune_prioritized_tbl_oi_all = multinichenetr::get_top_n_lr_pairs(
  Immune_multinichenet_output$prioritization_tables, 
  top_n = 30, 
  rank_per_group = T)

Immune_prioritized_tbl_oi = 
  Immune_multinichenet_output$prioritization_tables$group_prioritization_tbl %>%
  filter(id %in% Immune_prioritized_tbl_oi_all$id) %>%
  distinct(id, sender, receiver, ligand, receptor, group) %>% 
  left_join(Immune_prioritized_tbl_oi_all)
Immune_prioritized_tbl_oi$prioritization_score[is.na(Immune_prioritized_tbl_oi$prioritization_score)] = 0


senders_receivers = union(Immune_prioritized_tbl_oi$sender %>% unique(), Immune_prioritized_tbl_oi$receiver %>% unique()) %>% sort()

colors_sender = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)
colors_receiver = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)
Immune_circos_list = multinichenetr::make_circos_group_comparison(Immune_prioritized_tbl_oi, colors_sender, colors_receiver)
Immune_Neg <- Immune_circos_list$HIV.
Immune_3m <- Immune_circos_list$HIV..ART.3.Months
Immune_1y <- Immune_circos_list$HIV..ART.1.Year

library(gridGraphics)
recordedplot_to_grob <- function(recorded_plot) {
  gridGraphics::grid.echo(recorded_plot)
  grid::grid.grab()
}

Immune_Neg <- recordedplot_to_grob(Immune_Neg)
Immune_3m <- recordedplot_to_grob(Immune_3m)
Immune_1y <- recordedplot_to_grob(Immune_1y)

library(ggplot2)

# Example: Convert grob to ggplot
ggplot_from_grob <- function(grob) {
  ggplot() +
    annotation_custom(grob) +
    theme_void()  # Remove axes, grid lines, etc.
}

# Convert your grob
Immune_Neg <- ggplot_from_grob(Immune_Neg)
Immune_3m <- ggplot_from_grob(Immune_3m)
Immune_1y <- ggplot_from_grob(Immune_1y)
Immune_legend <- Immune_circos_list$legend

# Save Extended Fig5a
ggsave(Immune_Neg,filename="HIV-PAPER/Figures/Extended_Figure 4/Extended_Fig4a1.png",
       width = 12,height = 19,dpi = 300)
ggsave(Immune_3m,filename="HIV-PAPER/Figures/Extended_Figure 4/Extended_Fig4a2.png",
       width = 12,height = 19,dpi = 300)
ggsave(Immune_1y,filename="HIV-PAPER/Figures/Extended_Figure 4/Extended_Fig4a3.png",
       width = 12,height = 19,dpi = 300)
png("HIV-PAPER/Figures/Extended_Figure 4/Extended_Fig4a4.png", width = 800, height = 800)
replayPlot(Immune_legend)
dev.off()


Extended_Fig5a <- magick::image_read("New_Figures/Extended_Figure 5/Extended_Fig5a.png")
Extended_Fig5a <- ggdraw()+ 
  draw_image(Extended_Fig5a, scale = 1.0)+
  labs(title = "a.")+
  theme(
    plot.title = element_text(hjust = 0.0, vjust = .01, size = 50, face = "bold"),  # Customize title
    plot.margin = unit(c(0,0,0,0), "cm")  # Optional: adjust margins
  )
Extended_Fig5a

ggsave(filename = "New_Figures/Extended_Figure 5/Extended_Fig5.pdf",
       plot = ((Extended_Fig5a|plot_spacer()|plot_layout(ncol = 2,nrow = 1, width = c(1,0.0001))))/ 
         (plot_spacer()|plot_spacer()|plot_spacer()|plot_spacer()|plot_layout(width = c(1,1,1,1)))/
         (plot_spacer()|plot_layout(width = c(1,1,1,1))),
       width = 32, height = 29, units = "in", dpi = 300)

ggsave(filename = "New_Figures/Extended_Figure 5/Extended_Fig5.png",
       plot = ((Extended_Fig5a|plot_spacer()|plot_layout(ncol = 2,nrow = 1, width = c(1,0.0001))))/ 
         (plot_spacer()|plot_spacer()|plot_spacer()|plot_spacer()|plot_layout(width = c(1,1,1,1)))/
         (plot_spacer()|plot_layout(width = c(1,1,1,1))),
       width = 32, height = 29, units = "in", dpi = 300)
````

#################################### FIGURE 5 ##################################
# Differential gene expression using MAST

# Load libraries
library(Seurat)
library(dplyr)

load("Data/Single_Cell_Data/Immune_cells.RData")
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
Neutrophils <- subset(Immune_cells, 
                      idents = "Neutrophils",
                      invert=F)
# Add all module scores in one call
Immune_cells <- AddModuleScore(
  object = Immune_cells,
  features = module_features,
  name = names(module_features),  # Names will be appended with numbers (e.g., "Glycan_biosynthesis1")
  slot = "data"
)


Fig5g <- Immune_cells@meta.data %>%
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


Fig5h <- Immune_cells@meta.data %>%
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


Fig5i <- Immune_cells@meta.data %>%
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
Fig5j <- Immune_cells@meta.data %>%
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
# MFI of Phagocytosis
Fig5k_counts <- Phagocytosis_data %>%
  dplyr::filter(Time=="45min",
                #Visit=="Week 1",
                Stimulant=="LPS") %>%
  dplyr::mutate(Oxidation_SI = (`MFI+`-`MFI-`)/(2*`rSD-`)) %>%
  dplyr::mutate(Product_Phagocytosis_MF1_Frequency = `MFI phagocytosis`*`Proportion of Phagocytosis`) %>%
  dplyr::select(Product_Phagocytosis_MF1_Frequency,`HIV Status`) %>%
  dplyr::group_by(`HIV Status`) %>%
  dplyr::summarize(n = dplyr::n(), .groups = "drop")
Fig5k_counts


Fig5k <- Phagocytosis_data %>%
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
Fig5l_counts <- Phagocytosis_data %>%
  dplyr::filter(Time=="45min",
                #Visit=="Week 1",
                Stimulant=="LPS") %>%
  dplyr::mutate(Oxidation_SI = (`MFI+`-`MFI-`)/(2*`rSD-`)) %>%
  dplyr::mutate(Product_Phagocytosis_MF1_Frequency = `MFI phagocytosis`*`Proportion of Phagocytosis`) %>%
  dplyr::select(Oxidation_SI,`HIV Status`) %>%
  dplyr::group_by(`HIV Status`) %>%
  dplyr::summarize(n = dplyr::n(), .groups = "drop")
Fig5l_counts

Fig5l <- Phagocytosis_data %>%
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



#################################### EXTENDED FIGURE 5 #########################
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
all_merged_subset_labelled_new <- readRDS("Data/Single_Cell_Data/all_merged_subset_labelled_new.rds")
CD3_T_cells <- subset(all_merged_subset_labelled_new, idents = c("CD3+ T cells"))
Idents(CD3_T_cells)<-CD3_T_cells$HIV_Status
load("data/Single_Cell_data/Immune_cells.RData")
Tcells <- subset(Immune_cells, subset = Clusters == c("CD8+ T cells","NK T cells"))
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


all_merged_subset_labelled_new <- readRDS('Data/Single_Cell_Data/all_merged_subset_labelled_new.rds')
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



#################################### EXTENDED FIGURE 6 #########################
Extended_Fig6a <- Tcell_cem@ora %>%
  dplyr::filter(Module=="M1") %>%
  dplyr::mutate(Significance=ifelse(p.adjust<.05 & FoldEnrichment>0,"Upregulated",
                                    ifelse(p.adjust<.05 & FoldEnrichment<0,"Downregulated","Not significant")),
                label = ifelse(Significance=="Upregulated",Description,NA)) %>%
  ggplot(aes(FoldEnrichment,-log10(p.adjust), color=Significance))+
  geom_point(size=2.5,alpha=.5)+
  geom_hline(yintercept = -log10(0.05), linetype="dashed")+
  geom_text_repel(aes(label = label),size=5,color='black')+
  scale_color_manual(values = c('grey','red'))+
  labs(x="Fold-Enrichment",
       y="-Log"[10]~"(padj)",
       title = 'a.',
       subtitle = "M1 T cell enriched modules\nPLHIV-ART>1yr")+
  theme_classic()+
  theme(plot.title = element_text(size = 50, face = "bold"),
        legend.position = 'none',
        axis.ticks.length = unit(0.5,'cm'),
        plot.subtitle = element_text(size = 30, hjust = .5),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 35))
Extended_Fig6a

ggsave(Extended_Fig6a,filename="HIV-PAPER/Figures/Extended_Figure 6/Extended_Fig6a.png",
       width = 12,height = 15,dpi = 300,units = "in")
ggsave(Extended_Fig6a,filename="HIV-PAPER/Figures/Extended_Figure 6/Extended_Fig6a.pdf",
       width = 12,height = 15,dpi = 300,units = "in")

Extended_Fig6b <- Tcell_cem@ora %>%
  dplyr::filter(Module=="M2") %>%
  dplyr::mutate(Significance=ifelse(p.adjust<.05 & FoldEnrichment>0,"Upregulated",
                                    ifelse(p.adjust<.05 & FoldEnrichment<0,"Downregulated","Not significant")),
                label = ifelse(Significance=="Upregulated",Description,NA)) %>%
  ggplot(aes(FoldEnrichment,-log10(p.adjust), color=Significance))+
  geom_point(size=2.5,alpha=.5)+
  geom_hline(yintercept = -log10(0.05), linetype="dashed")+
  geom_text_repel(aes(label = label),size=5,color='black')+
  scale_color_manual(values = c('grey','red'))+
  labs(x="Fold-Enrichment",
       y="-Log"[10]~"(padj)",
       title = 'b.',
       subtitle = "M2 T cell enriched modules\n(PLHIV-ART<3m)")+
  theme_classic()+
  theme(plot.title = element_text(size = 50, face = "bold"),
        legend.position = 'none',
        axis.ticks.length = unit(0.5,'cm'),
        plot.subtitle = element_text(size = 30,hjust = .5),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 35))
Extended_Fig6b
# Save Figure 6b
ggsave(Extended_Fig6b,filename="HIV-PAPER/Figures/Extended_Figure 6/Extended_Fig6b.png",
       width = 12,height = 15,dpi = 300,units = "in")
ggsave(Extended_Fig6b,filename="HIV-PAPER/Figures/Extended_Figure 6/Extended_Fig6b.pdf",
       width = 12,height = 15,dpi = 300,units = "in")


# Save Extended Figure 6
ggsave(filename = "HIV-PAPER/Figures/Extended_Figure 6/Extended_Fig6.png",
       plot = ((Extended_Fig6a|Extended_Fig6b|plot_layout(ncol = 3,nrow = 1, width = c(1,1,1))))/
         (plot_spacer()|plot_layout(width = c(3)))/
         (plot_spacer()|plot_layout(width = c(3))),
       width = 32, height = 29, units = "in", dpi = 300,limitsize = F)

ggsave(filename = "HIV-PAPER/Figures/Extended_Figure 6/Extended_Fig6.pdf",
       plot = ((Extended_Fig6a|Extended_Fig6b|plot_layout(ncol = 3,nrow = 1, width = c(1,1,1))))/
         (plot_spacer()|plot_layout(width = c(3)))/
         (plot_spacer()|plot_layout(width = c(3))),
       width = 32, height = 29, units = "in", dpi = 300,limitsize = F)


#################################### FIGURE 7 ##################################
# Relationship between Nasopharyngeal pneumococcal carriage and Neutrophil dynamics
Fig7a_counts <- Masterfile %>%
  filter(`Epithelial count`>250,
         `CD45+ count`>500,
         `HIV Status`!="NA",
         `Visit` %in% c(#"Week 1"#,
                        #"Week 2",
                        #"Week 3",
                        #"Week 4",
                        "Week 5"
         )
  ) %>%
  dplyr::select(NER,`HIV Status`, `Carriage Status`) %>%
  dplyr::group_by(`HIV Status`, `Carriage Status`) %>%
  dplyr::summarize(n = dplyr::n(), .groups = "drop")
Fig7a_counts

Fig7a <- Masterfile %>%
  filter(`Epithelial count`>250,
         `CD45+ count`>500,
         `HIV Status`!="NA",
         `Visit` %in% c("Week 1",
                        #"Week 2",
                        #"Week 3",
                        #"Week 4",
                        "Week 5"
                        )
  ) %>%
  ggplot(aes(x=factor(`Carriage Status`,levels=c('Carriage Negative','Carriage Positive')),
             y=as.numeric(log10(`NER`)),
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
  #geom_boxplot(width=0.1,notch = T, outlier.shape = NA, fill='white')+
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
  #scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)), 
  #limits = c(10^-2.5,10^1.5))+
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

Fig7b_counts <- Masterfile %>%
  filter(`Epithelial count`>250,
         `CD45+ count`>500,
         `HIV Status`!="NA",
         `Visit` %in% c(#"Week 1"#,
           #"Week 2",
           #"Week 3",
           #"Week 4",
           "Week 5"
         )
  ) %>%
  dplyr::select(MER,`HIV Status`, `Carriage Status`) %>%
  dplyr::group_by(`HIV Status`, `Carriage Status`) %>%
  dplyr::summarize(n = dplyr::n(), .groups = "drop")
Fig7b_counts

Fig7b <- Masterfile %>%
  filter(`Epithelial count`>250,
         `CD45+ count`>500,
         `HIV Status`!="NA",
         `Visit` %in% c(
           "Week 1",
           #"Week 2"#,
           #"Week 3",
           #"Week 4",
           "Week 5"
         )
  ) %>%
  ggplot(aes(x=factor(`Carriage Status`,levels=c('Carriage Negative','Carriage Positive')),
             y=as.numeric(log10(`MER`)),
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
  #geom_boxplot(width=0.1,notch = T, outlier.shape = NA, fill='white')+
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
  #scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)), 
  #limits = c(10^-2.5,10^1.5))+
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


Fig7c_counts <- Masterfile %>%
  filter(`Epithelial count`>500,
         `CD45+ count`>200,
         Visit %in% c(
           "Week 1",
           "Week 2",
           "Week 3",
           "Week 4",
           "Week 5"
           ),
         `HIV Status`!="NA",
         `Neutrophil CD11b++`>15 & `Neutrophil CD11b++`<90,
         `Carriage Status`=="Carriage Positive",
         `CD11b Staining`=="Stained") %>%
  dplyr::select(`Neutrophil CD11b++`,`HIV Status`) %>%
  dplyr::group_by(`HIV Status`) %>%
  dplyr::summarize(n = dplyr::n(), .groups = "drop")
Fig7c_counts


Fig7c <- Masterfile %>%
  filter(`Epithelial count`>500,
         `CD45+ count`>200,
         Visit %in% c("Week 1","Week 2","Week 3","Week 4", "Week 5"),
         `HIV Status`!="NA",
         `Neutrophil CD11b++`>15 & `Neutrophil CD11b++`<90,
         `Carriage Status`=="Carriage Positive",
         `CD11b Staining`=="Stained") %>%
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


Fig7d_counts <- Masterfile %>%
  filter(`Epithelial count`>500,
         `CD45+ count`>200,
         Visit %in% c(
           "Week 1",
           "Week 2",
           "Week 3",
           "Week 4", 
           "Week 5"),
         `HIV Status`!="NA",
         #`Neutrophil CD11b++`>15 & `Neutrophil CD11b++`<90,
         `Carriage Status`=="Carriage Positive",
         `CD11b Staining`=="Stained") %>%
  dplyr::select(MER,`HIV Status`) %>%
  dplyr::group_by(`HIV Status`) %>%
  dplyr::summarize(n = dplyr::n(), .groups = "drop")
Fig7d_counts

Fig7d <- Masterfile %>%
  filter(`Epithelial count`>500,
         `CD45+ count`>200,
         Visit %in% c("Week 1","Week 2","Week 3","Week 4", "Week 5"),
         `HIV Status`!="NA",
         #`Neutrophil CD11b++`>15 & `Neutrophil CD11b++`<90,
         `Carriage Status`=="Carriage Positive",
         `CD11b Staining`=="Stained") %>%
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


Phagocytosis_data <- read_csv("Data/Main_Files_Thesis/Neutrophil_Phagocytosis.csv")
Fig7e_counts <- Phagocytosis_data %>%
  dplyr::filter(Time=="45min",
                #Visit=="Week 1",
                Stimulant=="LPS") %>%
  dplyr::mutate(Oxidation_SI = (`MFI+`-`MFI-`)/(2*`rSD-`)) %>%
  dplyr::mutate(Product_Phagocytosis_MF1_Frequency = `MFI phagocytosis`*`Proportion of Phagocytosis`) %>%
  dplyr::select(Product_Phagocytosis_MF1_Frequency,`HIV Status`, `Carriage Status`) %>%
  dplyr::group_by(`HIV Status`, `Carriage Status`) %>%
  dplyr::summarize(n = dplyr::n(), .groups = "drop")
Fig7e_counts

# MFI of Oxidation
Fig7e <- Phagocytosis_data %>%
  dplyr::filter(Time=="45min",
                #Visit=="Week 1",
                Stimulant=="LPS") %>%
  dplyr::mutate(Oxidation_SI = (`MFI+`-`MFI-`)/(2*`rSD-`)) %>%
  dplyr::mutate(Product_Phagocytosis_MF1_Frequency = `MFI phagocytosis`*`Proportion of Phagocytosis`) %>%
  ggplot(aes(x=factor(`Carriage Status`,
                      levels=c('Carriage Negative','Carriage Positive')),
             y=log10(Product_Phagocytosis_MF1_Frequency),
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
  dplyr::filter(Time=="45min",
                #Visit=="Week 1",
                Stimulant=="LPS") %>%
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


Fig8b_new <- Neut_CD11b_data %>%
  ggplot(aes(x=factor(`Carriage Status`,levels=c('Carriage Negative','Carriage Positive')),
             y=as.numeric(`Neut CD11b+CD63-`),
             fill=factor(`Carriage Status`,
                         levels=c('Carriage Negative','Carriage Positive'))))+
  geom_violin(scale = 'width',trim = T)+
  #geom_boxplot(width=0.1,notch = T, outlier.shape = NA, fill='white')+
  #geom_jitter()+
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
       title = 'a.',
       y=expression("Proportion of CD11b"^"+"~"Neutrophils"),
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
Fig8b_new
ggsave(Fig8b_new,filename="New_Figures/Figure 8/Fig8b_new.png",
       width = 15,height = 8,dpi = 300,units = "in") 

#################################### TABLE 1 ###################################
# Demographics table for samples included in the analysis 
Fig1c_samples <- Fig1c$data %>% dplyr::pull(`LAB ID`)
Fig1d_samples <- Fig1d$data %>% dplyr::pull(`LAB ID`)
Fig1e_samples <- Fig1e$data %>% dplyr::pull(`LAB ID`)
Fig2a_samples <- Fig2a$data %>% dplyr::pull(`LAB ID`)
Micro <- read_csv("Data/Main_Files_Thesis/Micro.csv")
Fig2a_PIDs <- Micro %>% dplyr::filter(`LAB ID` %in% Fig2a_samples) %>% dplyr::pull(PID)
Fig2a_samples <- Micro %>% dplyr::filter(`PID` %in% Fig2a_PIDs,
                                         Visit=="Week 1") %>% dplyr::pull(`LAB ID`)
Fig3c_samples <- Fig3c$data %>% dplyr::filter(`HIV Status`!="NA",`Myeloperoxidase (pg/mL)`!="NA") %>% dplyr::pull(`LAB ID`)
Fig3d_samples <- Fig3d$data %>% filter(`Epithelial count`>100,
                                       `CD45+ count`>200,
                                       `Myeloperoxidase (pg/mL)`!="NA",
                                       `HIV Status`!="NA") %>% 
  dplyr::pull(`LAB ID`)

scRNA_seq_samples <- unique(all_merged_subset_labelled_new$sample)
scRNA_seq_PIDs <- Micro %>% dplyr::filter(`LAB ID` %in% scRNA_seq_samples) %>% dplyr::pull(PID)
scRNA_seq_samples <- Micro %>% dplyr::filter(`PID` %in% scRNA_seq_PIDs,
                                             Visit=="Week 1") %>% dplyr::pull(`LAB ID`)  


samples_included <- unique(c(Fig1c_samples,Fig1d_samples,Fig1e_samples,Fig2a_samples,Fig3c_samples,Fig3d_samples,scRNA_seq_samples))


Demograph_data <- Micro %>%
  dplyr::filter(Visit=="Week 1",
                `LAB ID` %in% samples_included) %>%
  merge(Lims_Data, by=c("LAB ID"),all=F) %>%
  merge(Clinical_Data, by=c("LAB ID"),all=T) %>%
  dplyr::select(PID,`LAB ID`,`Age`,`Sex`,`HIV Status`,`Carriage Status`, U5, `HIV duration`,
                `ART duration`, Location,`Carriage Status`,Density,Abs_CD4_Count,HIV_VL,`Possession index`) %>%
  dplyr::filter(`HIV Status`!="NA") %>%
  dplyr::mutate(Age = as.numeric(Age),
                Abs_CD4_Count = as.numeric(Abs_CD4_Count),
                `Carriage Status`= factor(`Carriage Status`),
                `HIV Status`=factor(`HIV Status`),
                `HIV duration`=as.numeric(`HIV duration`),
                `ART duration` = as.numeric(`ART duration`),
                Sex = factor(Sex),
                HIV_VL=as.numeric(HIV_VL),
                VL_category = ifelse(is.na(HIV_VL),"Not detected","Detected")) %>%
  dplyr::select(`LAB ID`,Age,Sex,`HIV Status`,U5,`Possession index`,`HIV duration`,`ART duration`,`Carriage Status`,
                Abs_CD4_Count,VL_category) %>%
  dplyr::filter(`LAB ID`!="CUF10Q")


group_counts <- Demograph_data %>%
  dplyr::group_by(`HIV Status`) %>%
  dplyr::summarize(n = dplyr::n(), .groups = "drop")
group_counts


Demographics_table <- Demograph_data %>%
  dplyr::select(-`LAB ID`) %>%
  gtsummary::tbl_summary(
    by = `HIV Status`, # Grouping variable (optional)
    statistic = list(all_continuous() ~ "{median} ({p25},{p75})",  # Mean (SD) for continuous variables
                     all_categorical() ~ "{n} ({p}%)"),   # Counts and percentages for categorical variables
    missing = "no" # Specify how to handle missing data
  ) %>%
  #add_ci() %>% # Add confidence intervals
  add_p() # Add p-values
Demographics_table  
Demographics_table %>%
  gtsummary::as_flex_table() %>%
  flextable::save_as_docx(path = "New_Figures/Table 1/Demographics_table.docx")


samples_with_vl <- Demograph_data %>%
  dplyr::filter(HIV_VL>0) %>%
  dplyr::group_by(`HIV Status`) %>%
  dplyr::summarise(n=dplyr::n(),.groups = "drop")
samples_with_vl

#################################### SUPPLEMENTARY TABLE 1 #####################
library(dplyr)
library(tidyr)

Single_cell_data <- read_csv("Data/Single_Cell_Data/Single_cell_data.csv")

# Prepare the data
scDemographics_data <- Single_cell_data %>%
  dplyr::select(-`LAB ID`,PID) %>%
  dplyr::mutate(Age=as.numeric(as.character(Age)),
                Sex=factor(Sex),
                `HIV Status`=factor(`HIV Status`),
                `Blood HIV viral load category`=factor(`Blood HIV viral load category`),
                `Nasal HIV Viral Load`=factor(`Nasal HIV Viral Load`),
                `Nasal Respiratory Viruses`=factor(`Nasal Respiratory Viruses`))


# Summarize Age
age_summary <- scDemographics_data %>%
  dplyr::group_by(`HIV Status`) %>%
  dplyr::summarise(
    Median_Age = median(Age, na.rm = T),
    P25_Age = quantile(Age, 0.25,na.rm=T),
    P75_Age = quantile(Age, 0.75, na.rm=T)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    Age_Summary = sprintf("%.1f (%.1f, %.1f)", Median_Age, P25_Age, P75_Age)
  ) %>%
  dplyr::select(`HIV Status`,Age_Summary)
# Summarize Age
age_summary <- scDemographics_data %>%
  group_by(`HIV Status`) %>%
  summarise(
    Median_Age = median(Age, na.rm = TRUE),
    P25_Age = quantile(Age, 0.25, na.rm = TRUE),
    P75_Age = quantile(Age, 0.75, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    Age_Summary = sprintf("%.1f (%.1f, %.1f)", Median_Age, P25_Age, P75_Age)
  ) %>%
  dplyr::select(`HIV Status`, Age_Summary)

# Summarize Sex (fixed pipeline)
sex_summary <- scDemographics_data %>%
  dplyr::group_by(`HIV Status`, Sex) %>%
  dplyr::summarise(n = n(), .groups = "drop_last") %>%
  dplyr::mutate(p = n / sum(n) * 100) %>%
  dplyr::ungroup() %>%  # Explicitly ungroup to avoid grouping issues
  dplyr::mutate(
    Sex_Summary = sprintf("%d (%.1f%%)", n, p)
  ) %>%
  dplyr::select(Sex, `HIV Status`, Sex_Summary) %>%  # Explicitly keep Sex
  tidyr::pivot_wider(names_from = `HIV Status`, values_from = Sex_Summary)

# Check the structure of sex_summary
print(colnames(sex_summary))
print(head(sex_summary))

# Compute p-values
age_pvalue <- kruskal.test(Age ~ `HIV Status`, data = scDemographics_data)$p.value
sex_pvalue <- chisq.test(table(scDemographics_data$`HIV Status`, scDemographics_data$Sex))$p.value

# Combine into a table
demographics_table <- bind_rows(
  age_summary %>% pivot_wider(names_from = `HIV Status`, values_from = Age_Summary) %>% mutate(Variable = "Age (median, IQR)", P_Value = sprintf("%.3f", age_pvalue)),
  sex_summary %>% mutate(Variable = paste("Sex:", Sex), P_Value = sprintf("%.3f", sex_pvalue))
) %>%
  dplyr::select(Variable, everything())

# Display the table
print(demographics_table)


# Demographics table
scDemographics_table <- scDemographics_data %>%
  dplyr::select(-`Age`, -PID) %>%
  gtsummary::tbl_summary(
    by = `HIV Status`, # Grouping variable (optional)
    statistic = list(all_continuous() ~ "{median} ({p25},{p75})",  # Mean (SD) for continuous variables
                     all_categorical() ~ "{n} ({p}%)"),   # Counts and percentages for categorical variables
    missing = "no" # Specify how to handle missing data
  ) %>%
  #add_ci() %>% # Add confidence intervals
  add_p() # Add p-values
scDemographics_table  
scDemographics_table <- dplyr::as_tibble(scDemographics_table)
scDemographics_table %>%
  gtsummary::as_flex_table() %>%
  flextable::save_as_docx(path = "New_Figures/Table 1/scDemographics_table.docx")

library(officer)
library(gtsummary)
library(flextable)
ft <- flextable(scDemographics_table) %>%
  fontsize(size = 10, part = "all") %>%
  align(align = "center", part = "all") %>%
  bold(part = "header") %>%
  autofit()

scDemographics_data <- read_docx() %>%
  body_add_flextable(ft)

# Save the document
print(scDemographics_data, target = "New_Figures/Table 1/scDemographics_table.docx")