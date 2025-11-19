
###################################### Figure 1b Data ##########################
Figure_1b_data <- Fig1b@data %>%
  dplyr:: filter(`HIV Status`=="HIV-",
         `Epithelial count`>100,
         `CD45+ count`>200,
         `Visit`=='Week 1',
         `CD3 Staining`=='Stained'
  ) %>%
  dplyr::mutate(
    NER = `Neutrophil count` / `Epithelial count`,
    MER = `Monocyte count` / `Epithelial count`,
    TER = `T cell count` / `Epithelial count`
  ) %>%
  dplyr::select(PID, `LAB ID`, NER, MER, TER, `HIV Status`)
Figure_1b_data
write_csv(Figure_1b_data, "HIV-PAPER/Data/SPN-PLHIV_2025_data/Figure_1b_data.csv")

###################################### Figure 1c Data ##########################
Figure_1c_data <- Fig1c@data %>%
  dplyr::select(PID, `LAB ID`, NER, MER, `HIV Status`)
Figure_1c_data
write_csv(Figure_1c_data, "HIV-PAPER/Data/SPN-PLHIV_2025_data/Figure_1c_data.csv")

###################################### Figure 1d Data ##########################
Figure_1d_data <- Fig1d@data %>%
  dplyr::select(PID, `LAB ID`, NER, MER, `HIV Status`)
Figure_1d_data
write_csv(Figure_1d_data, "HIV-PAPER/Data/SPN-PLHIV_2025_data/Figure_1d_data.csv")

###################################### Figure 1e Data ##########################
Figure_1e_data <- Fig1e@data %>%
  dplyr::select(PID, `LAB ID`,`CD4+`, `HIV Status`)
Figure_1e_data
write_csv(Figure_1e_data, "HIV-PAPER/Data/SPN-PLHIV_2025_data/Figure_1e_data.csv")

###################################### Figure 1f Data ##########################
Figure_1f_data <- Fig1f@data %>%
  dplyr::select(PID, `LAB ID`,`CD8+`, `HIV Status`) %>%
  na.omit()
Figure_1f_data
write_csv(Figure_1f_data, "HIV-PAPER/Data/SPN-PLHIV_2025_data/Figure_1f_data.csv")

###################################### Figure 1g Data ##########################
Figure_1g_data <- Fig1g@data %>%
  dplyr::select(PID, `LAB ID`,`CD3+Mait`, `HIV Status`) %>%
  na.omit()
Figure_1g_data
write_csv(Figure_1g_data, "HIV-PAPER/Data/SPN-PLHIV_2025_data/Figure_1g_data.csv")

###################################### Figure 1h Data ##########################
Figure_1h_data <- Fig1h@data %>%
  dplyr::select(PID, `LAB ID`,`CD3+TCRgd+`, `HIV Status`) %>%
  na.omit()
Figure_1h_data
write_csv(Figure_1h_data, "HIV-PAPER/Data/SPN-PLHIV_2025_data/Figure_1h_data.csv")


###################################### Figure 1i Data ##########################
Figure_1i_data <- Fig1i@data %>%
  dplyr::select(PID, `LAB ID`,`CD3+CD56+`, `HIV Status`) %>%
  na.omit()
Figure_1i_data
write_csv(Figure_1i_data, "HIV-PAPER/Data/SPN-PLHIV_2025_data/Figure_1i_data.csv")

###################################### Figure 1 Data ###########################

Figure_1_data <- Figure_1b_data %>%
  merge(Figure_1c_data, by=c("PID", "LAB ID", "NER", "MER", "HIV Status"), all=T) %>%
  merge(Figure_1d_data, by=c("PID", "LAB ID", "NER", "MER", "HIV Status"), all=T) %>%
  merge(Figure_1e_data, by=c("PID", "LAB ID", "HIV Status"), all=T) %>%
  merge(Figure_1f_data, by=c("PID", "LAB ID", "HIV Status"), all=T) %>%
  merge(Figure_1g_data, by=c("PID", "LAB ID", "HIV Status"), all=T) %>%
  merge(Figure_1h_data, by=c("PID", "LAB ID", "HIV Status"), all=T) %>%
  merge(Figure_1i_data, by=c("PID", "LAB ID", "HIV Status"), all=T) %>%
  dplyr::distinct(`LAB ID`, .keep_all = T)
write_csv(Figure_1_data, "HIV-PAPER/Data/SPN-PLHIV_2025_data/Figure_1_data.csv")

################################### Figure 2 data ##############################
Cytokines <- read_csv("Data/Main_Files_Thesis/Cytokines_averages_zeros_replaced_with_minimum_values.csv")
Micro <- read_csv("Data/Main_Files_Thesis/Micro.csv")

Neutrophil_Monocytes_data <- Masterfile %>%
  dplyr::select(`LAB ID`,NER, MER, `CD45+ count`, `Epithelial count`)

colnames(Neutrophil_Monocytes_data) <- c("LAB ID", "Neutrophils", "Monocytes", "CD45+ count", "Epithelial count")

Cytokines <- Cytokines %>%
  merge(Micro, by=c("LAB ID"),all=F) %>%
  merge(Neutrophil_Monocytes_data, by=c("LAB ID"),all=F)

colnames(Cytokines)

Cytokines <- Cytokines %>%
  dplyr::select(
    -`Epithelial count`,
    -`CD45+ count`, 
    -`Carriage Status`,
    -Serotype,
    -Density,
    -serotype_Kit_latex_reaction,
    -Visit
  ) %>%
  na.omit()

write_csv(Cytokines, "HIV-PAPER/Data/SPN-PLHIV_2025_data/Figure_2_data.csv")

Figure_2c_data <- Fig2c@data %>%
  dplyr::select(PID, `LAB ID`, `Myeloperoxidase (pg/mL)`, `HIV Status`) %>%
  na.omit()
Figure_2c_data
write_csv(Figure_2c_data, "HIV-PAPER/Data/SPN-PLHIV_2025_data/Figure_2c_data.csv")

Figure_2d_data <- Fig2d@data %>%
  dplyr::select(PID, `LAB ID`, `Myeloperoxidase (pg/mL)`, NER, `HIV Status`) %>%
  na.omit()
Figure_2d_data
write_csv(Figure_2d_data, "HIV-PAPER/Data/SPN-PLHIV_2025_data/Figure_2d_data.csv")  


Figure_2e_data <- Fig2e@data %>%
  dplyr::select(PID, `LAB ID`, `Myeloperoxidase (pg/mL)`, MER, `HIV Status`) %>%
  na.omit()
Figure_2e_data
write_csv(Figure_2e_data, "HIV-PAPER/Data/SPN-PLHIV_2025_data/Figure_2e_data.csv")  

# Supplementary Figure 2 ####
Supp_Figure_2_data <- Figure_2_data %>%
  dplyr::select(
    `LAB ID`,
    `TNF-a`, 
    `IL-6`, 
    `IL-8`, 
    `IL-1B`, 
    `IFN-y`, 
    `IL-2`, 
    `IL-12p70`,
    `IL-13`, 
    `IL-4`, 
    `IL-10`, 
    `HIV Status`) %>%
  dplyr::mutate(
    Visit = ifelse(grepl("CUH", `LAB ID`),
                    "Week 3", "Week 4")
  ) %>%
  dplyr::filter(
    Visit == "Week 3"
  )

write_csv(Supp_Figure_2_data, file="HIV-PAPER/Data/SPN-PLHIV_2025_data/Supp_Figure_2_data.csv")
################################### Figure 3 data ##############################
all_merged_subset_labelled_new <- readRDS("Data/Single_Cell_Data/all_merged_subset_labelled_new.rds")

Figure_3_data <- all_merged_subset_labelled_new

saveRDS(Figure_3_data, "HIV-PAPER/Data/SPN-PLHIV_2025_data/Figure_3_data.rds")
################################### Figure 4 data ##############################
saveRDS(Figure_4_data, "HIV-PAPER/Data/SPN-PLHIV_2025_data/Figure_4_data.rds")
 # Figure 7 #############
Figure_7a_data <- Fig7a@data %>%
  dplyr::select(PID, `LAB ID`,`HIV Status`, NER, `Carriage Status`)
Figure_7a_data
write_csv(Figure_7a_data, "HIV-PAPER/Data/SPN-PLHIV_2025_data/Figure_7a_data.csv")

Figure_7b_data <- Fig7b@data %>%
  dplyr::select(PID, `LAB ID`,`HIV Status`, MER, `Carriage Status`)
Figure_7b_data
write_csv(Figure_7b_data, "HIV-PAPER/Data/SPN-PLHIV_2025_data/Figure_7b_data.csv")

Figure_7c_data <- Fig7c@data %>%
  dplyr::select(PID, `LAB ID`,`HIV Status`, `Neutrophil CD11b++`, Density, `Carriage Status`)
Figure_7c_data
write_csv(Figure_7c_data, "HIV-PAPER/Data/SPN-PLHIV_2025_data/Figure_7c_data.csv")

Figure_7d_data <- Fig7d@data %>%
  dplyr::select(PID, `LAB ID`,`HIV Status`, MER, Density, `Carriage Status`)
Figure_7d_data
write_csv(Figure_7d_data, "HIV-PAPER/Data/SPN-PLHIV_2025_data/Figure_7d_data.csv")
