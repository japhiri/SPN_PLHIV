# Table 1 ----
library(readr)
library(tidyverse)
library(gtsummary)
library(flextable)
# Demographics table for samples included in the analysis 
Micro <- read_csv("Data/Main_Files_Thesis/Micro.csv")
Samples_Included_SPN_PLHIV_Paper <- read_csv("Data/Main_Files_Thesis/Samples_Included_SPN_PLHIV_Paper.csv")
Clinical_Data <- read_csv("Data/Main_Files_Thesis/Clinical_Data.csv")
Lab_Data <- read_csv("Data/Main_Files_Thesis/Lab_Data.csv")
Demograph_data <- Samples_Included_SPN_PLHIV_Paper %>%
  merge(Lab_Data, by=c("LAB ID"),all=T) %>%
  merge(Clinical_Data, by=c("LAB ID"),all=T) %>%
  dplyr::filter(`LAB ID` %in% Samples_Included_SPN_PLHIV_Paper$`LAB ID`) %>%
  dplyr::rename(Watch=watch,
                Radio=radio,
                Bank_account=bank,
                Electric_iron=iron,
                Sewing_machine=sew,
                Mobile_phone=mobile,
                CD_player=cd,
                Electric_fan=fanelec,
                Mosquito_net=netyn,
                Mattress=mattress,
                Bed=bed,
                Bicycle=bike,
                Motorcycle=motor,
                Car=car,
                TV=tv) %>%
  dplyr::mutate(across(c(Watch, Electric_iron, Bank_account, Sewing_machine,
                         Mobile_phone, CD_player, Electric_fan, Mosquito_net, Mattress,
                         Bed, Bicycle, Motorcycle, Car, TV),
                       ~ as.numeric(.))) %>%
  dplyr::mutate(Possession_index = rowSums(across(c(Watch, Electric_iron, Bank_account, Sewing_machine,
                                        Mobile_phone, CD_player, Electric_fan, Mosquito_net, Mattress,
                                        Bed, Bicycle, Motorcycle, Car, TV)), na.rm = TRUE)) %>% 
  dplyr::mutate(
    HIV_VL = case_when(
      `HIV Status` == "HIV-" ~ NA_character_,                                       # HIV-negative → NA
      `HIV Status` %in% c("PLHIV ART <3m", "PLHIV ART >1y") & !is.na(as.numeric(HIV_VL)) ~ "Detected", # numeric → Detected
      `HIV Status` %in% c("PLHIV ART <3m", "PLHIV ART >1y") & is.na(as.numeric(HIV_VL))  ~ "Not Detected",
      TRUE ~ HIV_VL))

Demographics_table <- Demograph_data %>%
  dplyr::select(`HIV Status`, Age, Sex, Under_5, Possession_index,  `ART duration`, Abs_CD4_Count, HIV_VL) %>%
  dplyr::mutate(
    Age = as.numeric(Age),
    Abs_CD4_Count = as.numeric(Abs_CD4_Count),
    `ART duration` = as.numeric(`ART duration`),
    `HIV Status`= factor(`HIV Status`, levels = c("HIV-", "PLHIV ART <3m", "PLHIV ART >1y")),
    Possession_index = as.numeric(Possession_index),
    Sex = factor(Sex),
    HIV_VL = factor(HIV_VL,levels = c("Not Detected","Detected"))
  ) %>%
  gtsummary::tbl_summary(
    by = `HIV Status`, # Grouping variable (optional)
    statistic = list(all_continuous() ~ "{median} ({p25},{p75})",  # Mean (SD) for continuous variables
                     all_categorical() ~ "{n} ({p}%)"),   # Counts and percentages for categorical variables
    missing = "no" # Specify how to handle missing data
  ) %>%
  #add_ci() %>% # Add confidence intervals
  add_p(pvalue_fun = identity) # Add p-values
Demographics_table 
Demographics_table %>%
  gtsummary::as_flex_table() %>%
  flextable::save_as_docx(path = "New_Figures/Table 1/Demographics_table.docx")
