########### Install and Load required packages ##########

install.packages("emmeans")
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggpubr)
########## MIXED EFFECTS MODELS FOR NER ACROSS VISISTS #########
########## Ordering the Visit column is factor ordered 
data <- Masterfile %>%
  dplyr::mutate(
    Visit = factor(Visit,
                   levels = c("Week 1",
                              "Week 2",
                              "Week 3",
                              "Week 4",
                              "Week 5"),
                   ordered = TRUE)
  ) %>%
  dplyr::mutate(`HIV Status` = as.factor(`HIV Status`)
               ) %>%
  dplyr::select(
    `LAB ID`,
    `PID`,
    `Visit`,
    `Neutrophil proportion`,
    `Monocyte proportion`,
    `MER`,
    `NER`,
    `CD4+`,
    `CD8+`,
    `HIV Status`,
    `Carriage Status`,
    `Density`
  ) %>%
  dplyr::filter(
    NER!="NA",
    PID!="NA",
    `Neutrophil proportion`!="NA"
  )
########### Fit mixed-effects model (NER ~ Visit + random intercept for PID) 
fit <- lmer(NER ~ Visit + (1 | PID), data = data)
summary(fit)
########## Estimated marginal means (population-adjusted means per Visit) 
emm <- emmeans(fit, ~ Visit)
emm_df <- as.data.frame(emm)
########## Plot adjusted means with 95% CI 
p_emm <- ggplot(emm_df, aes(x = Visit, y = emmean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  labs(
    title = "Estimated Marginal Means of NER by Visit",
    x = "Visit",
    y = "Estimated Marginal Mean of NER"
  ) +
  theme_minimal()
print(p_emm)
########## Remove subject-level offsets for boxplot 
ranefs <- ranef(fit)$PID %>%
  rownames_to_column("PID") %>%
  rename(b_intercept = `(Intercept)`)


data_adj <- data %>%
  left_join(ranefs, by = "PID") %>%
  dplyr::mutate(
    Corrected_NER = NER - b_intercept
    )
##########Boxplot of corrected values
p_box <- ggplot(data_adj, aes(x = Visit, y = Corrected_NER))+
  geom_boxplot(outlier.shape = NA, width = 0.6)+
  geom_jitter(width = 0.15, alpha = 0.35, size = 1)+
  geom_pwc(hide.ns = T)+
  labs(title = "Boxplots after removing subject-level offsets",
       y = "NER (subject corrected)")+
  scale_y_log10()+
  theme_classic()
print(p_box)
########## Optional: spaghetti plot of raw values by PID 
p_spag <- ggplot(data, aes(Visit, NER, group=PID))+
  geom_line(alpha = 0.25) +
  stat_summary(fun = mean, geom = "point", size = 0.25, alpha=0.35) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.15, alpha=0.35) +
  labs(title = "Within-subject NER trajectories",
       y = "Neutrophil to Epithelial cell Ratio") +
  scale_y_log10()+
  theme_classic()
print(p_spag)
########## Post-hoc tests between visits 
pairs(emm, adjust = "holm")
######### Sow plots 
print(p_emm)
print(p_box)
print(p_spag)

######### Fitting model on NER with HIV as a grouping factor #########
fit_NER <- lmer(NER ~ Visit * `HIV Status` * `Carriage Status` + (1 | PID), data = data)
summary(fit_NER)
########## Estimated marginal means (adjusted means) 
emm_NER <- emmeans(fit_NER, ~ Visit | `HIV Status`|`Carriage Status`)
emm_df <- as.data.frame(emm_NER)
########## Pairwise comparisons of HIV groups at each. visit
contrast(emm_NER, method = "pairwise", by = "Visit", adjust = "holm")
########## Plot adjusted means
ggplot(emm_df, aes(Visit, emmean, color = `HIV Status`, group = `HIV Status`))+
  geom_point(position = position_dodge(width = 0.3), size = 3)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position = position_dodge(width = 0.3), width = 0.15)+
  labs(
    title = "Adjusted NER by Visit and HIV Status",
    y= "Adjusted mean NER (95% CI)"
  )+
  facet_wrap(~`Carriage Status`)+
  theme_classic()
########## subject-corrected boxplot by HIV status
ranefs <- ranef(fit_NER)$PID %>%
  rownames_to_column("PID") %>%
  rename(b_intercept = `(Intercept)`)

data_adj <- data %>%
  left_join(ranefs, by = "PID") %>%
  dplyr::mutate(NER_corrected = NER - b_intercept)

NER_plot <- ggplot(data_adj, aes(Visit, log10(NER_corrected), color = `HIV Status`)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, position = position_dodge(width = 0.7)) +
  geom_jitter(alpha = 0.3, 
              size = 1, 
              position = position_jitterdodge(jitter.width=0.15, dodge.width=0.7)
              ) +
  labs(title = "NER (subject-corrected) by Visit and HIV Status",
       y = "Neutrophil to Epithelial Cell Ratio") +
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  geom_pwc(
    hide.ns = T
  )+
  #facet_wrap(~`Carriage Status`)+
  #scale_y_log10() +
  theme_classic()
print(NER_plot)
######### Fitting model on MER with HIV as a grouping factor #########
fit_MER <- lmer(MER ~ Visit * `HIV Status` + (1 | PID), data = data)
summary(fit_MER)
########## Estimated marginal means (adjusted means) 
emm_MER <- emmeans(fit_MER, ~ Visit | `HIV Status`)
emm_df <- as.data.frame(emm_MER)
########## Pairwise comparisons of HIV groups at each. visit
contrast(emm_MER, method = "pairwise", by = "Visit", adjust = "holm")
########## Plot adjusted means
ggplot(emm_df, aes(Visit, emmean, color = `HIV Status`, group = `HIV Status`))+
  geom_point(position = position_dodge(width = 0.3), size = 3)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position = position_dodge(width = 0.3), width = 0.15)+
  labs(
    title = "Adjusted MER by Visit and HIV Status",
    y= "Adjusted mean MER (95% CI)"
  )+
  theme_classic()
########## subject-corrected boxplot by HIV status
ranefs <- ranef(fit_MER)$PID %>%
  rownames_to_column("PID") %>%
  rename(b_intercept = `(Intercept)`)

data_adj <- data %>%
  left_join(ranefs, by = "PID") %>%
  dplyr::mutate(MER_corrected = MER - b_intercept)

MER_plot <- ggplot(data_adj, aes(Visit, log10(MER_corrected), color = `HIV Status`)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, position = position_dodge(width = 0.7)) +
  geom_jitter(alpha = 0.3, 
              size = 1, 
              position = position_jitterdodge(jitter.width=0.15, dodge.width=0.7)
  ) +
  labs(title = "MER (subject-corrected) by Visit and HIV Status",
       y = "Monocyte Epithelial Ratio") +
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  geom_pwc(
    hide.ns = T
  )+
  #scale_y_log10() +
  theme_bw()+
  theme(legend.position = "none")
print(MER_plot)
######### Fitting model on Neutrophil proportion with HIV as a grouping factor #########
fit_Neut <- lmer(`Neutrophil proportion` ~ Visit * `HIV Status` + (1 | PID), data = data)
summary(fit_Neut)
########## Estimated marginal means (adjusted means) 
emm_Neut <- emmeans(fit_Neut, ~ Visit | `HIV Status`)
emm_df <- as.data.frame(emm_Neut)
########## Pairwise comparisons of HIV groups at each. visit
contrast(emm_Neut, method = "pairwise", by = "Visit", adjust = "holm")
########## Plot adjusted means
ggplot(emm_df, aes(Visit, emmean, color = `HIV Status`, group = `HIV Status`))+
  geom_point(position = position_dodge(width = 0.3), size = 3)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position = position_dodge(width = 0.3), width = 0.15)+
  labs(
    title = "Adjusted Neutrophil Proportion by Visit and HIV Status",
    y= "Adjusted mean Neutrophil Proportion (95% CI)"
  )+
  theme_classic()
########## subject-corrected boxplot by HIV status
ranefs <- ranef(fit_Neut)$PID %>%
  rownames_to_column("PID") %>%
  rename(b_intercept = `(Intercept)`)

data_adj <- data %>%
  left_join(ranefs, by = "PID") %>%
  dplyr::mutate(Neut_corrected = `Neutrophil proportion` - b_intercept)

Neut_plot <- ggplot(data_adj, aes(Visit, Neut_corrected, color = `HIV Status`)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, position = position_dodge(width = 0.7)) +
  geom_jitter(alpha = 0.3, 
              size = 1, 
              position = position_jitterdodge(jitter.width=0.15, dodge.width=0.7)
  ) +
  labs(title = "Neutrophil proportion (subject-corrected) by Visit and HIV Status",
       y = "Neutrophil Proportion") +
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  geom_pwc(
    hide.ns = T
  )+
  #scale_y_log10() +
  theme_bw()+
  theme(legend.position = "none")
print(Neut_plot)
######### Fitting model on Monocyte proportion with HIV as a grouping factor #########
fit_Mono <- lmer(`Monocyte proportion` ~ Visit * `HIV Status` + (1 | PID), data = data)
summary(fit_Mono)
########## Estimated marginal means (adjusted means) 
emm_mono <- emmeans(fit_Mono, ~ Visit | `HIV Status`)
emm_df <- as.data.frame(emm_mono)
########## Pairwise comparisons of HIV groups at each. visit
contrast(emm_mono, method = "pairwise", by = "Visit", adjust = "holm")
########## Plot adjusted means
ggplot(emm_df, aes(Visit, emmean, color = `HIV Status`, group = `HIV Status`))+
  geom_point(position = position_dodge(width = 0.3), size = 3)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position = position_dodge(width = 0.3), width = 0.15)+
  labs(
    title = "Adjusted Monocyte Proportion by Visit and HIV Status",
    y= "Adjusted mean Monocyte Proportion (95% CI)"
  )+
  theme_classic()
########## subject-corrected boxplot by HIV status
ranefs <- ranef(fit_Mono)$PID %>%
  rownames_to_column("PID") %>%
  rename(b_intercept = `(Intercept)`)

data_adj <- data %>%
  left_join(ranefs, by = "PID") %>%
  dplyr::mutate(Mono_corrected = `Monocyte proportion` - b_intercept)

Mono_plot <- ggplot(data_adj, aes(Visit, Mono_corrected, color = `HIV Status`)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, position = position_dodge(width = 0.7)) +
  geom_jitter(alpha = 0.3, 
              size = 1, 
              position = position_jitterdodge(jitter.width=0.15, dodge.width=0.7)
  ) +
  labs(title = "Monocyte proportion (subject-corrected) by Visit and HIV Status",
       y = "Monocyte Proportion") +
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  geom_pwc(
    hide.ns = T
  )+
  #scale_y_log10() +
  theme_bw()+
  theme(legend.position = "none")
print(Mono_plot)
######### Fitting model on CD4+ T cell proportion with HIV as a grouping factor #########
setwd("/Users/471288/OneDrive - Malawi-Liverpool Wellcome Trust/Projects/PhD Work/Nasomune/Thesis Paper Publications/Nasopharynx vs inferior turbinate")
########### Read data
Micro <- read_csv("Data/Micro.csv")
########### Load Week 1
Innate_like_Week_1_2025 <- read_csv("Data/Innate_like_Week_1_2025.csv") %>%
  dplyr::select(
    `LAB ID`,`Sample Type`,`Sample Quality`,
    `CD4+`,`CD8+`,`DN-`,
    `CD4+CXCR3-CCR6+`,`CD4+CXCR3+CCR6+`,`CD4+CXCR3+CCR6-`,
    `CD8+CXCR3-CCR6+`,`CD8+CXCR3+CCR6+`,`CD8+CXCR3+CCR6-`,
    `CD3+ MAIT cells`,`CD3+ gd T cells`,`CD3+ CD56+ T cells`
  ) %>%
  dplyr::filter(`Sample Quality` == "Good",
                `Sample Type` %in% c("Swab","Scrape"))
########### Load Week 2
Innate_like_Week_2_2025 <- read_csv("Data/Innate_like_Week_2_2025.csv") %>%
  dplyr::select(
    `LAB ID`,`Sample Type`,`Sample Quality`,
    `CD4+`,`CD8+`,`DN-`,
    `CD4+CXCR3-CCR6+`,`CD4+CXCR3+CCR6+`,`CD4+CXCR3+CCR6-`,
    `CD8+CXCR3-CCR6+`,`CD8+CXCR3+CCR6+`,`CD8+CXCR3+CCR6-`,
    `CD3+ MAIT cells`,`CD3+ gd T cells`,`CD3+ CD56+ T cells`
  )
########### Load Week 3
Innate_like_Week_3_2025 <- read_csv("Data/Innate_like_Week_3_2025.csv") %>%
  dplyr::select(
    `LAB ID`,`Sample Type`,`Sample Quality`,
    `CD4+`,`CD8+`,`DN-`,
    `CD4+CXCR3-CCR6+`,`CD4+CXCR3+CCR6+`,`CD4+CXCR3+CCR6-`,
    `CD8+CXCR3-CCR6+`,`CD8+CXCR3+CCR6+`,`CD8+CXCR3+CCR6-`,
    `CD3+ MAIT cells`,`CD3+ gd T cells`,`CD3+ CD56+ T cells`
  )
########### Load Week 4
Innate_like_Week_4_2025 <- read_csv("Data/Innate_like_Week_4_2025.csv") %>%
  dplyr::select(
    `LAB ID`,`Sample Type`,`Sample Quality`,
    `CD4+`,`CD8+`,`DN-`,
    `CD4+CXCR3-CCR6+`,`CD4+CXCR3+CCR6+`,`CD4+CXCR3+CCR6-`,
    `CD8+CXCR3-CCR6+`,`CD8+CXCR3+CCR6+`,`CD8+CXCR3+CCR6-`,
    `CD3+ MAIT cells`,`CD3+ gd T cells`,`CD3+ CD56+ T cells`
  )
########### Load Week 5
Innate_like_Week_5_2025 <- read_csv("Data/Innate_like_Week_5_2025.csv") %>%
  dplyr::select(
    `LAB ID`,`Sample Type`,`Sample Quality`,
    `CD4+`,`CD8+`,`DN-`,
    `CD4+CXCR3-CCR6+`,`CD4+CXCR3+CCR6+`,`CD4+CXCR3+CCR6-`,
    `CD8+CXCR3-CCR6+`,`CD8+CXCR3+CCR6+`,`CD8+CXCR3+CCR6-`,
    `CD3+ MAIT cells`,`CD3+ gd T cells`,`CD3+ CD56+ T cells`
  )
########### Combine all weeks data and Merge with micro data
Innate_like_Tcells <- rbind(
  Innate_like_Week_1_2025,
  Innate_like_Week_2_2025,
  Innate_like_Week_3_2025,
  Innate_like_Week_4_2025,
  Innate_like_Week_5_2025
) %>%
  merge(Micro, by = "LAB ID") %>%
  dplyr::mutate(
    Visit = ifelse(grepl("CUF",`LAB ID`),"Week 1",
                   ifelse(grepl("CUG",`LAB ID`),"Week 2",
                          ifelse(grepl("CUH",`LAB ID`),"Week 3",
                                 ifelse(grepl("CUI",`LAB ID`),"Week 4","Week 5"))))
  ) %>%
  dplyr::filter(
    `Sample Quality`=="Good"
  ) %>%
  dplyr::mutate(
    `CD4+CXCR3+` = `CD4+CXCR3+CCR6-` + `CD4+CXCR3+CCR6+`,
    `CD4+CCR6+` = `CD4+CXCR3-CCR6+` + `CD4+CXCR3+CCR6+`,
    `CD8+CXCR3+` = `CD8+CXCR3+CCR6-` + `CD8+CXCR3+CCR6+`,
    `CD8+CCR6+` = `CD8+CXCR3-CCR6+` + `CD8+CXCR3+CCR6+`
  )


data <- Innate_like_Tcells %>%
  # Only turn "NA" strings into real NA **if** the column is character
  mutate(
    across(c(`CD4+`, `CD8+`), ~ if (is.character(.x)) na_if(.x, "NA") else .x),
    `HIV Status` = na_if(`HIV Status`, "NA"),
    Visit        = na_if(Visit, "NA")
  ) %>%
  mutate(
    Visit = factor(str_trim(Visit),
                   levels  = c("Week 1","Week 2","Week 3","Week 4","Week 5"),
                   ordered = TRUE),
    `HIV Status` = factor(str_trim(`HIV Status`)),
    # Coerce to numeric (safe even if already numeric)
    `CD4+` = suppressWarnings(as.numeric(`CD4+`)),
    `CD8+` = suppressWarnings(as.numeric(`CD8+`))
  ) %>%
  filter(
    !is.na(`CD4+`),
    !is.na(`CD8+`),
    `Sample Type` == "Scrape",
    !is.na(`HIV Status`),
    !is.na(Visit)
  ) %>%
  droplevels()


fit_CD4 <- lmer(`CD4+` ~ Visit * `HIV Status` + (1 | PID), data = data)
summary(fit_CD4)
########## Estimated marginal means (adjusted means) 
emm_CD4 <- emmeans(fit_CD4, ~ Visit | `HIV Status`)
emm_df <- as.data.frame(emm_CD4)
########## Pairwise comparisons of HIV groups at each. visit
contrast(emm_CD4, method = "pairwise", by = "Visit", adjust = "holm")
########## Plot adjusted means
ggplot(emm_df, aes(Visit, emmean, color = `HIV Status`, group = `HIV Status`))+
  geom_point(position = position_dodge(width = 0.3), size = 3)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position = position_dodge(width = 0.3), width = 0.15)+
  labs(
    title = "Adjusted CD4+ Proportion by Visit and HIV Status",
    y= "Adjusted mean CD4+ Proportion (95% CI)"
  )+
  theme_classic()
########## subject-corrected boxplot by HIV status
ranefs <- ranef(fit_CD4)$PID %>%
  rownames_to_column("PID") %>%
  rename(b_intercept = `(Intercept)`)

data_adj <- data %>%
  left_join(ranefs, by = "PID") %>%
  dplyr::mutate(CD4_corrected = `CD4+` - b_intercept)

CD4_plot <- ggplot(data_adj, aes(Visit, CD4_corrected, color = `HIV Status`)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, position = position_dodge(width = 0.7)) +
  geom_jitter(alpha = 0.4, 
              size = 7, 
              position = position_jitterdodge(jitter.width=0.15, dodge.width=0.7)
  ) +
  labs(#title = "CD4 proportion (subject-corrected) by Visit and HIV Status",
       y = "Proportion of CD4+ T cells (%)") +
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  geom_pwc(method = 'wilcox.test',
           label = "{ifelse(p < 0.0001, 'p < 0.0001', sprintf('p = %.4f', p))}",
           label.size = 8,
           tip.length = 0.01,
           hide.ns = F)+
  #scale_y_continuous(limits = c(0,150), breaks=c(0,20,40,60,80,100,120))+
  theme_bw()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 50, face = "bold"),
        axis.ticks.length = unit(0.5,"cm"),
        legend.text = element_text(size=25),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 30))
print(CD4_plot)
ggsave(CD4_plot,filename="Mixed_Models/Figure 1/CD4_plot.png",
       width = 17,height = 9,dpi = 300)
######### Fitting model on CD8+ T cell proportion with HIV as a grouping factor #########
fit_CD8 <- lmer(`CD8+` ~ Visit * `HIV Status` + (1 | PID), data = data)
summary(fit_CD8)
########## Estimated marginal means (adjusted means) 
emm_CD8 <- emmeans(fit_CD8, ~ Visit | `HIV Status`)
emm_df <- as.data.frame(fit_CD8)
########## Pairwise comparisons of HIV groups at each. visit
contrast(emm_CD8, method = "pairwise", by = "Visit", adjust = "holm")
########## Plot adjusted means
ggplot(emm_df, aes(Visit, emmean, color = `HIV Status`, group = `HIV Status`))+
  geom_point(position = position_dodge(width = 0.3), size = 3)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position = position_dodge(width = 0.3), width = 0.15)+
  labs(
    title = "Adjusted CD8 Proportion by Visit and HIV Status",
    y= "Adjusted mean CD8 Proportion (95% CI)"
  )+
  theme_classic()
########## subject-corrected boxplot by HIV status
ranefs <- ranef(fit_CD8)$PID %>%
  rownames_to_column("PID") %>%
  rename(b_intercept = `(Intercept)`)

data_adj <- data %>%
  left_join(ranefs, by = "PID") %>%
  dplyr::mutate(CD8_corrected = `CD8+` - b_intercept)

CD8_plot <- ggplot(data_adj, aes(Visit, CD8_corrected, color = `HIV Status`)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, position = position_dodge(width = 0.7)) +
  geom_jitter(alpha = 0.4, 
              size = 7, 
              position = position_jitterdodge(jitter.width=0.15, dodge.width=0.7)
  ) +
  labs(#title = "CD4 proportion (subject-corrected) by Visit and HIV Status",
    y = "Proportion of CD8+ T cells (%)") +
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  geom_pwc(method = 'wilcox.test',
           label = "{ifelse(p < 0.0001, 'p < 0.0001', sprintf('p = %.4f', p))}",
           label.size = 8,
           tip.length = 0.01,
           hide.ns = F)+
  scale_y_continuous(limits = c(0,130), breaks=c(0,20,40,60,80,100,120))+
  theme_bw()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 50, face = "bold"),
        axis.ticks.length = unit(0.5,"cm"),
        legend.text = element_text(size=25),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 30))
print(CD8_plot)
ggsave(CD8_plot,filename="Mixed_Models/Figure 1/CD8_plot.png",
       width = 17,height = 9,dpi = 300)

######### Fitting model on CD4+CXCR3+ T cell proportion with HIV as a grouping factor #########
fit_CD4CXCR3 <- lmer(`CD4+CXCR3+` ~ Visit * `HIV Status` + (1 | PID), data = data)
summary(fit_CD4CXCR3)
########## Estimated marginal means (adjusted means) 
emm_CD4CXCR3 <- emmeans(fit_CD4CXCR3, ~ Visit | `HIV Status`)
emm_df <- as.data.frame(emm_CD4CXCR3)
########## Pairwise comparisons of HIV groups at each. visit
contrast(emm_CD4CXCR3, method = "pairwise", by = "Visit", adjust = "holm")
########## Plot adjusted means
ggplot(emm_df, aes(Visit, emmean, color = `HIV Status`, group = `HIV Status`))+
  geom_point(position = position_dodge(width = 0.3), size = 3)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position = position_dodge(width = 0.3), width = 0.15)+
  labs(
    title = "Adjusted CD4+CXCR3+ Proportion by Visit and HIV Status",
    y= "Adjusted mean CD4+CXCR3+ Proportion (95% CI)"
  )+
  theme_classic()
########## subject-corrected boxplot by HIV status
ranefs <- ranef(fit_CD4CXCR3)$PID %>%
  rownames_to_column("PID") %>%
  rename(b_intercept = `(Intercept)`)

data_adj <- data %>%
  left_join(ranefs, by = "PID") %>%
  dplyr::mutate(CD4CXCR3_corrected = `CD4+CXCR3+` - b_intercept)

CD4CXCR3_plot <- ggplot(data_adj, aes(Visit, CD4CXCR3_corrected, color = `HIV Status`)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, position = position_dodge(width = 0.7)) +
  geom_jitter(alpha = 0.4, 
              size = 7, 
              position = position_jitterdodge(jitter.width=0.15, dodge.width=0.7)
  ) +
  labs(#title = "CD4 proportion (subject-corrected) by Visit and HIV Status",
    y = "Proportion of CD4+CXCR3+ T cells (%)") +
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  geom_pwc(method = 'wilcox.test',
           label = "{ifelse(p < 0.0001, 'p < 0.0001', sprintf('p = %.4f', p))}",
           label.size = 8,
           tip.length = 0.01,
           hide.ns = F)+
  #scale_y_continuous(limits = c(0,150), breaks=c(0,20,40,60,80,100,120))+
  theme_bw()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 50, face = "bold"),
        axis.ticks.length = unit(0.5,"cm"),
        legend.text = element_text(size=25),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 30))
print(CD4CXCR3_plot)
ggsave(CD4CXCR3_plot,filename="Mixed_Models/Figure 1/CD4CXCR3_plot.png",
       width = 17,height = 9,dpi = 300)
######### Fitting model on CD8+CCR6+ T cell proportion with HIV as a grouping factor #########
fit_CD8CCR6 <- lmer(`CD8+CCR6+` ~ Visit * `HIV Status` + (1 | PID), data = data)
summary(fit_CD8CCR6)
########## Estimated marginal means (adjusted means) 
emm_CD8CCR6 <- emmeans(fit_CD8CCR6, ~ Visit | `HIV Status`)
emm_df <- as.data.frame(emm_CD8CCR6)
########## Pairwise comparisons of HIV groups at each. visit
contrast(emm_CD8CCR6, method = "pairwise", by = "Visit", adjust = "holm")
########## Plot adjusted means
ggplot(emm_df, aes(Visit, emmean, color = `HIV Status`, group = `HIV Status`))+
  geom_point(position = position_dodge(width = 0.3), size = 3)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position = position_dodge(width = 0.3), width = 0.15)+
  labs(
    title = "Adjusted CD8+CCR6+ Proportion by Visit and HIV Status",
    y= "Adjusted mean CD8+CCR6+ Proportion (95% CI)"
  )+
  theme_classic()
########## subject-corrected boxplot by HIV status
ranefs <- ranef(fit_CD8CCR6)$PID %>%
  rownames_to_column("PID") %>%
  rename(b_intercept = `(Intercept)`)

data_adj <- data %>%
  left_join(ranefs, by = "PID") %>%
  dplyr::mutate(CD8CCR6_corrected = `CD8+CCR6+` - b_intercept)

CD8CCR6_plot <- ggplot(data_adj, aes(Visit, CD8CCR6_corrected, color = `HIV Status`)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, position = position_dodge(width = 0.7)) +
  geom_jitter(alpha = 0.4, 
              size = 7, 
              position = position_jitterdodge(jitter.width=0.15, dodge.width=0.7)
  ) +
  labs(#title = "CD4 proportion (subject-corrected) by Visit and HIV Status",
    y = "Proportion of CD8+CCR6+ T cells (%)") +
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  geom_pwc(method = 'wilcox.test',
           label = "{ifelse(p < 0.0001, 'p < 0.0001', sprintf('p = %.4f', p))}",
           label.size = 8,
           tip.length = 0.01,
           hide.ns = F)+
  #scale_y_continuous(limits = c(0,150), breaks=c(0,20,40,60,80,100,120))+
  theme_bw()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 50, face = "bold"),
        axis.ticks.length = unit(0.5,"cm"),
        legend.text = element_text(size=25),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 30))
print(CD8CCR6_plot)
ggsave(CD8CCR6_plot,filename="Mixed_Models/Figure 1/CD8CCR6_plot.png",
       width = 17,height = 9,dpi = 300)
######### Fitting model on CD8+CXCR3+ T cell proportion with HIV as a grouping factor #########
fit_CD8CXCR3 <- lmer(`CD8+CXCR3+` ~ Visit * `HIV Status` + (1 | PID), data = data)
summary(fit_CD8CXCR3)
########## Estimated marginal means (adjusted means) 
emm_CD8CXCR3 <- emmeans(fit_CD8CXCR3, ~ Visit | `HIV Status`)
emm_df <- as.data.frame(emm_CD8CXCR3)
########## Pairwise comparisons of HIV groups at each. visit
contrast(emm_CD8CXCR3, method = "pairwise", by = "Visit", adjust = "holm")
########## Plot adjusted means
ggplot(emm_df, aes(Visit, emmean, color = `HIV Status`, group = `HIV Status`))+
  geom_point(position = position_dodge(width = 0.3), size = 3)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position = position_dodge(width = 0.3), width = 0.15)+
  labs(
    title = "Adjusted CD8+CXCR3+ Proportion by Visit and HIV Status",
    y= "Adjusted mean CD8+CXCR3+ Proportion (95% CI)"
  )+
  theme_classic()
########## subject-corrected boxplot by HIV status
ranefs <- ranef(fit_CD8CXCR3)$PID %>%
  rownames_to_column("PID") %>%
  rename(b_intercept = `(Intercept)`)

data_adj <- data %>%
  left_join(ranefs, by = "PID") %>%
  dplyr::mutate(CD8CXCR3_corrected = `CD8+CXCR3+` - b_intercept)

CD8CXCR3_plot <- ggplot(data_adj, aes(Visit, CD8CXCR3_corrected, color = `HIV Status`)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, position = position_dodge(width = 0.7)) +
  geom_jitter(alpha = 0.4, 
              size = 7, 
              position = position_jitterdodge(jitter.width=0.15, dodge.width=0.7)
  ) +
  labs(#title = "CD4 proportion (subject-corrected) by Visit and HIV Status",
    y = "Proportion of CD8+CXCR3+ T cells (%)") +
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  geom_pwc(method = 'wilcox.test',
           label = "{ifelse(p < 0.0001, 'p < 0.0001', sprintf('p = %.4f', p))}",
           label.size = 8,
           tip.length = 0.01,
           hide.ns = F)+
  #scale_y_continuous(limits = c(0,150), breaks=c(0,20,40,60,80,100,120))+
  theme_bw()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 50, face = "bold"),
        axis.ticks.length = unit(0.5,"cm"),
        legend.text = element_text(size=25),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 30))
print(CD8CXCR3_plot)
ggsave(CD8CXCR3_plot,filename="Mixed_Models/Figure 1/CD8CXCR3_plot.png",
       width = 17,height = 9,dpi = 300)
######### Fitting model on CD8+CCR6+ T cell proportion with HIV as a grouping factor #########
fit_CD8CCR6 <- lmer(`CD8+CCR6+` ~ Visit * `HIV Status` + (1 | PID), data = data)
summary(fit_CD8CCR6)
########## Estimated marginal means (adjusted means) 
emm_CD8CCR6 <- emmeans(fit_CD8CCR6, ~ Visit | `HIV Status`)
emm_df <- as.data.frame(emm_CD8CCR6)
########## Pairwise comparisons of HIV groups at each. visit
contrast(emm_CD8CCR6, method = "pairwise", by = "Visit", adjust = "holm")
########## Plot adjusted means
ggplot(emm_df, aes(Visit, emmean, color = `HIV Status`, group = `HIV Status`))+
  geom_point(position = position_dodge(width = 0.3), size = 3)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position = position_dodge(width = 0.3), width = 0.15)+
  labs(
    title = "Adjusted CD8+CCR6+ Proportion by Visit and HIV Status",
    y= "Adjusted mean CD8+CCR6+ Proportion (95% CI)"
  )+
  theme_classic()
########## subject-corrected boxplot by HIV status
ranefs <- ranef(fit_CD8CCR6)$PID %>%
  rownames_to_column("PID") %>%
  rename(b_intercept = `(Intercept)`)

data_adj <- data %>%
  left_join(ranefs, by = "PID") %>%
  dplyr::mutate(CD8CCR6_corrected = `CD8+CCR6+` - b_intercept)

CD8CCR6_plot <- ggplot(data_adj, aes(Visit, CD8CCR6_corrected, color = `HIV Status`)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, position = position_dodge(width = 0.7)) +
  geom_jitter(alpha = 0.4, 
              size = 7, 
              position = position_jitterdodge(jitter.width=0.15, dodge.width=0.7)
  ) +
  labs(#title = "CD4 proportion (subject-corrected) by Visit and HIV Status",
    y = "Proportion of CD8+CCR6+ T cells (%)") +
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  geom_pwc(method = 'wilcox.test',
           label = "{ifelse(p < 0.0001, 'p < 0.0001', sprintf('p = %.4f', p))}",
           label.size = 8,
           tip.length = 0.01,
           hide.ns = F)+
  #scale_y_continuous(limits = c(0,150), breaks=c(0,20,40,60,80,100,120))+
  theme_bw()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 50, face = "bold"),
        axis.ticks.length = unit(0.5,"cm"),
        legend.text = element_text(size=25),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 30))
print(CD8CCR6_plot)
ggsave(CD8CCR6_plot,filename="Mixed_Models/Figure 1/CD8CCR6_plot.png",
       width = 17,height = 9,dpi = 300)
######### Fitting model on MAIT cells proportion with HIV as a grouping factor #########
fit_MAIT <- lmer(`CD3+ MAIT cells` ~ Visit * `HIV Status` + (1 | PID), data = data)
summary(fit_MAIT)
########## Estimated marginal means (adjusted means) 
emm_MAIT <- emmeans(fit_MAIT, ~ Visit | `HIV Status`)
emm_df <- as.data.frame(emm_MAIT)
########## Pairwise comparisons of HIV groups at each. visit
contrast(emm_MAIT, method = "pairwise", by = "Visit", adjust = "holm")
########## Plot adjusted means
ggplot(emm_df, aes(Visit, emmean, color = `HIV Status`, group = `HIV Status`))+
  geom_point(position = position_dodge(width = 0.3), size = 3)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position = position_dodge(width = 0.3), width = 0.15)+
  labs(
    title = "Adjusted CD3+ MAIT cells Proportion by Visit and HIV Status",
    y= "Adjusted mean CD3+ MAIT cells Proportion (95% CI)"
  )+
  theme_classic()
########## subject-corrected boxplot by HIV status
ranefs <- ranef(fit_MAIT)$PID %>%
  rownames_to_column("PID") %>%
  rename(b_intercept = `(Intercept)`)

data_adj <- data %>%
  left_join(ranefs, by = "PID") %>%
  dplyr::mutate(MAIT_corrected = `CD3+ MAIT cells` - b_intercept)

MAIT_plot <- ggplot(data_adj, aes(Visit, MAIT_corrected, color = `HIV Status`)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, position = position_dodge(width = 0.7)) +
  geom_jitter(alpha = 0.4, 
              size = 7, 
              position = position_jitterdodge(jitter.width=0.15, dodge.width=0.7)
  ) +
  labs(#title = "CD4 proportion (subject-corrected) by Visit and HIV Status",
    y = "Proportion of CD3+ MAIT cells (%)") +
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  geom_pwc(method = 'wilcox.test',
           label = "{ifelse(p < 0.0001, 'p < 0.0001', sprintf('p = %.4f', p))}",
           label.size = 8,
           tip.length = 0.01,
           hide.ns = F)+
  #scale_y_continuous(limits = c(0,150), breaks=c(0,20,40,60,80,100,120))+
  theme_bw()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 50, face = "bold"),
        axis.ticks.length = unit(0.5,"cm"),
        legend.text = element_text(size=25),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 30))
print(MAIT_plot)
ggsave(MAIT_plot,filename="Mixed_Models/Figure 1/MAIT_plot.png",
       width = 17,height = 9,dpi = 300)
######### Fitting model on Gamma Delta T cell proportion with HIV as a grouping factor #########
fit_dg <- lmer(`CD3+ gd T cells` ~ Visit * `HIV Status` + (1 | PID), data = data)
summary(fit_dg)
########## Estimated marginal means (adjusted means) 
emm_gd <- emmeans(fit_dg, ~ Visit | `HIV Status`)
emm_df <- as.data.frame(emm_gd)
########## Pairwise comparisons of HIV groups at each. visit
contrast(emm_gd, method = "pairwise", by = "Visit", adjust = "holm")
########## Plot adjusted means
ggplot(emm_df, aes(Visit, emmean, color = `HIV Status`, group = `HIV Status`))+
  geom_point(position = position_dodge(width = 0.3), size = 3)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position = position_dodge(width = 0.3), width = 0.15)+
  labs(
    title = "Adjusted CD3+ gamma Delta T cells Proportion by Visit and HIV Status",
    y= "Adjusted mean CD3+ gamma Delta T cells Proportion (95% CI)"
  )+
  theme_classic()
########## subject-corrected boxplot by HIV status
ranefs <- ranef(fit_dg)$PID %>%
  rownames_to_column("PID") %>%
  rename(b_intercept = `(Intercept)`)

data_adj <- data %>%
  left_join(ranefs, by = "PID") %>%
  dplyr::mutate(gd_corrected = `CD3+ gd T cells` - b_intercept)

gd_plot <- ggplot(data_adj, aes(Visit, gd_corrected, color = `HIV Status`)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, position = position_dodge(width = 0.7)) +
  geom_jitter(alpha = 0.4, 
              size = 7, 
              position = position_jitterdodge(jitter.width=0.15, dodge.width=0.7)
  ) +
  labs(#title = "CD4 proportion (subject-corrected) by Visit and HIV Status",
    y = "Proportion of CD3+ Gamma Delta cells (%)") +
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  geom_pwc(method = 'wilcox.test',
           label = "{ifelse(p < 0.0001, 'p < 0.0001', sprintf('p = %.4f', p))}",
           label.size = 8,
           tip.length = 0.01,
           hide.ns = T)+
  #scale_y_continuous(limits = c(0,150), breaks=c(0,20,40,60,80,100,120))+
  theme_bw()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 50, face = "bold"),
        axis.ticks.length = unit(0.5,"cm"),
        legend.text = element_text(size=25),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 30))
print(gd_plot)
ggsave(gd_plot,filename="Mixed_Models/Figure 1/gd_plot.png",
       width = 17,height = 9,dpi = 300)
######### Fitting model on CD3+ CD56+ T cell proportion with HIV as a grouping factor #########
fit_nkt <- lmer(`CD3+ CD56+ T cells` ~ Visit * `HIV Status` + (1 | PID), data = data)
summary(fit_nkt)
########## Estimated marginal means (adjusted means) 
emm_nkt <- emmeans(fit_nkt, ~ Visit | `HIV Status`)
emm_df <- as.data.frame(emm_nkt)
########## Pairwise comparisons of HIV groups at each. visit
contrast(emm_nkt, method = "pairwise", by = "Visit", adjust = "holm")
########## Plot adjusted means
ggplot(emm_df, aes(Visit, emmean, color = `HIV Status`, group = `HIV Status`))+
  geom_point(position = position_dodge(width = 0.3), size = 3)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position = position_dodge(width = 0.3), width = 0.15)+
  labs(
    title = "Adjusted CD3+ CD56+ T cells Proportion by Visit and HIV Status",
    y= "Adjusted mean CD3+ CD56+ T cells Proportion (95% CI)"
  )+
  theme_classic()
########## subject-corrected boxplot by HIV status
ranefs <- ranef(fit_nkt)$PID %>%
  rownames_to_column("PID") %>%
  rename(b_intercept = `(Intercept)`)

data_adj <- data %>%
  left_join(ranefs, by = "PID") %>%
  dplyr::mutate(nkt_corrected = `CD3+ CD56+ T cells` - b_intercept)

nkt_plot <- ggplot(data_adj, aes(Visit, nkt_corrected, color = `HIV Status`)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, position = position_dodge(width = 0.7)) +
  geom_jitter(alpha = 0.4, 
              size = 7, 
              position = position_jitterdodge(jitter.width=0.15, dodge.width=0.7)
  ) +
  labs(#title = "CD4 proportion (subject-corrected) by Visit and HIV Status",
    y = "Proportion of CD3+ CD56+ T cells (%)") +
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  geom_pwc(method = 'wilcox.test',
           label = "{ifelse(p < 0.0001, 'p < 0.0001', sprintf('p = %.4f', p))}",
           label.size = 8,
           tip.length = 0.01,
           hide.ns = T)+
  #scale_y_continuous(limits = c(0,150), breaks=c(0,20,40,60,80,100,120))+
  theme_bw()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 50, face = "bold"),
        axis.ticks.length = unit(0.5,"cm"),
        legend.text = element_text(size=25),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 30))
print(nkt_plot)
ggsave(nkt_plot,filename="Mixed_Models/Figure 1/nkt_plot.png",
       width = 17,height = 9,dpi = 300)
######### Fitting model on CD4+CXCR3-CCR6+ T cell proportion with HIV as a grouping factor #########
fit_4p3n6p <- lmer(`CD4+CXCR3-CCR6+` ~ Visit * `HIV Status` + (1 | PID), data = data)
summary(fit_4p3n6p)
########## Estimated marginal means (adjusted means) 
emm_4p3n6p <- emmeans(fit_4p3n6p, ~ Visit | `HIV Status`)
emm_df <- as.data.frame(emm_4p3n6p)
########## Pairwise comparisons of HIV groups at each. visit
contrast(emm_4p3n6p, method = "pairwise", by = "Visit", adjust = "holm")
########## Plot adjusted means
ggplot(emm_df, aes(Visit, emmean, color = `HIV Status`, group = `HIV Status`))+
  geom_point(position = position_dodge(width = 0.3), size = 3)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position = position_dodge(width = 0.3), width = 0.15)+
  labs(
    title = "Adjusted CD4+CXCR3-CCR6+ T cells Proportion by Visit and HIV Status",
    y= "Adjusted mean CD4+CXCR3-CCR6+ T cells Proportion (95% CI)"
  )+
  theme_classic()
########## subject-corrected boxplot by HIV status
ranefs <- ranef(fit_4p3n6p)$PID %>%
  rownames_to_column("PID") %>%
  rename(b_intercept = `(Intercept)`)

data_adj <- data %>%
  left_join(ranefs, by = "PID") %>%
  dplyr::mutate(CD4p3n6p_corrected = `CD4+CXCR3-CCR6+` - b_intercept)

CD4p3n6p_plot <- ggplot(data_adj, aes(Visit, CD4p3n6p_corrected, color = `HIV Status`)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, position = position_dodge(width = 0.7)) +
  geom_jitter(alpha = 0.4, 
              size = 7, 
              position = position_jitterdodge(jitter.width=0.15, dodge.width=0.7)
  ) +
  labs(#title = "CD4 proportion (subject-corrected) by Visit and HIV Status",
    y = "Proportion of CD4+CXCR3-CCR6+ T cells (%)") +
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  geom_pwc(method = 'wilcox.test',
           label = "{ifelse(p < 0.0001, 'p < 0.0001', sprintf('p = %.4f', p))}",
           label.size = 8,
           tip.length = 0.01,
           hide.ns = T)+
  #scale_y_continuous(limits = c(0,150), breaks=c(0,20,40,60,80,100,120))+
  theme_bw()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 50, face = "bold"),
        axis.ticks.length = unit(0.5,"cm"),
        legend.text = element_text(size=25),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 30))
print(CD4p3n6p_plot)
ggsave(CD4p3n6p_plot,filename="Mixed_Models/Figure 1/CD4p3n6p_plot.png",
       width = 17,height = 9,dpi = 300)
######### Fitting model on CD4+CXCR3+CCR6+ T cell proportion with HIV as a grouping factor #########
fit_4p3p6p <- lmer(`CD4+CXCR3+CCR6+` ~ Visit * `HIV Status` + (1 | PID), data = data)
summary(fit_4p3p6p)
########## Estimated marginal means (adjusted means) 
emm_4p3p6p <- emmeans(fit_4p3p6p, ~ Visit | `HIV Status`)
emm_df <- as.data.frame(emm_4p3p6p)
########## Pairwise comparisons of HIV groups at each. visit
contrast(emm_4p3p6p, method = "pairwise", by = "Visit", adjust = "holm")
########## Plot adjusted means
ggplot(emm_df, aes(Visit, emmean, color = `HIV Status`, group = `HIV Status`))+
  geom_point(position = position_dodge(width = 0.3), size = 3)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position = position_dodge(width = 0.3), width = 0.15)+
  labs(
    title = "Adjusted CD4+CXCR3+CCR6+ T cells Proportion by Visit and HIV Status",
    y= "Adjusted mean CD4+CXCR3+CCR6+ T cells Proportion (95% CI)"
  )+
  theme_classic()
########## subject-corrected boxplot by HIV status
ranefs <- ranef(fit_4p3p6p)$PID %>%
  rownames_to_column("PID") %>%
  rename(b_intercept = `(Intercept)`)

data_adj <- data %>%
  left_join(ranefs, by = "PID") %>%
  dplyr::mutate(CD4p3p6p_corrected = `CD4+CXCR3+CCR6+` - b_intercept)

CD4p3p6p_plot <- ggplot(data_adj, aes(Visit, CD4p3p6p_corrected, color = `HIV Status`)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, position = position_dodge(width = 0.7)) +
  geom_jitter(alpha = 0.4, 
              size = 7, 
              position = position_jitterdodge(jitter.width=0.15, dodge.width=0.7)
  ) +
  labs(#title = "CD4 proportion (subject-corrected) by Visit and HIV Status",
    y = "Proportion of CD4+CXCR3+CCR6+ T cells (%)") +
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  geom_pwc(method = 'wilcox.test',
           label = "{ifelse(p < 0.0001, 'p < 0.0001', sprintf('p = %.4f', p))}",
           label.size = 8,
           tip.length = 0.01,
           hide.ns = T)+
  #scale_y_continuous(limits = c(0,150), breaks=c(0,20,40,60,80,100,120))+
  theme_bw()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 50, face = "bold"),
        axis.ticks.length = unit(0.5,"cm"),
        legend.text = element_text(size=25),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 30))
print(CD4p3p6p_plot)
ggsave(CD4p3p6p_plot,filename="Mixed_Models/Figure 1/CD4p3p6p_plot.png",
       width = 17,height = 9,dpi = 300)
