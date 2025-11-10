library(ggpubr)
library(tidyverse)
Nasomune <- read_csv("Data/Main_Files_Thesis/Samples_Included_SPN_PLHIV_Paper.csv")
Masterfile <- read_csv("Results/Masterfile.csv")


Files_in_Paper <- Masterfile %>%
  filter(`LAB ID` %in% Nasomune$`LAB ID`)

Fig1a <- Files_in_Paper %>%
  dplyr::filter(
    `Epithelial count`>100,
    `CD45+ count`>200,
    `Visit`=="Week 1",
    `HIV Status`!="HIV-"
  ) %>%
  dplyr::mutate(Viral_load_Cat = case_when(
    grepl("[0-9]", HIV_VL) ~ "Detected",
          TRUE ~ "Undetected"
    )
  ) %>%
  ggplot(aes(y=log10(NER), x=factor(Viral_load_Cat,
                                    levels = c("Undetected","Detected")), 
             color = factor(Viral_load_Cat,
                            levels = c("Undetected","Detected"))))+
  geom_boxplot(width=0.5,notch = F, outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.2),size=7) +
  geom_pwc(method = 'wilcox.test',
           label = "{ifelse(p < 0.0001, 'p < 0.0001', sprintf('p = %.4f', p))}",
           label.size = 10,
           tip.length = 0.01,
           hide.ns = F)+
  labs(y="Neutrophil to epithelial cell ratio",
       x="Viremic status",
       title = "a.")+
  scale_color_manual(values = c('#E97132','#00B0F0'))+
  theme_bw()+
  #facet_wrap(~`HIV Status`)+
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
Fig1a
ggsave(Fig1a,filename="New_Figures/Reviewers_Comments/Fig1a.png",
       width = 8,height = 10,dpi = 300)
``````
Fig1b <- Files_in_Paper %>%
  dplyr::filter(
    `Epithelial count`>100,
    `CD45+ count`>200,
    `Visit`=="Week 1",
    `HIV Status`!="HIV-"
  ) %>%
  dplyr::mutate(Viral_load_Cat = case_when(
    grepl("[0-9]", HIV_VL) ~ "Detected",
    TRUE ~ "Undetected"
  )
  ) %>%
  ggplot(aes(y=log10(MER), x=factor(Viral_load_Cat,
                                    levels = c("Undetected","Detected")), 
             color = factor(Viral_load_Cat,
                            levels = c("Undetected","Detected"))))+
  geom_boxplot(width=0.5,notch = F, outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.2),size=7) +
  geom_pwc(method = 'wilcox.test',
           label = "{ifelse(p < 0.0001, 'p < 0.0001', sprintf('p = %.4f', p))}",
           label.size = 10,
           tip.length = 0.01,
           hide.ns = F)+
  labs(y="Monocyte to epithelial cell ratio",
       x="Viremic status",
       title = "b.")+
  scale_color_manual(values = c('#E97132','#00B0F0'))+
  #facet_wrap(~`HIV Status`)+
  theme_bw()+
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
ggsave(Fig1b,filename="New_Figures/Reviewers_Comments/Fig1b.png",
       width = 8,height = 10,dpi = 300)
``````
ggsave(filename = "New_Figures/Reviewers_Comments/Fig1.png",
       plot = ((Fig1a|Fig1b|plot_layout(ncol = 3,nrow = 1, width = c(1,1,1))))/ 
         (plot_spacer()|plot_layout(width = c(1,1,1)))/
         (plot_spacer()|plot_layout(width = c(1))),
       width = 32, height = 29, units = "in", dpi = 300)

````````

Fig1c <- Files_in_Paper %>%
  dplyr::filter(
    `Epithelial count`>100,
    `CD45+ count`>200,
    `Visit`=="Week 1",
    `HIV Status`!="HIV-"
  ) %>%
  ggplot(aes(y=Abs_CD4_Count, x=log10(NER)))+
  geom_point(size=7, alpha=0.7, aes(color=`HIV Status`))+
  geom_smooth(method = 'lm', color='black', size=2)+
  stat_cor(method = "spearman",
             size = 10)+
  #facet_wrap(~`HIV Status`)+
  scale_color_manual(values = c('#941100','#005493'))+
  labs(y="CD4 count (cells/ul)",
       x="Neutrophil to epithelial cell ratio",
       title = "a.")+
  theme_bw()+
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
ggsave(Fig1c,filename="New_Figures/Reviewers_Comments/Fig1c.png",
       width = 8,height = 10,dpi = 300)


````````
Fig1d <- Files_in_Paper %>%
  dplyr::filter(
    `Epithelial count`>100,
    `CD45+ count`>200,
    `Visit`=="Week 1",
    `HIV Status`!="HIV-"
  ) %>%
  ggplot(aes(y=Abs_CD4_Count, x=log10(MER)))+
  geom_point(size=7, alpha=0.7, aes(color=`HIV Status`))+
  geom_smooth(method = 'lm', color='black', size=2)+
  stat_cor(method = "spearman",
           size = 10)+
  #facet_wrap(~`HIV Status`)+
  scale_color_manual(values = c('#941100','#005493'))+
  labs(y="CD4 count (cells/ul)",
       x="Monocyte to epithelial cell ratio",
       title = "b.")+
  theme_bw()+
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
ggsave(Fig1d,filename="New_Figures/Reviewers_Comments/Fig1d.png",
       width = 8,height = 10,dpi = 300)

ggsave(filename = "New_Figures/Reviewers_Comments/Fig2.png",
       plot = ((Fig1c|Fig1d|plot_layout(ncol = 3,nrow = 1, width = c(1,1,1))))/ 
         (plot_spacer()|plot_layout(width = c(1,1,1)))/
         (plot_spacer()|plot_layout(width = c(1))),
       width = 32, height = 29, units = "in", dpi = 300)
````````

Fig1e <- Files_in_Paper %>%
  dplyr::filter(
    `Epithelial count`>100,
    `CD45+ count`>200,
    `Visit`=="Week 1",
    `HIV Status`!="HIV-"
  ) %>%
  ggplot(aes(y=`ART duration`, x=log10(NER)))+
  geom_point(size=7, aes(color=`HIV Status`), alpha=0.7)+
  geom_smooth(method = 'lm', color='black', size=2)+
  stat_cor(method = "spearman",
           size = 10)+
  #facet_wrap(~`HIV Status`)+
  scale_color_manual(values = c('#941100','#005493'))+
  labs(y="ART duration (months)",
       x="Neutrophil to epithelial cell ratio",
       title = "a.")+
  theme_bw()+
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
ggsave(Fig1e,filename="New_Figures/Reviewers_Comments/Fig1e.png",
       width = 8,height = 10,dpi = 300)

``````

Fig1f <- Files_in_Paper %>%
  dplyr::filter(`Epithelial count`>100,
                `CD45+ count`>200,
                `Visit`=="Week 1",
                `HIV Status`!="HIV-") %>%
  ggplot(aes(y=`ART duration`, x=log10(MER)))+
  geom_point(size=7, aes(color=`HIV Status`), alpha=0.7)+
  geom_smooth(method = 'lm', color='black', size=2)+
  stat_cor(method = "spearman",
           size = 10)+
  #facet_wrap(~`HIV Status`)+
  scale_color_manual(values = c('#941100','#005493'))+
  labs(y="ART duration (months)",
       x="Monocyte to epithelial cell ratio",
       title = "b.")+
  theme_bw()+
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
ggsave(Fig1f,filename="New_Figures/Reviewers_Comments/Fig1f.png",
       width = 8,height = 10,dpi = 300)
``````

ggsave(filename = "New_Figures/Reviewers_Comments/Fig3.png",
       plot = ((Fig1e|Fig1f|plot_layout(ncol = 3,nrow = 1, width = c(1,1,1))))/ 
         (plot_spacer()|plot_layout(width = c(1,1,1)))/
         (plot_spacer()|plot_layout(width = c(1))),
       width = 32, height = 29, units = "in", dpi = 300)

``````
# Logistic regression
# ---- Packages ----
# Install once if needed:
install.packages(c("tidyverse","broom","car","pROC","yardstick","rsample","gt","performance"))

library(tidyverse)
library(broom)
library(car)         # VIF
library(pROC)        # ROC/AUC
library(yardstick)   # metrics
library(rsample)     # train/test split (optional)
library(gt)          # pretty tables
library(performance) # quick diagnostics

# ---- Load your data ----
# df <- read.csv("your_file.csv")
# For this script we assume df already exists in the workspace.

# Column names expected:
# "LAB ID", "ART duration", "Carriage Status", "NER", "MER",
# "HIV Status", "HIV duration", "HIV_VL"

# ---- Rename columns to syntactically safe names ----
Regression_data <- Files_in_Paper %>%
  dplyr::select(
    `LAB ID`,
    `ART duration`,
    `Carriage Status`,
    `HIV Status`,
    `HIV duration`,
    NER,
    MER,
    Abs_CD4_Count,
    HIV_VL
  ) %>%
  dplyr::rename(
    LAB_ID        = `LAB ID`,
    ART_duration  = `ART duration`,
    CarriageStatus= `Carriage Status`,
    HIV_Status    = `HIV Status`,
    HIV_duration  = `HIV duration`
  ) %>%
  dplyr::mutate(ifelse(grepl()))
  dplyr::filter(NER>0, MER>0)

# ---- Basic cleaning / types ----
Regression_data_clean <- Regression_data %>%
  dplyr::mutate
  dplyr::mutate(
    # Outcome: ensure binary factor with "Carriage Positive" as the event level
    CarriageStatus = factor(CarriageStatus),
    CarriageStatus = forcats::fct_relevel(CarriageStatus, "Carriage Positive"),
    
    # Predictor: HIV_Status as a factor with HIV- as reference
    HIV_Status = factor(HIV_Status),
    HIV_Status = forcats::fct_relevel(HIV_Status, "HIV-"),
    
    # Continuous predictors (ensure numeric)
    ART_duration = as.numeric(ART_duration),
    HIV_duration = as.numeric(HIV_duration),
    NER          = as.numeric(NER),
    MER          = as.numeric(MER)
  ) %>%
  # Drop rows with missing values in variables used (simple approach)
  drop_na(CarriageStatus, HIV_Status, ART_duration, HIV_duration, HIV_VL, NER, MER)

# ---- (Optional) Train/Test split for honest evaluation ----
set.seed(123)
split <- initial_split(df_clean, prop = 0.8, strata = CarriageStatus)
train <- training(split)
test  <- testing(split)

# ---- Fit logistic regression (main effects) ----
# NOTE: glm uses the first level of the outcome as "success".
# We set "Carriage Positive" first above, so the model estimates
# log-odds of being Carriage Positive.
fit <- glm(
  CarriageStatus ~ ART_duration + HIV_duration + HIV_VL + NER + MER + HIV_Status,
  data = train,
  family = binomial(link = "logit")
)

summary(fit)

# ---- Multicollinearity check ----
car::vif(fit)

# ---- Odds ratios with 95% CIs ----
or_tbl <- broom::tidy(fit, conf.int = TRUE, conf.level = 0.95, exponentiate = TRUE) %>%
  mutate(
    term = recode(term,
                  "(Intercept)" = "Intercept",
                  "ART_duration" = "ART duration (per unit)",
                  "HIV_duration" = "HIV duration (per unit)",
                  "HIV_VL"       = "HIV_VL (per unit)",
                  "NER"          = "NER (per unit)",
                  "MER"          = "MER (per unit)"
    )
  ) %>%
  arrange(desc(estimate))

or_tbl %>%
  mutate(across(c(estimate, conf.low, conf.high, p.value), ~round(.x, 3))) %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  gt::gt() %>%
  gt::fmt_markdown(columns = term) %>%
  gt::tab_header(title = "Logistic Regression: Odds Ratios (Carriage Positive as event)",
                 subtitle = "Exponential of coefficients with 95% CIs")

# ---- Quick model diagnostics ----
performance::check_model(fit)  # residuals, linearity, outliers (opens plots)

# ---- Discrimination: ROC/AUC on test set ----
test$pred_prob <- predict(fit, newdata = test, type = "response")

roc_obj <- pROC::roc(response = test$CarriageStatus,
                     predictor = test$pred_prob,
                     levels = levels(test$CarriageStatus),
                     direction = ">")  # higher prob = more likely Positive
auc_val <- pROC::auc(roc_obj)
print(roc_obj)
print(paste("AUC:", round(as.numeric(auc_val), 3)))

# ---- Choose a threshold (Youden's J) and compute confusion metrics ----
opt <- pROC::coords(roc_obj, x = "best", best.method = "youden", ret = c("threshold","sensitivity","specificity"))
thr <- as.numeric(opt["threshold"])
thr

test$pred_class <- ifelse(test$pred_prob >= thr, "Carriage Positive", "Carriage Negative") %>% factor(levels = levels(test$CarriageStatus))

# Accuracy, sensitivity, specificity, etc.
conf_mat <- yardstick::conf_mat(test, truth = CarriageStatus, estimate = pred_class)
conf_mat
yardstick::summary(conf_mat, event_level = "first")  # event is "Carriage Positive"

# ---- Calibration (how well predicted probs match observed) ----
# Group into bins and compare mean predicted vs. observed rates
calib_df <- test %>%
  mutate(bin = ntile(pred_prob, 10)) %>%
  group_by(bin) %>%
  summarise(
    mean_pred = mean(pred_prob),
    obs_rate  = mean(CarriageStatus == "Carriage Positive"),
    n = n(),
    .groups = "drop"
  )

print(calib_df)

# ---- Optional: Standardize continuous predictors for comparability ----
train_std <- train %>%
  mutate(across(c(ART_duration, HIV_duration, HIV_VL, NER, MER), scale))

fit_std <- glm(
  CarriageStatus ~ ART_duration + HIV_duration + HIV_VL + NER + MER + HIV_Status,
  data = train_std,
  family = binomial()
)

broom::tidy(fit_std, exponentiate = TRUE, conf.int = TRUE) %>%
  arrange(desc(estimate)) %>%
  print(n = Inf)

# ---- Optional: Check nonlinearity with restricted cubic splines (rms) ----
# install.packages("rms")
# library(rms)
# dd <- datadist(train); options(datadist = "dd")
# fit_spline <- lrm(CarriageStatus ~ rcs(ART_duration, 4) + rcs(HIV_duration, 4) +
#                                  rcs(HIV_VL, 4) + rcs(NER, 4) + rcs(MER, 4) + HIV_Status,
#                   data = train, x = TRUE, y = TRUE)
# anova(fit_spline)
# plot(Predict(fit_spline, ART_duration))

# ---- Optional: Add domain-driven interactions (example) ----
# fit_int <- glm(
#   CarriageStatus ~ ART_duration + HIV_duration + HIV_VL + NER + MER + HIV_Status +
#                    HIV_Status:HIV_VL,
#   data = train, family = binomial()
# )
# broom::tidy(fit_int, exponentiate = TRUE, conf.int = TRUE)
