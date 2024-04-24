#fixing the matrix package to go back to 6.1.2... nothing will work if it is updated
db <- available.packages(repos = "https://cran.r-project.org/")
tools::package_dependencies("Matrix", db = db, which = "LinkingTo", reverse = 
                              TRUE)[[1L]]
#load packeges
library(lme4)
library(plyr)
library(psych)
library(lmerTest)

library(TMB)
library(glmmTMB)

library(tidyr)
library(tidyverse)
library(dplyr)
library (haven)
library (rio)

library (sjstats)

library (bmlm) #to easily get the between and within effects
library (MuMIn) # to get R2 for linear models

#read in the data and call it data
#data <- read_csv(ENTER PATH)

#Recoding arm and burst as factors
data$arm<-as.factor(data$arm)
data$arm = relevel(factor(data$arm), ref="0")   ## CC as usual = reference

data$burst<-as.factor(data$burst)
data$burst = relevel(factor(data$burst), ref="1") ## pre-intervention = reference

###############Mindfulness analyses#####################################

#Calculating the random slopes

#I need to pull out just baseline, then mid intervention then post intervention
data_pre <- data %>% filter(burst == 1)

#Null model
mod1 <- glmmTMB(maas_avg~stressfulevent + (1| day) +  (1| base_Id),
                data = data_pre, family = Gamma(link = "log"))

summary(mod1)

#with random slope
mod2 <- glmmTMB(maas_avg ~stressfulevent + (1| day) + (1 + stressfulevent| base_Id),data = data_pre, family = Gamma(link = "log"))

summary(mod2)

#Print baseline random effects
random_effect_pre <- ranef(mod2)
print (random_effect_pre)

#Creating the data frame with the baseline random slope 
df_pre <- as.data.frame(random_effect_pre[["cond"]][["base_Id"]])
df_pre$base_Id<- rownames(random_effect_pre[["cond"]][["base_Id"]])
names(df_pre)[names(df_pre) == "stressfulevent"] <- "maas_slope_1"
names(df_pre)[names(df_pre) == "(Intercept)"] <- "maas_intercept_1"

#Mid intervention
data_mid <- data %>% filter(burst == 2)

mod3 <- glmmTMB(maas_avg ~stressfulevent + (1| day) + (1 + stressfulevent| base_Id),data = data_mid, family = Gamma(link = "log"))

summary(mod3)

#Print mid random effects
random_effect_mid <- ranef(mod3)
print (random_effect_mid)

#Creating the data frame with the mid-intervention random slope 
df_mid <- as.data.frame(random_effect_mid[["cond"]][["base_Id"]])
df_mid$base_Id<- rownames(random_effect_mid[["cond"]][["base_Id"]])
names(df_mid)[names(df_mid) == "stressfulevent"] <- "maas_slope_2"
names(df_mid)[names(df_mid) == "(Intercept)"] <- "maas_intercept_2"

#Post intervention
data_post <- data %>% filter(burst == 3)

mod4 <- glmmTMB(maas_avg ~stressfulevent + (1| day) + (1 + stressfulevent| base_Id),data = data_post, family = Gamma(link = "log"))
summary(mod4)

#Printpost random effects
random_effect_post <- ranef(mod4)
print (random_effect_post)

#Creating the data frame with the post-intervention random slope 
df_post <- as.data.frame(random_effect_post[["cond"]][["base_Id"]])
df_post$base_Id<- rownames(random_effect_post[["cond"]][["base_Id"]])
names(df_post)[names(df_post) == "stressfulevent"] <- "maas_slope_3"
names(df_post)[names(df_post) == "(Intercept)"] <- "maas_intercept_3"

#merging random effects all together
merged_random_premid <- df_pre %>%
  left_join(df_mid, by = c("base_Id"))
#merging data file
merged_re <- merged_random_premid  %>%
  left_join(df_post, by = c("base_Id"))

names(merged_re)

#merging all
long_merged_re <- merged_re %>%
  pivot_longer(cols = starts_with("maas_"), 
               names_to = c("variable", "burst"), 
               names_pattern = "maas_(\\w+)_(\\d+)$") %>%
  pivot_wider(names_from = "variable", values_from = "value")

# View the long format data frame
print(long_merged_re)

names(long_merged_re)[names(long_merged_re) == "intercept"] <- "maas_intercept"
names(long_merged_re)[names(long_merged_re) == "slope"] <- "maas_slope"

#must specfy burst as the same type.
long_merged_re$burst<-as.factor(long_merged_re$burst)
long_merged_re$base_Id<-as.double(long_merged_re$base_Id)

#merging data file
merged_data <- data %>%
  left_join(long_merged_re, by = c("base_Id", "burst"))

#moving variables up
merged_data <-(merged_data %>%
                 relocate (maas_slope, .after=burst))

merged_data <-(merged_data %>%
                 relocate (maas_intercept, .after=burst))

#Recode -999 to NA.. AKA recoding missing- nothing should be -999 but to
merged_data[merged_data=="-999"] <-NA

#Recoding arm and burst as factors
merged_data$arm<-as.factor(merged_data$arm)
merged_data$arm = relevel(factor(merged_data$arm), ref="0")  # CC as usual = reference

#making sure slope is a numeric value
merged_data$maas_slope<-as.numeric(merged_data$maas_slope)

#armXburst predicting the random slope of life stressors predicting mindfulness
#remember that baseline is the reference group here

mod6 <- lmer(maas_slope ~ arm + burst + arm*burst + (1 |base_Id), REML= F, data = merged_data)
summary(mod6)

#Now calculating the effect size 

#run the similar model without the interaction
mod6_noint <- lmer(maas_slope ~ arm + burst + (1 |base_Id), REML= F, data = merged_data)

summary(mod6_noint)

R2_mod6<-as.data.frame(r.squaredGLMM(mod6))
R2_mod6_null<-as.data.frame(r.squaredGLMM(mod6_noint))

f2_maas <- ((R2_mod6$R2m-R2_mod6_null$R2m)/(1-R2_mod6$R2m))
f2_maas

# Create Estimated Marginal Means for the model
library(emmeans)
emm_options(pbkrtest.limit = 3178)

emm <- emmeans(mod6, ~ arm:burst)
summary(emm)

# Extract the estimated marginal means as a data frame
emm_means_mindfulness <- as.data.frame(emm)

# Plotting the bar chart
p <- ggplot(emm_means_mindfulness, aes(x = factor(burst), y = emmean, fill = factor(arm))) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_manual(values = c("0" = "grey50", "1" = "black"), name = "Arm", labels = c("Mentoring-alone", "MBI+mentoring")) +
  scale_x_discrete(name = "Burst", labels = c("1" = "Baseline", "2" = "Mid-intervention", "3" = "Post-intervention")) +
  labs(
    title = "EMM Random Slope of Life Stressors predicting Mindlessness (Attention)
            by Burst and Arm",
    x = "Burst",
    y = "EMM of Random Slopes of Life Stressors 
    on Mindlessness (Attention)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(color = "black"),
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom"  # Move legend to the bottom
  )

p

library(gridExtra) #for adding table into figure

annotation_data <- as.data.frame(emm_means_mindfulness) %>%
  select(burst, arm, Estimate = emmean, SE, lower.CL, upper.CL)

# Formatting 95% CI as a single column
annotation_data$`95% CI` <- sprintf("(%0.2f, %0.2f)", round(c(annotation_data$lower.CL),2), round(c(annotation_data$upper.CL),2))

annotation_data$Estimate <- round(annotation_data$Estimate,3)
annotation_data$SE<- round(annotation_data$SE,2)

# Creating a table with EMM, SE, and CI values
table <- tableGrob(annotation_data, rows = NULL)

# Creating a table with modified column names and formatting
table <- tableGrob(annotation_data[, c("burst", "arm", "Estimate", "SE", "95% CI")], rows = NULL)

# Arranging the plot and table side by side
grid.arrange(p, table, ncol = 2, widths = c(3, 1))

# Conduct pairwise comparisons
emm_pairs <- pairs(emm)
emm_pairs <- as.data.frame(emm_pairs)

# Display the results
summary(emm_pairs)

#now looking to add info about mid vs. post intervention. given that burst is dummy coded- we need to run the model again to get all of the important information
merged_data$burst = relevel(factor(merged_data$burst), ref="2")

mod7 <- lmer(maas_slope ~ arm + burst + arm*burst + (1 |base_Id), REML= F, data = merged_data) 

summary(mod7)

###The same code listed above was used to investigate mindful non-judgment (scs2) and difficulites with emotion regulation (ders_avg_r)


##################Missing Data Exploration######################
pacman::p_load(
  rio,           # import/export
  tidyverse,     # data mgmt and viz
  naniar,        # assess and visualize missingness
  mice           # missing data imputation
)
#percent missing of all data in data frame- .43% missing across entire data frame
myvars_key <- c("maas_avg", "ders_avg_r", "scs2","stressfulevent","arm",
                "burst")
newdata_key<- data[myvars_key]
pct_miss(newdata_key)
pct_miss(newdata_key$maas_avg)
pct_miss(newdata_key$scs2)
pct_miss(newdata_key$ders_avg_r)
pct_miss(newdata_key$stressfulevent)

