#--------------------------------------------------------------------------------
# Baleen SIA in fin and blue whale hybrids 
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# 05 - LMM
#--------------------------------------------------------------------------------

# Load packages

library(lme4)
library(lmerTest)
library(effects)
library(lvmisc)

setwd("/Users/marcruizisagales/Documents/GitHub/Blue-fin-hybrids-isotopes/2code")
source("03_Suess_Laws_effect_correction.R")
df

whale <- df[,c(3,34,5,7,19)] # select only necessary columns
whale1 <-whale[complete.cases(whale), ] # keep rows where all columns have non-missing values
whale1$species<- factor(whale1$species, levels = c('fin','hybrid','blue'))

#test for normality
whale1 %>% dplyr::group_by(species) %>% shapiro_test(d15n)
whale1 %>% dplyr::group_by(species) %>% shapiro_test(d13c)
whale1 %>% dplyr::group_by(species) %>% shapiro_test(d34s)

#test for homoscedasticity
leveneTest(d15n ~ species, data = whale1)
leveneTest(d13c ~ species, data = whale1)
leveneTest(d34s ~ species, data = whale1)

#d15n
lmer_dN <- lmer(d15n ~ 1 + species + (1|Whale), data=whale1)
summary(lmer_dN)
plot(allEffects(lmer_dN))

# Model validation
plot_model(lmer_dN) # diagnostics
cv_m1 <- lvmisc::loo_cv(lmer_dN, whale1, "Whale", keep = "used") # cross-validation
accuracy(lmer_dN) # accuracy
compare_accuracy(lmer_dN, cv_m1) # compare accuracy
confint(lmer_dN) # confint

# Model validation with boostraping
qqPlot(resid(lmer_dN)) # normality of the residuals
qqPlot(c(ranef(lmer_dN)$Whale$`(Intercept)`)) # normality of the random effects
plot(lmer_dN) # variance homogeneity

f = function(m) {res= fixef(m)[2]; names(res) = names(fixef(m)[2]); res}
boost <- lme4::bootMer(lmer_dN, f, nsim=200, type="semiparametric", use.u=TRUE)
plot(boost, index= 1)

pred_dN <- ggpredict(lmer_dN, terms = "species")
pred_dN$x <- factor(pred_dN$x, levels = c("blue", "hybrid", "fin"))

pN <- ggplot(pred_dN, aes(x = x, y = predicted)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high, color = x), size = 0.8) +
  scale_color_manual(values = species_colors) +  # your predefined palette
  labs(
    #title = expression(paste(delta^15, "N")), 
    x = "", 
    y = expression(delta^15*N~"(\u2030)")
  ) +
  theme_article(base_size = 15, base_family = "Optima") +
  theme(legend.position = "none", aspect.ratio = 1, axis.text.x = element_blank()) + ylim(7.2,12.4) 

pN

#d13c
lmer_dC <- lmer(d13c ~ 1 + species + (1|Whale), data=whale1)
summary(lmer_dC)
plot(allEffects(lmer_dC))

# Model validation
plot_model(lmer_dC) # diagnostics
cv_m1 <- lvmisc::loo_cv(lmer_dC, whale1, "Whale", keep = "used") # cross-validation
accuracy(lmer_dC) # accuracy
compare_accuracy(lmer_dC, cv_m1) # compare accuracy
confint(lmer_dC) # confint

# Model validation with boostraping
qqPlot(resid(lmer_dC)) # normality of the residuals
qqPlot(c(ranef(lmer_dC)$Whale$`(Intercept)`)) # normality of the random effects
plot(lmer_dC) # variance homogeneity

f = function(m) {res= fixef(m)[2]; names(res) = names(fixef(m)[2]); res}
boost <- lme4::bootMer(lmer_dC, f, nsim=200, type="semiparametric", use.u=TRUE)
plot(boost, index= 1)

pred_dC <- ggpredict(lmer_dC, terms = "species")
pred_dC$x <- factor(pred_dC$x, levels = c("blue", "hybrid", "fin"))

pC <- ggplot(pred_dC, aes(x = x, y = predicted)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high, color = x), size = 0.8) +
  scale_color_manual(values = species_colors) +  # your predefined palette
  labs(
    #title = expression(paste(delta^13, "C")), 
    x = "", 
    y = expression(delta^13*C~"(\u2030)")
  ) +
  theme_article(base_size = 15, base_family = "Optima") +
  theme(legend.position = "none", aspect.ratio = 1, axis.text.x = element_blank()) + ylim(-20.9,-18.1) 

pC

#d34s
lmer_dS <- lmer(d34s ~ 1 + species + (1|Whale), data=whale1)
summary(lmer_dS)
plot(allEffects(lmer_dS))

# Model validation
plot_model(lmer_dS) # diagnostics
cv_m1 <- loo_cv(lmer_dS, df, "Whale", keep = "used") # cross-validation
accuracy(lmer_dS) # accuracy
compare_accuracy(lmer_dS, cv_m1) # compare accuracy
confint(lmer_dS) # confint

# Model validation with boostraping
qqPlot(resid(lmer_dS)) # normality of the residuals
qqPlot(c(ranef(lmer_dS)$Whale$`(Intercept)`)) # normality of the random effects
plot(lmer_dS) # variance homogeneity

f = function(m) {res= fixef(m)[2]; names(res) = names(fixef(m)[2]); res}
boost <- lme4::bootMer(lmer_dS, f, nsim=200, type="semiparametric", use.u=TRUE)
plot(boost, index= 1)

pred_dS <- ggpredict(lmer_dS, terms = "species")
pred_dS$x <- factor(pred_dS$x, levels = c("blue", "hybrid", "fin"))

pS <- ggplot(pred_dS, aes(x = x, y = predicted)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high, color = x), size = 0.8) +
  scale_color_manual(values = species_colors) +  # your predefined palette
  labs(
    #title = expression(paste(delta^34, "S")), 
    x = "Species", 
    y = expression(delta^34*S~"(\u2030)")
  ) +
  theme_article(base_size = 15, base_family = "Optima") +
  theme(legend.position = "none", aspect.ratio = 1) + scale_x_discrete(labels = c("blue" = "Blue", "hybrid" = "Hybrid", "fin" = "Fin")) + ylim(16.1,19.9) 

pS

ptot <- pN / pC / pS
ptot

#ggsave("/Users/marcruizisagales/Documents/GitHub/Blue-fin-hybrids-isotopes/03figures/predictions_N_C_S.svg", plot = last_plot(), dpi = 300, width = 10, height = 15, units = "cm", device = "svg")

