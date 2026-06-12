# 
library(ggplot2)
traits_thmeme = theme_classic() + 
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title = element_text(colour='black', size=13),
        axis.text.x = element_text(colour='black',size=11, angle = 35, vjust = 1, hjust = 1),
        axis.text.y = element_text(colour='black',size=11),
        plot.tag = element_text(size = 14, face = "bold"),
        #plot.margin = margin(1,1,2,1,unit="cm"),
        legend.position = "none")


# Effects of species' origin and species identity on the traits measured in situ
# Height
field_traits1 = read.xlsx("field traits 260409.xlsx", sheet = "Sheet2", rowNames = F, colNames = T)
field_traits1$Origin2 = ifelse(field_traits1$Origin2 == "Exotic", "Alien", field_traits1$Origin2)
unique(field_traits1$Origin2)
field_mean_Height_sp = Rmisc::summarySE(subset(field_traits1, Time %in% c(5,17,29,40)), 
                                        groupvars = c("species"), measurevar = "Height", na.rm = TRUE)

# X2 species
field_traits_X2 = subset(field_traits1, type  != "R" & Time %in% c(5,17,29,40))
unique(field_traits_X2$species)
dim(field_traits_X2)
field_traits_X2$Height = log10(field_traits_X2$Height)
field_traits_X2$species = factor(field_traits_X2$species, levels = c("BF", "BB", "PA", "SN", "SOC", "ST", "CAR", "AA", "DA", "CAL"))

# Height
library(emmeans)
mod_Height1 = lmer(Height ~ Origin2 + SR2 + 
                     (1|Time) + 
                     (1|species) + 
                     (1|Origin2:SR2:pot) + 
                     (1|res_sp_list) + 
                     (1|add_sp_list), data = field_traits_X2)
summary(mod_Height1)
shapiro.test(residuals(mod_Height1))
car::Anova(mod_Height1)
effect1 = emmeans(mod_Height1, pairwise ~ Origin2, type = "response")
emm1_multi1 = as.tibble(multcomp::cld(effect1, alpha=0.05, Letters=LETTERS, adjust = "none", decreasing = T))
emm1_multi1$.group <- trimws(emm1_multi1$.group)
colnames(emm1_multi1)[1] = "Group"

ggplot(emm1_multi1, aes(x = Group, y = emmean, fill = Group, color = Group)) + 
  #geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = emmean - 1.96*SE, ymax = emmean + 1.96*SE), width = 0) + 
  geom_point(size = 3.5, pch = 21) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  scale_color_manual(values = c("Native" = "#3C8FB6", "Alien" = "#A9405F")) + 
  scale_fill_manual(values = c("Native" = "#3C8FB6", "Alien" = "#A9405F")) + 
  traits_thmeme + 
  labs(x = NULL, 
       y = "Plant Height (cm)",
       #y = expression("Plant Height (cm, log"[10]*")"),
       tag = "A") -> p1; p1

# 与单株地上生物量的关系
field_traits_X2_filter = field_traits_X2 %>% tidyr::drop_na(Height, Rel_abun)
field_traits_X2_filter = subset(field_traits_X2_filter, Rel_abun != 0)
field_traits_X2_filter$Rel_abun = log10(field_traits_X2_filter$Rel_abun*100)
field_traits_X2_filter$Ind_mass = log10(field_traits_X2_filter$Ind_mass)
#field_traits_X2_filter$Height = log10(field_traits_X2_filter$Height)

mod_Height_lme = lmer(Rel_abun ~ Height + 
                        (1|Time) + 
                        #(1|species) + 
                        (1|Origin2:SR2:pot) + 
                        (1|res_sp_list) + 
                        (1|add_sp_list), data = field_traits_X2_filter)
#shapiro.test(residuals(mod_Height_lme))
summary(mod_Height_lme)
std_mod_Height_lme = effectsize::effectsize(mod_Height_lme)
r2_values <- MuMIn::r.squaredGLMM(mod_Height_lme)
r2_marginal <- sprintf("%.3f", r2_values[1])

library(AICcmodavg) 

field_traits_X2_filter$F0 = predictSE(mod_Height_lme, field_traits_X2_filter, level = 0)$fit
field_traits_X2_filter$SE <- predictSE(mod_Height_lme, field_traits_X2_filter, level = 0)$se.fit

ggplot(field_traits_X2_filter, aes(x=Height, y=Rel_abun)) +
  #geom_point(size = 2.5, pch = 21, fill = "grey80") + 
  geom_ribbon(fill = "#769D89", aes(x=Height,
                                    ymin = F0 - 1.96 * SE,
                                    ymax = F0 + 1.96 * SE), alpha = 0.2) + 
  geom_line(aes(y=F0), size=1.2, color = "#769D89") + 
  #geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = F) + 
  theme_classic() +
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title = element_text(colour='black', size=11),
        axis.text = element_text(colour='black',size=11),
        plot.tag = element_text(size = 14, face = "bold"),
        #plot.margin = margin(1,1,2,1,unit="cm"),
        legend.position = "none") +
  ylim(-3.5, 3.1) + 
  labs(x = expression(Plant~height~("cm,"~log[10])),
       y = expression("Relative abundance\nof introduced species (%, log"[10]*")"),
       tag = "A") -> pp1; pp1


# LA
field_traits2 = read.xlsx("field traits 260409.xlsx", sheet = "Sheet1", rowNames = F, colNames = T)
field_traits2$Origin2 = ifelse(field_traits2$Origin2 == "Exotic", "Alien", field_traits2$Origin2)

# X2 species
field_traits_X2 = subset(field_traits2, type != "R")
unique(field_traits_X2$species)
field_traits_X2$LA = log10(field_traits_X2$LA)

field_traits_X2$species = factor(field_traits_X2$species, levels = c("BF", "BB", "PA", "SN", "SOC", "ST", "CAR", "AA", "DA", "CAL"))

# LA
mod_LA1 = lmer(LA ~ Origin2 + 
                 (1|Origin2:SR2:pot) + 
                 (1|species) + 
                 (1|res_sp_list) + 
                 (1|add_sp_list), data = field_traits_X2)
summary(mod_LA1)
shapiro.test(residuals(mod_LA1))
car::Anova(mod_LA1)
effect1 = emmeans(mod_LA1, pairwise ~ Origin2, type = "response")
emm1_multi1 = as.tibble(multcomp::cld(effect1, alpha=0.05, Letters=LETTERS, adjust = "none", decreasing = T))
emm1_multi1$.group <- trimws(emm1_multi1$.group)
colnames(emm1_multi1)[1] = "Group"

ggplot(emm1_multi1, aes(x = Group, y = emmean, fill = Group, color = Group)) + 
  #geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0) + 
  geom_point(size = 3.5, pch = 21) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  scale_color_manual(values = c("Native" = "#3C8FB6", "Alien" = "#A9405F")) + 
  scale_fill_manual(values = c("Native" = "#3C8FB6", "Alien" = "#A9405F")) + 
  traits_thmeme + 
  labs(x = NULL, 
       y = expression("Leaf area ("*cm^2*")"),
       tag = "B") -> p2; p2

# 与单株地上生物量的关系
colnames(field_traits_X2)
field_traits_X2_filter = field_traits_X2 %>% tidyr::drop_na(LA, Rel_abun)
dim(field_traits_X2_filter)
field_traits_X2_filter = subset(field_traits_X2_filter, Rel_abun != 0)
field_traits_X2_filter$Rel_abun = log10(field_traits_X2_filter$Rel_abun*100)
#field_traits_X2_filter$LA = log10(field_traits_X2_filter$LA)

mod_LA_area_lme = lmer(Rel_abun ~ LA + 
                         #(1|species) + 
                         (1|Origin2:SR2:pot) + 
                         (1|res_sp_list) + 
                         (1|add_sp_list), data = field_traits_X2_filter)
#shapiro.test(residuals(mod_Height_lme))
summary(mod_LA_area_lme)
std_mod_LA_area_lme = effectsize::effectsize(mod_LA_area_lme)
r2_values <- MuMIn::r.squaredGLMM(mod_LA_area_lme)
r2_marginal <- sprintf("%.3f", r2_values[1])

field_traits_X2_filter$F0 = predictSE(mod_LA_area_lme, field_traits_X2_filter, level = 0)$fit
field_traits_X2_filter$SE <- predictSE(mod_LA_area_lme, field_traits_X2_filter, level = 0)$se.fit

ggplot(field_traits_X2_filter, aes(x=LA, y=Rel_abun, ymax=0.5, ymin=0)) +
  #geom_point(size = 2.5, pch = 21, fill = "grey80") + 
  geom_ribbon(fill = "#EFA961", aes(x=LA,
                                    ymin = F0 - 1.96 * SE,
                                    ymax = F0 + 1.96 * SE), alpha = 0.2) + 
  geom_line(aes(y=F0), size=1, color = "#EFA961") + 
  #geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = F) + 
  theme_classic() +
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title = element_text(colour='black', size=13),
        axis.text = element_text(colour='black',size=11),
        plot.tag = element_text(size = 14, face = "bold"),
        #plot.margin = margin(1,1,2,1,unit="cm"),
        legend.position = "none") +
  ylim(-3.5, 3.1) + 
  labs(x = expression(Leaf~area~(cm^2~", "~log[10])),
       y = expression("Relative abundance\nof introduced species (%, log"[10]*")"),
       tag = "B") -> pp2; pp2


# SLA
field_traits2 = read.xlsx("field traits 260409.xlsx", sheet = "Sheet1", rowNames = F, colNames = T)

field_traits2$Origin2 = ifelse(field_traits2$Origin2 == "Exotic", "Alien", field_traits2$Origin2)

# X2 species
field_traits_X2 = subset(field_traits2, type != "R")
unique(field_traits_X2$species)
field_traits_X2$SLA = log10(field_traits_X2$SLA)

field_traits_X2$species = factor(field_traits_X2$species, levels = c("BF", "BB", "PA", "SN", "SOC", "ST", "CAR", "AA", "DA", "CAL"))

# SLA
mod_SLA1 = lmer(SLA ~ Origin2 + 
                  (1|Origin2:SR2:pot) + 
                  (1|species) + 
                  (1|res_sp_list) + 
                  (1|add_sp_list), data = field_traits_X2)
summary(mod_SLA1)
#shapiro.test(residuals(mod_SLA1))
car::Anova(mod_SLA1)
effect1 = emmeans(mod_SLA1, pairwise ~ Origin2, type = "response")
emm1_multi1 = as.tibble(multcomp::cld(effect1, alpha=0.05, Letters=LETTERS, adjust = "none", decreasing = T))
emm1_multi1$.group <- trimws(emm1_multi1$.group)
colnames(emm1_multi1)[1] = "Group"

ggplot(emm1_multi1, aes(x = Group, y = emmean, fill = Group, color = Group)) + 
  #geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0) + 
  geom_point(size = 3.5, pch = 21) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  scale_color_manual(values = c("Native" = "#3C8FB6", "Alien" = "#A9405F")) + 
  scale_fill_manual(values = c("Native" = "#3C8FB6", "Alien" = "#A9405F")) + 
  traits_thmeme + 
  labs(x = NULL, 
       y = expression("Specific leaf area (cm"^2*" g"^{-1}*")"),
       tag = "C") -> p3; p3

# 与单株地上生物量的关系
field_traits_X2_filter = field_traits_X2 %>% tidyr::drop_na(SLA, Rel_abun)
field_traits_X2_filter = subset(field_traits_X2_filter, Rel_abun != 0)
field_traits_X2_filter$Rel_abun = log10(field_traits_X2_filter$Rel_abun*100)
#field_traits_X2_filter$SLA = log10(field_traits_X2_filter$SLA)

mod_SLA_lme = lmer(Rel_abun ~ SLA + 
                     #(1|species) + 
                     (1|Origin2:SR2:pot) + 
                     (1|res_sp_list) + 
                     (1|add_sp_list), data = field_traits_X2_filter)
#shapiro.test(residuals(mod_SLA_lme))
summary(mod_SLA_lme)

std_mod_SLA_lme = effectsize::effectsize(mod_SLA_lme)
r2_values <- MuMIn::r.squaredGLMM(mod_SLA_lme)
r2_marginal <- sprintf("%.3f", r2_values[1])

field_traits_X2_filter$F0 = predictSE(mod_SLA_lme, field_traits_X2_filter, level = 0)$fit
field_traits_X2_filter$SE = predictSE(mod_SLA_lme, field_traits_X2_filter, level = 0)$se.fit

ggplot(field_traits_X2_filter, aes(x=SLA, y=Rel_abun, ymax=0.5, ymin=0)) +
  #geom_point(size = 2.5, pch = 21, fill = "grey80") + 
  geom_ribbon(fill = "#EFF0F0", aes(x=SLA,
                                    ymin = F0 - 1.96 * SE,
                                    ymax = F0 + 1.96 * SE)) + 
  geom_line(aes(y=F0), size=1.2, linetype = 1) + 
  #geom_smooth(method = "lm", formula = y ~ x, se = F) + 
  theme_classic() +
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title = element_text(colour='black', size=13),
        axis.text = element_text(colour='black',size=11),
        plot.tag = element_text(size = 14, face = "bold"),
        #plot.margin = margin(1,1,2,1,unit="cm"),
        legend.position = "none") +
  ylim(-3.5, 3.1) + 
  labs(x = expression("Specific leaf area (cm"^2*" g"^{-1}*", log"[10]*")"),
       y = expression("Relative abundance\nof introduced species (%, log"[10]*")"),
       tag = "C") -> pp3; pp3


# LDMC
field_traits2 = read.xlsx("field traits 260409.xlsx", sheet = "Sheet1", rowNames = F, colNames = T)

field_traits2$Origin2 = ifelse(field_traits2$Origin2 == "Exotic", "Alien", field_traits2$Origin2)

# X2 species
field_traits_X2 = subset(field_traits2, type != "R")
unique(field_traits_X2$species)
field_traits_X2$LDMC = log10(field_traits_X2$LDMC)

field_traits_X2$species = factor(field_traits_X2$species, levels = c("BF", "BB", "PA", "SN", "SOC", "ST", "CAR", "AA", "DA", "CAL"))

mod_LDMC1 = lmer(LDMC ~ Origin2 + 
                   (1|Origin2:SR2:pot) + 
                   (1|species) + 
                   (1|res_sp_list) + 
                   (1|add_sp_list), data = field_traits_X2)
summary(mod_LDMC1)
#shapiro.test(residuals(mod_LDMC1))
car::Anova(mod_LDMC1)
effect1 = emmeans(mod_LDMC1, pairwise ~ Origin2, type = "response")
emm1_multi1 = as.tibble(multcomp::cld(effect1, alpha=0.05, Letters=LETTERS, adjust = "none", decreasing = T))
emm1_multi1$.group <- trimws(emm1_multi1$.group)
colnames(emm1_multi1)[1] = "Group"

ggplot(emm1_multi1, aes(x = Group, y = emmean, fill = Group, color = Group)) + 
  #geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0) + 
  geom_point(size = 3.5, pch = 21) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  scale_color_manual(values = c("Native" = "#3C8FB6", "Alien" = "#A9405F")) + 
  scale_fill_manual(values = c("Native" = "#3C8FB6", "Alien" = "#A9405F")) + 
  traits_thmeme + 
  labs(x = NULL, 
       #y = "Leaf dry matter content (g g-1)",
       y = expression("Leaf dry matter content (g g"^{-1}*")"),
       tag = "D") -> p4; p4

# 与单株地上生物量的关系
field_traits_X2_filter = field_traits_X2 %>% tidyr::drop_na(LDMC, Rel_abun)
field_traits_X2_filter = subset(field_traits_X2_filter, Rel_abun != 0)
field_traits_X2_filter$Rel_abun = log10(field_traits_X2_filter$Rel_abun*100)
#field_traits_X2_filter$LDMC = log10(field_traits_X2_filter$LDMC)

mod_LDMC_lme = lmer(Rel_abun ~ LDMC + 
                      #(1|species) + 
                      (1|Origin2:SR2:pot) + 
                      (1|res_sp_list) + 
                      (1|add_sp_list), data = field_traits_X2_filter)
shapiro.test(residuals(mod_LDMC_lme))
summary(mod_LDMC_lme)

std_mod_LDMC_lme = effectsize::effectsize(mod_LDMC_lme)
r2_values <- MuMIn::r.squaredGLMM(mod_LDMC_lme)
r2_marginal <- sprintf("%.3f", r2_values[1])

field_traits_X2_filter$F0 = predictSE(mod_LDMC_lme, field_traits_X2_filter, level = 0)$fit
field_traits_X2_filter$SE <- predictSE(mod_LDMC_lme, field_traits_X2_filter, level = 0)$se.fit

ggplot(field_traits_X2_filter, aes(x=LDMC, y=Rel_abun, ymax=0.5, ymin=0)) +
  #geom_point(size = 2.5, pch = 21, fill = "grey80") + 
  geom_ribbon(fill = "#8E333A", aes(x=LDMC,
                                    ymin = F0 - 1.96 * SE,
                                    ymax = F0 + 1.96 * SE), alpha=0.2) + 
  geom_line(aes(y=F0), size=1, linetype = 2, color = "#8E333A") + 
  #geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = F) + 
  theme_classic() +
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title = element_text(colour='black', size=13),
        axis.text = element_text(colour='black',size=11),
        plot.tag = element_text(size = 14, face = "bold"),
        #plot.margin = margin(1,1,2,1,unit="cm"),
        legend.position = "none") +
  ylim(-3.5, 3.1) + 
  labs(x = expression("Leaf dry matter content (g g"^{-1}*", log"[10]*")"),
       y = expression("Relative abundance\nof introduced species (%, log"[10]*")"),
       tag = "D") -> pp4; pp4


# germinate
field_traits2 = read.xlsx("field traits 260409.xlsx", sheet = "Sheet7", rowNames = F, colNames = T)

field_mean_germinate_sp = Rmisc::summarySE(subset(field_traits2, year %in% c(2024, 2025)), 
                                           groupvars = c("species"), measurevar = "germinate", na.rm = TRUE)
#write.xlsx(field_mean_germinate_sp, "field_mean_germinate_sp.xlsx")

field_traits2$Origin2 = ifelse(field_traits2$Origin2 == "Exotic", "Alien", field_traits2$Origin2)

# X2 species
field_traits_X2 = subset(field_traits2, type != "R")
unique(field_traits_X2$species)
#field_traits_X2$germinate = sqrt(field_traits_X2$germinate)

field_traits_X2$species = factor(field_traits_X2$species, levels = c("BF", "BB", "PA", "SN", "SOC", "ST", "CAR", "AA", "DA", "CAL"))

mod_germinate1 = lmer(germinate ~ Origin2 + 
                        (1|year) + 
                        (1|species) + 
                        (1|Origin2:SR2:pot) +
                        (1|res_sp_list) + 
                        (1|add_sp_list), data = field_traits_X2)
summary(mod_germinate1)
shapiro.test(residuals(mod_germinate1))
car::Anova(mod_germinate1)
effect1 = emmeans(mod_germinate1, pairwise ~ Origin2, type = "response")
emm1_multi1 = as.tibble(multcomp::cld(effect1, alpha=0.05, Letters=LETTERS, adjust = "none", decreasing = T))
emm1_multi1$.group <- trimws(emm1_multi1$.group)
colnames(emm1_multi1)[1] = "Group"

ggplot(emm1_multi1, aes(x = Group, y = emmean, fill = Group, color = Group)) + 
  #geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0) + 
  geom_point(size = 3.5, pch = 21) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  scale_color_manual(values = c("Native" = "#3C8FB6", "Alien" = "#A9405F")) + 
  scale_fill_manual(values = c("Native" = "#3C8FB6", "Alien" = "#A9405F")) + 
  traits_thmeme + 
  labs(x = NULL, 
       y = "First germination data (Julian)",
       tag = "E") -> p5; p5

# 与单株地上生物量的关系
field_traits_X2_filter = field_traits_X2 %>% tidyr::drop_na(germinate, Rel_abun)
field_traits_X2_filter = subset(field_traits_X2_filter, Rel_abun != 0)
field_traits_X2_filter$Rel_abun = log10(field_traits_X2_filter$Rel_abun*100)

mod_germinate_lme = lmer(Rel_abun ~ germinate + 
                           (1|year) +    
                           #(1|species) + 
                           (1|Origin2:SR2:pot) + 
                           (1|res_sp_list) + 
                           (1|add_sp_list), data = field_traits_X2_filter)
#shapiro.test(residuals(mod_germinate_lme))
summary(mod_germinate_lme)

std_mod_germinate_lme = effectsize::effectsize(mod_germinate_lme)
r2_values <- MuMIn::r.squaredGLMM(mod_germinate_lme)
r2_marginal <- sprintf("%.3f", r2_values[1])

field_traits_X2_filter$F0 = predictSE(mod_germinate_lme, field_traits_X2_filter, level = 0)$fit
field_traits_X2_filter$SE <- predictSE(mod_germinate_lme, field_traits_X2_filter, level = 0)$se.fit

ggplot(field_traits_X2_filter, aes(x=germinate, y=Rel_abun, ymax=0.5, ymin=0)) +
  #geom_point(size = 2.5, pch = 21, fill = "grey80") + 
  geom_ribbon(fill = "#6EA3C5", aes(x=germinate,
                                    ymin = F0 - 1.96 * SE,
                                    ymax = F0 + 1.96 * SE), alpha = 0.2) + 
  geom_line(aes(y=F0), size=1.2, linetype = 1, color = "#6EA3C5") + 
  #geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = F) + 
  theme_classic() +
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title = element_text(colour='black', size=13),
        axis.text = element_text(colour='black',size=11),
        plot.tag = element_text(size = 14, face = "bold"),
        #plot.margin = margin(1,1,2,1,unit="cm"),
        legend.position = "none") +
  ylim(-3.5, 3.1) + 
  labs(x = "First germination data (Julian)",
       y = expression("Relative abundance\nof introduced species (%, log"[10]*")"),
       tag = "E") -> pp5; pp5


# flowering
field_traits2 = read.xlsx("field traits 260409.xlsx", sheet = "Sheet7", rowNames = F, colNames = T)

field_mean_flowering_sp = Rmisc::summarySE(subset(field_traits2, year %in% c(2023, 2024)), 
                                           groupvars = c("species"), measurevar = "flowering", na.rm = TRUE)

field_traits2$Origin2 = ifelse(field_traits2$Origin2 == "Exotic", "Alien", field_traits2$Origin2)

# X2 species
field_traits_X2 = subset(field_traits2, type != "R")
unique(field_traits_X2$species)
#field_traits_X2$flowering = sqrt(field_traits_X2$flowering)

field_traits_X2$species = factor(field_traits_X2$species, levels = c("BF", "BB", "PA", "SN", "SOC", "ST", "CAR", "AA", "DA", "CAL"))

# flowering
mod_flowering1 = lmer(flowering ~ Origin2 + 
                        (1|year) + 
                        (1|species) + 
                        (1|Origin2:SR2:pot) +
                        (1|res_sp_list) + 
                        (1|add_sp_list), data = field_traits_X2)
summary(mod_flowering1)
shapiro.test(residuals(mod_flowering1))
car::Anova(mod_flowering1)
effect1 = emmeans(mod_flowering1, pairwise ~ Origin2, type = "response")
emm1_multi1 = as.tibble(multcomp::cld(effect1, alpha=0.05, Letters=LETTERS, adjust = "none", decreasing = T))
emm1_multi1$.group <- trimws(emm1_multi1$.group)
colnames(emm1_multi1)[1] = "Group"


ggplot(emm1_multi1, aes(x = Group, y = emmean, fill = Group, color = Group)) + 
  #geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0) + 
  geom_point(size = 3.5, pch = 21) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  scale_color_manual(values = c("Native" = "#3C8FB6", "Alien" = "#A9405F")) + 
  scale_fill_manual(values = c("Native" = "#3C8FB6", "Alien" = "#A9405F")) + 
  traits_thmeme + 
  labs(x = NULL, 
       y = "First flowering data (Julian)",
       tag = "F") -> p6; p6

# 与单株地上生物量的关系
field_traits_X2_filter = field_traits_X2 %>% tidyr::drop_na(flowering, Rel_abun)
field_traits_X2_filter = subset(field_traits_X2_filter, Rel_abun != 0)
field_traits_X2_filter$Rel_abun = log10(field_traits_X2_filter$Rel_abun*100)
#field_traits_X2_filter$flowering = log10(field_traits_X2_filter$flowering)

mod_flowering_lme = lmer(Rel_abun ~ flowering + 
                           (1|year) + 
                           #(1|species) + 
                           (1|Origin2:SR2:pot) + 
                           (1|res_sp_list) + 
                           (1|add_sp_list), data = field_traits_X2_filter)
#shapiro.test(residuals(mod_flowering_lme))
summary(mod_flowering_lme)
Anova(mod_flowering_lme)
std_mod_flowering_lme = effectsize::effectsize(mod_flowering_lme)
r2_values <- MuMIn::r.squaredGLMM(mod_flowering_lme)
r2_marginal <- sprintf("%.3f", r2_values[1])

field_traits_X2_filter$F0 = predictSE(mod_flowering_lme, field_traits_X2_filter, level = 0)$fit
field_traits_X2_filter$SE <- predictSE(mod_flowering_lme, field_traits_X2_filter, level = 0)$se.fit

ggplot(field_traits_X2_filter, aes(x=flowering, y=Rel_abun, ymax=0.5, ymin=0)) +
  #geom_point(size = 2.5, pch = 21, fill = "grey80") + 
  geom_ribbon(fill = "#2F4590", aes(x=flowering,
                                    ymin = F0 - 1.96 * SE,
                                    ymax = F0 + 1.96 * SE), alpha = 0.2) + 
  geom_line(aes(y=F0), size=1.2, linetype = 1, color = "#2F4590") + 
  #geom_smooth(method = "lm", formula = y ~ x, se = F) + 
  theme_classic() +
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title = element_text(colour='black', size=13),
        axis.text = element_text(colour='black',size=11),
        plot.tag = element_text(size = 14, face = "bold"),
        #plot.margin = margin(1,1,2,1,unit="cm"),
        legend.position = "none") +
  ylim(-3.5, 3.1) + 
  labs(x = "First flowering data (Julian)",
       y = expression("Relative abundance\nof introduced species (%, log"[10]*")"),
       tag = "F") -> pp6; pp6


# Seed_plant
field_traits2 = read.xlsx("field traits 260409.xlsx", sheet = "Sheet4", rowNames = F, colNames = T)

field_traits2$Origin2 = ifelse(field_traits2$Origin2 == "Exotic", "Alien", field_traits2$Origin2)

field_mean_seed_sp = Rmisc::summarySE(subset(field_traits2, Time %in% c(5,17,29,40)), 
                                      groupvars = c("species"), measurevar = "Seed_plant", na.rm = TRUE)
#write.xlsx(field_mean_seed_sp, "field_mean_seed_sp.xlsx")

# X2 species
field_traits_X2 = subset(field_traits2, type != "R")
unique(field_traits_X2$species)
field_traits_X2$Seed_plant = log10(field_traits_X2$Seed_plant)

field_traits_X2$species = factor(field_traits_X2$species, levels = c("BF", "BB", "PA", "SN", "SOC", "ST", "CAR", "AA", "DA", "CAL"))

# Seed_plant
mod_Seed_plant1 = lmer(Seed_plant ~ Origin2 + 
                         (1|Time) + 
                         (1|species) + 
                         (1|Origin2:SR2:pot) +
                         (1|res_sp_list) + 
                         (1|add_sp_list), data = field_traits_X2)
summary(mod_Seed_plant1)
shapiro.test(residuals(mod_Seed_plant1))
car::Anova(mod_Seed_plant1)
effect1 = emmeans(mod_Seed_plant1, pairwise ~ Origin2, type = "response")
emm1_multi1 = as.tibble(multcomp::cld(effect1, alpha=0.05, Letters=LETTERS, adjust = "none", decreasing = T))
emm1_multi1$.group <- trimws(emm1_multi1$.group)
colnames(emm1_multi1)[1] = "Group"

ggplot(emm1_multi1, aes(x = Group, y = emmean, fill = Group, color = Group)) + 
  #geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0) + 
  geom_point(size = 3.5, pch = 21) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  scale_color_manual(values = c("Native" = "#3C8FB6", "Alien" = "#A9405F")) + 
  scale_fill_manual(values = c("Native" = "#3C8FB6", "Alien" = "#A9405F")) + 
  traits_thmeme + 
  labs(x = NULL, 
       y = "Seed production per individual",
       tag = "G") -> p7; p7

# 与单株地上生物量的关系
field_traits_X2_filter = field_traits_X2 %>% tidyr::drop_na(Seed_plant, Rel_abun)
field_traits_X2_filter = subset(field_traits_X2_filter, Rel_abun != 0)
field_traits_X2_filter$Rel_abun = log10(field_traits_X2_filter$Rel_abun*100)
#field_traits_X2_filter$Seed_plant = log10(field_traits_X2_filter$Seed_plant)

mod_Seed_plant_lme = lmer(Rel_abun ~ Seed_plant + 
                            (1|Time) + 
                            #(1|species) + 
                            (1|Origin2:SR2:pot) + 
                            (1|res_sp_list) + 
                            (1|add_sp_list), data = field_traits_X2_filter)
summary(mod_Seed_plant_lme)

std_mod_Seed_plant_lme = effectsize::effectsize(mod_Seed_plant_lme)
# Standardization method: refit
#Parameter   | Std. Coef. |        95% CI
#----------------------------------------
#(Intercept) |       0.07 | [-0.39, 0.52]
#Seed_plant  |       0.47 | [ 0.42, 0.52]

r2_values <- MuMIn::r.squaredGLMM(mod_Seed_plant_lme)
r2_marginal <- sprintf("%.3f", r2_values[1])

field_traits_X2_filter$F0 = predictSE(mod_Seed_plant_lme, field_traits_X2_filter, level = 0)$fit
field_traits_X2_filter$SE <- predictSE(mod_Seed_plant_lme, field_traits_X2_filter, level = 0)$se.fit

ggplot(field_traits_X2_filter, aes(x=Seed_plant, y=Rel_abun, ymax=0.5, ymin=0)) +
  #geom_point(size = 2.5, pch = 21, fill = "grey80") + 
  geom_ribbon(fill = "#BC5546", aes(x=Seed_plant,
                                    ymin = F0 - 1.96 * SE,
                                    ymax = F0 + 1.96 * SE), alpha = 0.2) + 
  geom_line(aes(y=F0), size=1.2, linetype = 1, color = "#BC5546") + 
  #geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = F) + 
  theme_classic() +
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title = element_text(colour='black', size=13),
        axis.text = element_text(colour='black',size=11),
        plot.tag = element_text(size = 14, face = "bold"),
        #plot.margin = margin(1,1,2,1,unit="cm"),
        legend.position = "none") +
  ylim(-3.5, 3) + 
  labs(x = expression("Number of seed production (log"[10]*")"),
       y = expression("Relative abundance\nof introduced species (%, log"[10]*")"),
       tag = "G") -> pp7; pp7


(pp1/pp4/pp7)|(pp2/pp5/pp7)|(pp3/pp6/pp7) -> Figure_4B

(p1/p5)|(p2/p6)|(p3/p7)|(p4/p7) -> Figure_S5

#ggsave("Figure_S6-0416.pdf", plot = Figure_S5, width = 13, height = 11, units = "in", dpi = 300)

#ggsave("Figure_S7-.pdf", plot = Figure_S7, width = 13, height = 11, units = "in", dpi = 300)

