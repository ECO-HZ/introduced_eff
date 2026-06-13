# Loading packages
library(openxlsx)
library(car)
library(lme4)
library(lmerTest) 
library(emmeans)   
library(multcomp) 
library(ggplot2)
library(glmmTMB)
library(tibble)
library(dplyr)
library(patchwork)
library(ggtext)

# load data
persistence_all_data = read.xlsx("field experiment cover biomass survival 260409.xlsx", sheet = "surv_correct", rowNames = F, colNames = T)
colnames(persistence_all_data)

id_vars <- setdiff(colnames(persistence_all_data), c("pot", "cover_mass", "Time", "species"))

persistence_all_data$Time <- as.character(persistence_all_data$Time)

persistence_wide <- persistence_all_data[,c("pot", "cover_mass", "Time", "species")] %>%
  tidyr::pivot_wider(
    id_cols = c("pot", "species"),
    names_from = Time,
    values_from = cover_mass,
    names_prefix = "cover_mass|")

persistence_wide$present_y1 = ifelse((persistence_wide$`cover_mass|1` + persistence_wide$`cover_mass|5`) > 0, 1, 0)
persistence_wide$present_y2 = ifelse((persistence_wide$`cover_mass|12` + persistence_wide$`cover_mass|14` + persistence_wide$`cover_mass|17`) > 0, 1, 0)
persistence_wide$present_y3 = ifelse((persistence_wide$`cover_mass|23` + persistence_wide$`cover_mass|25` + persistence_wide$`cover_mass|29`) > 0, 1, 0)
persistence_wide$present_y4 = ifelse((persistence_wide$`cover_mass|36` + persistence_wide$`cover_mass|39` + persistence_wide$`cover_mass|40`) > 0, 1, 0)

# Add community species treatment information.
colnames(persistence_all_data)
pot_infor = unique(persistence_all_data[,c(1:3,19:47)])
colnames(pot_infor)
dim(pot_infor)
#unique(pot_infor$pot)

# Add species identity treatment information.
species_infor = read.xlsx("field experiment cover biomass survival 260409.xlsx", sheet = "traits_mean", rowNames = F, colNames = T)
species_infor = species_infor[,c("species", "Latin_name", "type")]

persistence_wide_add = persistence_wide %>% left_join(pot_infor, by = "pot") %>% 
  left_join(species_infor, by = "species")

persistence_wide_add$Origin2 = as.factor(persistence_wide_add$Origin2)
persistence_wide_add$SR2 = as.factor(persistence_wide_add$SR2)
persistence_wide_add$pot = as.factor(persistence_wide_add$pot)
persistence_wide_add$add_sp_list = as.factor(persistence_wide_add$add_sp_list)
persistence_wide_add$res_sp_list = as.factor(persistence_wide_add$res_sp_list)

# 
colnames(persistence_wide_add)
persistence_wide_add_long = persistence_wide_add[,c(1:2,14:50)]
persistence_wide_add_long <- persistence_wide_add_long %>%
  tidyr::pivot_longer(
    cols = c("present_y1", "present_y2", "present_y3", "present_y4"), 
    names_to = "year",
    values_to = "present")

persistence_wide_add_long$Time[persistence_wide_add_long$year == "present_y1"] <- "5"
persistence_wide_add_long$Time[persistence_wide_add_long$year == "present_y2"] <- "17"
persistence_wide_add_long$Time[persistence_wide_add_long$year == "present_y3"] <- "29"
persistence_wide_add_long$Time[persistence_wide_add_long$year == "present_y4"] <- "40"
persistence_wide_add_long$Time = as.factor(persistence_wide_add_long$Time)

################################################################################
################################# Figure S2 ####################################
################################################################################

# Effects on the persistence probability of resident species
persistence_wide_R = subset(persistence_wide_add_long, type == "R")

model_all_y_R = glmer(present ~ species * Time + 
                        (1|Origin2:SR2:pot) + 
                        (1|res_sp_list) + 
                        (1|add_sp_list),
                      control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                      #control = glmmTMBControl(optimizer = "nlminb", optCtrl = list(iter.max = 1000, eval.max = 2000)),
                      family = binomial(link = "logit"), data = persistence_wide_R)

car::Anova(model_all_y_R, type = 2)
surv_all_y_eff_R = as_tibble(emmeans(model_all_y_R, specs = ~ Time|species, type = "response"))

# add y0 data
model_y0_data = data.frame(Time = "0", species = unique(surv_all_y_eff_R$species),
                           prob = 1,SE = 0,df = "Inf", asymp.LCL = 1,asymp.UCL = 1)

surv_all_y_eff_R = rbind(model_y0_data, surv_all_y_eff_R)

surv_all_y_eff_R$Time = factor(surv_all_y_eff_R$Time, levels = c("0","5","17","29","40"))

ggplot(data=surv_all_y_eff_R, aes(x = Time, y = prob, group=species, color=species)) + 
  geom_line(size=0.8) + 
  ggrepel::geom_text_repel(data = subset(surv_all_y_eff_R, Time == 40),
                           aes(label = species, color = species),direction = "y", hjust = 0, nudge_x = 0.5,segment.size = 0.2, 
                           #segment.color = "inherit",  # 连接线颜色
                           box.padding = 0.3,  point.padding = 0.3,min.segment.length = 0,size = 3.5,show.legend = FALSE) +
  theme_classic() + 
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title = element_text(colour='black', size=13),
        axis.text.x = element_text(colour='black',size=11),
        axis.text.y = element_text(colour='black',size=11),
        plot.tag = element_text(size = 14, face = "bold"),
        plot.title = ggtext::element_textbox(
          size = 14, color = "black", fill = "#E6E5E5",
          box.color = "black",padding = margin(5, 5, 5, 5), margin = margin(b = 0),       
          halign = 0.5, width = grid::unit(1, "npc")), #r = unit(3, "pt")     
        legend.position = "none") + 
  scale_color_manual(values = c("SV" = "#92BEA7", "PL" = "#A13E38", "MD" = "#D7B91B", "PP" = "#6699CC",
                                "PF" = "#7B7776", "MC" = "#3A55A4", "EP" = "#557278", "EC" = "#653028",
                                "SOR" = "#C15924", "EI" = "#637A34", "UL" = "#BFC525", "LH" = "#AC871F")) + 
  scale_x_discrete(breaks = unique(surv_all_y_eff_R$Time),
                   #labels = c("Starting", "2021", "2022", "2023", "2024"),
                   labels = c("0" = "Starting", "5" = "5", "17" = "17", "29" = "29", "40" = "40"),  
                   expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Sampling time (months)", 
       y = "Proportion of plots persisting\nacross four years",
       #title = "Resident species",tag = "A",
       color = "Species") -> Figure_S2; Figure_S2

################################################################################
################################# Figure 2A ####################################
################################################################################

# Effects on the persistence probability of added species
persistence_wide_X2 = subset(persistence_wide_add_long, type != "R")

model_all_y_X2 = glmer(present ~ species * Time + 
                         (1|Origin2:SR2:pot) + 
                         (1|res_sp_list) + 
                         (1|add_sp_list),
                       control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                       #control = glmmTMBControl(optimizer = "nlminb", optCtrl = list(iter.max = 1000, eval.max = 2000)),
                       family = binomial(link = "logit"), data = persistence_wide_X2)

car::Anova(model_all_y_X2, type = 2)
surv_all_y_eff_X2 = as_tibble(emmeans(model_all_y_X2, specs = ~ Time|species, type = "response"))

# add y0 data
model_y0_data = data.frame(Time = "0", species = unique(surv_all_y_eff_X2$species),
                           prob = 1,SE = 0,df = "Inf", asymp.LCL = 1,asymp.UCL = 1)

surv_all_y_eff_X2 = rbind(model_y0_data, surv_all_y_eff_X2)

surv_all_y_eff_X2$Time = factor(surv_all_y_eff_X2$Time, levels = c("0","5","17","29","40"))

ggplot(data=surv_all_y_eff_X2, aes(x = Time, y = prob, group=species, color=species)) + 
  geom_line(size=0.8) + 
  ggrepel::geom_text_repel(data = subset(surv_all_y_eff_X2, Time == 40),
                           aes(label = species, color = species),direction = "y", hjust = 0, nudge_x = 0.5,segment.size = 0.2, 
                           #segment.color = "inherit",  # 连接线颜色
                           box.padding = 0.3,  point.padding = 0.3,min.segment.length = 0,size = 3.5,show.legend = FALSE) +
  theme_classic() + 
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title = element_text(colour='black', size=13),
        axis.text.x = element_text(colour='black',size=11),
        axis.text.y = element_text(colour='black',size=11),
        plot.tag = element_text(size = 14, face = "bold"),
        plot.title = ggtext::element_textbox(
          size = 14, color = "black", fill = "#E6E5E5",
          box.color = "black",padding = margin(5, 5, 5, 5), margin = margin(b = 0),       
          halign = 0.5, width = grid::unit(1, "npc")), #r = unit(3, "pt")     
        legend.position = "none") + 
  scale_color_manual(values = c( "AA" = "#117732", "BB" = "#989932", "BF" = "#43AA9A", "CAL"  = "#CD6677", "CAR" = "#872254",
                                 "DA" = "#AB4499", "PA" = "#33228A", "SN" = "#6699CC", "SOC" = "#88CCEE", "ST"  = "#DDCB76")) + 
  scale_fill_manual(values = c( "AA" = "#117732", "BB" = "#989932", "BF" = "#43AA9A", "CAL"  = "#CD6677", "CAR" = "#872254",
                                "DA" = "#AB4499", "PA" = "#33228A", "SN" = "#6699CC", "SOC" = "#88CCEE", "ST"  = "#DDCB76")) + 
  scale_x_discrete(breaks = unique(surv_all_y_eff_X2$Time),
                   #labels = c("Starting", "2021", "2022", "2023", "2024"),
                   labels = c("0" = "Starting", "5" = "5", "17" = "17", "29" = "29", "40" = "40"),  
                   expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Sampling time (months)", 
       y = "Proportion of plots persisting\nacross four years",
       #title = "Introduced species",
       tag = "A", 
       color = "Species") -> Figure_2A; Figure_2A

################################################################################
######################### Insets in (Figure 2A)  ###############################
################################################################################

# select the 4th year data (for E2 AND N2 species)
persistence_4y_X2 = subset(persistence_wide_add, type != "R")
persistence_4y_X2$Origin2 = as.factor(persistence_4y_X2$Origin2)
persistence_4y_X2$SR2 = as.factor(persistence_4y_X2$SR2)
persistence_4y_X2$pot = as.factor(persistence_4y_X2$pot)
persistence_4y_X2$species = as.factor(persistence_4y_X2$species)
persistence_4y_X2$res_sp_list = as.factor(persistence_4y_X2$res_sp_list)
persistence_4y_X2$add_sp_list = as.factor(persistence_4y_X2$add_sp_list)

unique(persistence_4y_X2$species)
mod_present = glmer(present_y4 ~ SR2 * Origin2 + 
                      (1|species) + 
                      (1|Origin2:SR2:pot) + 
                      (1|res_sp_list) + 
                      (1|add_sp_list),
                    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                    #control = glmerControl(optimizer = nlminb, optCtrl = list(iter.max = 1000, eval.max = 2000)),
                    family = binomial(link = "logit"), data = persistence_4y_X2)
summary(mod_present)
car::Anova(mod_present, type = 2)

# effect of origin on the persistence of species added in the second planting
second_sp_origin = emmeans(mod_present, specs = pairwise ~ Origin2, type = "response")
second_sp_origin_data = as_tibble(second_sp_origin$emmeans)

ggplot(data = second_sp_origin_data, aes(x = Origin2, y = prob, fill = Origin2, color = Origin2)) + 
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0, size = 0.6) +
  geom_point(size = 3.5, pch = 21) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,1)) + 
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        legend.text = element_text(colour='black',size=11),
        legend.title = element_text(colour='black',size=12),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title = element_text(colour='black', size=13),
        axis.text.x = element_text(colour='black',size=11, angle = 35, vjust = 1, hjust = 1),
        #axis.text.x = element_text(colour='black',size=11),
        axis.text.y = element_text(colour='black',size=11),
        plot.tag = element_text(size = 14, face = "bold"),
        legend.position = "none") + 
  #scale_color_manual(values = c("Control" = "#5C6572", "Native" = "#36A9E1", "Exotic" = "#EB5B25")) + 
  scale_color_manual(values = c("Control" = "#70319D", "Native" = "#3C8FB6", "Exotic" = "#A9405F")) +
  scale_fill_manual(values = c("Control" = "#70319D", "Native" = "#3C8FB6", "Exotic" = "#A9405F")) +
  labs(x = NULL, y = NULL)

################################################################################
################################# Figure 4A  ###################################
################################################################################

traits_mean = read.xlsx("field experiment cover biomass survival 260409.xlsx", sheet = "traits_mean", rowNames = F, colNames = T)
traits_mean$F_flowT = as.numeric(traits_mean$F_flowT)

persistence_4y_X2 = persistence_4y_X2 %>% left_join(traits_mean)
library(effects)

colnames(persistence_4y_X2_multi)
persistence_4y_X2_multi <- persistence_4y_X2
persistence_4y_X2_multi$F_hgt <- log10(persistence_4y_X2_multi$F_hgt)
persistence_4y_X2_multi$F_LA <- log10(persistence_4y_X2_multi$F_LA)
persistence_4y_X2_multi$F_SLA <- log10(persistence_4y_X2_multi$F_SLA)
persistence_4y_X2_multi$F_LDMC <- log10(persistence_4y_X2_multi$F_LDMC)
persistence_4y_X2_multi$F_seeds <- log10(persistence_4y_X2_multi$F_seeds)

# Consider normalizing your data
var_select <- c("F_hgt","F_LA","F_SLA","F_LDMC","F_seeds","F_germT","F_flowT")
pd_attributes_variable <- attributes(scale(persistence_4y_X2_multi[var_select]))
total_data <- persistence_4y_X2_multi
total_data[var_select] <- scale(total_data[var_select])

cor.test(total_data$F_SLA, total_data$F_LDMC)

library(corrplot)
library(Hmisc)
corr_matrix <- rcorr(as.matrix((total_data[, c("F_hgt","F_LA","F_SLA","F_LDMC","F_seeds","F_germT","F_flowT")])), type = 'spearman')
corr_matrix$r   
corr_matrix$P  

p.mat = corr_matrix$P
diag(p.mat) = 0

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(corr_matrix$r, p.mat = p.mat, sig.level = 0.05, insig = 'blank', method = 'number',type = 'lower',
         diag = F, col=col(200), tl.cex = 0.8,tl.col = "black", number.cex = 0.8, order = "original",tl.srt = 45)

corrplot(corr_matrix$r, p.mat = p.mat, sig.level = 0.05, insig = 'blank', method = 'square',
         add = TRUE, type = 'lower', diag = F, col=col(200), tl.pos = 'n', cl.pos = 'n',outline = F, order = "original",
         addCoef.col = "black", number.cex = 0.8)

# performed the global linear regression model
total_data = subset(total_data, total_data$F_LA != "NA")

options(na.action = "na.fail")
fm1 <- glmer(present_y4 ~ F_hgt + F_LA + F_LDMC + F_germT + F_flowT + F_seeds + 
               (1|Origin2:SR2:pot) + 
               (1|res_sp_list) + 
               (1 |add_sp_list), 
             control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
             family = binomial(link = "logit"), data = total_data)
vif(fm1)
summary(fm1)
performance::r2(fm1)
MuMIn::r.squaredGLMM(fm1)

library(MuMIn)
dd12 <- dredge(fm1, trace = 2, rank = "AICc")
                 
# get the coefficient and the standard error
resultModel <- summary(object = MuMIn::model.avg(object = dd12, 
                                                 subset = delta < 4))$coefmat.subset[c(-1),]

resultModel <- as.data.frame(resultModel)
colnames(resultModel)[c(1,2)] <- c("coeff", "stdError")
rownames(resultModel)[which(rownames(resultModel) == "F_flowT")] <- "Flower data"
rownames(resultModel)[which(rownames(resultModel) == "F_germT")] <- "Germ data"
rownames(resultModel)[which(rownames(resultModel) == "F_hgt")] <- "Height"
rownames(resultModel)[which(rownames(resultModel) == "F_LDMC")] <- "LDMC"
rownames(resultModel)[which(rownames(resultModel) == "F_seeds")] <- "# of seeds"
rownames(resultModel)[which(rownames(resultModel) == "F_LA")] <- "LA"

resultModel$label <- rownames(resultModel)
resultModel$label <- factor(resultModel$label, levels = rev(c("Height", "LA", "LDMC", "Slope", 
                                                          "Germ data", "Flower data", "# of seeds")))

resultModel$explained_var = abs(resultModel$coeff)/sum(abs(resultModel$coeff))*100

ggplot(resultModel, aes(x=label, y=coeff,  color = label)) + 
  geom_errorbar(aes(ymin=coeff-1.96*stdError, ymax=coeff+1.96*stdError, color = label), 
                width=0, size = 0.8)+
  geom_point(size = 3.5, pch = 16) + 
  labs(x = '', 
       y = 'Standard regression coefficients') +  
  theme_minimal() + coord_flip() +
  scale_x_discrete(position = "top") + 
  scale_color_manual(values = c("# of seeds" = "#BC5546", "Flower data" = "#2F4590", "Germ data" = "#6EA3C5",
                                "LDMC" = "#8E333A", "LA" = "#EFA961", "Height" = "#769D89")) +
  scale_fill_manual(values = c("# of seeds" = "#BC5546", "Flower data" = "#2F4590", "Germ data" = "#6EA3C5",
                               "LDMC" = "#8E333A", "LA" = "#EFA961", "Height" = "#769D89")) +
  theme(axis.text = element_text(color = "black", size = 12),
        axis.title =  element_text(color = "black", size = 14),
        legend.text = element_text(size = 9, color = "black"),
        #plot.title = element_textbox(size = 12, color = "black", fill = "white",     
        #                             box.color = "black", width = grid::unit(1, "npc"),padding = margin(5, 5, 5, 5),  
        #                             margin = margin(b = 5), halign = 0.5,linetype = "solid"),
        #plot.margin = margin(0.5,1.5,0.5,1.5, unit = "cm"),
        axis.line.x = element_line(color = "black"),
        axis.ticks.x = element_line(color='black'),
        legend.position = 'none',
        panel.grid = element_blank(), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(), 
        plot.tag = element_text(size = 14, face = "bold")) +
  geom_hline(yintercept = 0,linetype=2) -> Fig4_ab
  

resultModel$label <- factor(resultModel$label, levels = (c("Height", "LA", "LDMC", 
                                                              "Germ data", "Flower data", "# of seeds")))

ggplot(resultModel, aes(x = "Variables", y = explained_var, fill = label)) +
  geom_col(width = 1, color = "white") +
  #scale_fill_viridis(option = "D", direction = -1) + # + 
  theme_classic()+ 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = c("# of seeds" = "#BC5546", "Flower data" = "#2F4590", "Germ data" = "#6EA3C5",
                              "LDMC" = "#8E333A", "LA" = "#EFA961", "Height" = "#769D89")) +
  theme(panel.grid = element_blank(), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(), 
        axis.title =  element_text(color = "black", size = 14),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        axis.text = element_text(size = 12, color = "black"),
        legend.position = "none") + 
  labs(x = '', y = "Relative effect of estimates") -> Fig4_aa

library(patchwork)
(Fig4_aa|Fig4_ab) + plot_layout(widths = c(0.3,0.7))

