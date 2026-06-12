# Effects of species' origin and species identity on the traits measured in situ

field_traits1 = read.xlsx("field traits.xlsx", sheet = "Sheet2", rowNames = F, colNames = T)
field_traits1$Origin2 = ifelse(field_traits1$Origin2 == "Exotic", "Alien", field_traits1$Origin2)
unique(field_traits1$Origin2)

# X2 species
field_traits_X2 = subset(field_traits1, type  != "R" & Time %in% c(5,17,29,40))
unique(field_traits_X2$species)
dim(field_traits_X2)

#field_traits_X2 = subset(field_traits_X2, Rel_abun != 0)
field_traits_X2$Rel_abun = (field_traits_X2$Rel_abun*100)
field_traits_X2_sum = Rmisc::summarySE(data = field_traits_X2, groupvars = c("species", "Time"),
                                       measurevar = "Rel_abun", na.rm = T)

field_traits_X2_order = Rmisc::summarySE(data = field_traits_X2, groupvars = c("species"),
                                         measurevar = "Rel_abun", na.rm = T) %>%
  arrange(Rel_abun)

field_traits_X2_sum$species = factor(field_traits_X2_sum$species, levels = field_traits_X2_order$species)
field_traits_X2_sum$Time = as.factor(field_traits_X2_sum$Time)

ggplot(field_traits_X2_sum, aes(x = Time, y = Rel_abun, group = species, fill = species, color = species)) + 
  #geom_errorbar(aes(ymin = Rel_abun-ci, ymax = Rel_abun+ci), 
  #              width = 0, position = position_dodge(width = 5), size = 0.8) + 
  geom_line(aes(color = species), size = 0.8) +
  #geom_point(size = 3, position = position_dodge(width = 5), pch = 21, color = "black") + 
  scale_x_discrete(breaks = unique(field_traits_X2_sum$Time),
                   #labels = c("Starting", "2021", "2022", "2023", "2024"),
                   labels = c("5" = "5", "17" = "17", "29" = "29", "40" = "40"),  
                   expand = expansion(mult = c(0.1, 0.1))) +
  #scale_y_continuous(limits = c(-0, 12),
  #                   expand = expansion(mult = 0),
  #                   breaks = seq(-0, 12, by = 2)) +
  theme_classic() +
  scale_color_manual(values = c( "AA" = "#117732", "BB" = "#989932", "BF" = "#43AA9A", "CAL"  = "#CD6677", "CAR" = "#872254",
                                 "DA" = "#AB4499", "PA" = "#33228A", "SN" = "#6699CC", "SOC" = "#88CCEE", "ST"  = "#DDCB76")) + 
  scale_fill_manual(values = c( "AA" = "#117732", "BB" = "#989932", "BF" = "#43AA9A", "CAL"  = "#CD6677", "CAR" = "#872254",
                                "DA" = "#AB4499", "PA" = "#33228A", "SN" = "#6699CC", "SOC" = "#88CCEE", "ST"  = "#DDCB76")) + 
  #scale_shape_manual(values = c(1:10)) + 
  #facet_wrap(~ SR2, labeller = labeller(SR2 = c("1" = "No. of species = 1", 
  #                                              "2" = "No. of species = 2", 
  #                                              "4" = "No. of species = 4"))) +
  ggrepel::geom_text_repel(data = subset(field_traits_X2_sum, Time == 40),
                           aes(label = species, color = species),direction = "y", hjust = 0, nudge_x = 0.5,segment.size = 0.2, 
                           #segment.color = "inherit",  # 连接线颜色
                           box.padding = 0.3,  point.padding = 0.3,min.segment.length = 0,size = 3.5,show.legend = FALSE) +
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
        strip.background = element_rect(color = "black", size = 0.5, fill = "white"),
        strip.text = element_text(colour='black', size=11), 
        #plot.margin = margin(1,1,1,2,unit="cm"),
        legend.position = "none") + 
  labs(x = "Sampling time (months)", tag = "B",
       y = expression("Relative abundance of species\nin introduced assemblages (%)"))->p1





# Resident species
field_traits_R = subset(field_traits1, type  == "R" & Time %in% c(5,17,29,40))
unique(field_traits_R$species)
dim(field_traits_R)

#field_traits_R_no0 = subset(field_traits_R, Rel_abun != 0)
field_traits_R$Rel_abun_log = (field_traits_R$Rel_abun*100)
field_traits_R_sum = Rmisc::summarySE(data = field_traits_R, groupvars = c("species", "Time"),
                                      measurevar = "Rel_abun_log", na.rm = T)

field_traits_R_order = Rmisc::summarySE(data = field_traits_R, groupvars = c("species"),
                                        measurevar = "Rel_abun_log", na.rm = T) %>%
  arrange(Rel_abun_log)

field_traits_R_sum$species = factor(field_traits_R_sum$species, levels = field_traits_R_order$species)
field_traits_R_sum$Time = as.factor(field_traits_R_sum$Time)

ggplot(field_traits_R_sum, aes(x = Time, y = Rel_abun_log, group = species, fill = species, color = species)) + 
  #geom_errorbar(aes(ymin = Rel_abun_log-ci, ymax = Rel_abun_log+ci), 
  #              width = 0, position = position_dodge(width = 5), size = 0.8) + 
  geom_line(aes(color = species), size = 0.8) +
  #geom_point(size = 3, position = position_dodge(width = 5), pch = 21, color = "black") + 
  scale_x_discrete(breaks = unique(field_traits_R_sum$Time),
                   #labels = c("Starting", "2021", "2022", "2023", "2024"),
                   labels = c("5" = "5", "17" = "17", "29" = "29", "40" = "40"),  
                   expand = expansion(mult = c(0.1, 0.1))) +
  #scale_y_continuous(limits = c(-0, 12),
  #                   expand = expansion(mult = 0),
  #                   breaks = seq(-0, 12, by = 2)) +
  theme_classic() +
  scale_color_manual(values = c("SV" = "#92BEA7", "PL" = "#A13E38", "MD" = "#D7B91B", "PP" = "#6699CC",
                                "PF" = "#7B7776", "MC" = "#3A55A4", "EP" = "#557278", "EC" = "#653028",
                                "SOR" = "#C15924", "EI" = "#637A34", "UL" = "#BFC525", "LH" = "#AC871F")) + 
  #scale_shape_manual(values = c(1:10)) + 
  #facet_wrap(~ SR2, labeller = labeller(SR2 = c("1" = "No. of species = 1", 
  #                                              "2" = "No. of species = 2", 
  #                                              "4" = "No. of species = 4"))) +
  ggrepel::geom_text_repel(data = subset(field_traits_R_sum, Time == 40),
                           aes(label = species, color = species),direction = "y", hjust = 0, nudge_x = 0.5,segment.size = 0.2, 
                           #segment.color = "inherit",  # 连接线颜色
                           box.padding = 0.3,  point.padding = 0.3,min.segment.length = 0,size = 3.5,show.legend = FALSE) +
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
        strip.background = element_rect(color = "black", size = 0.5, fill = "white"),
        strip.text = element_text(colour='black', size=11), 
        #plot.margin = margin(1,1,1,2,unit="cm"),
        legend.position = "none") + 
  labs(x = "Sampling time (months)", tag = "A",
       y = expression("Relative abundance of species\nin recipient native communities (%)")) -> p2; p2

library(patchwork)
p2|p1




