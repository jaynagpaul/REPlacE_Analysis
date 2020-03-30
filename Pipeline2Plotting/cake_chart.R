setwd("C:/Users/User/Documents/Lab/LongBCRs/bioinformatics/eric_faa/")

library(ggplot2)
library(tidyr)
library(dplyr)
pie <- read.table("pie_chart_fractions.txt", header = TRUE,
                  sep = "\t", stringsAsFactors = FALSE)

pie <- pie %>% mutate(method = factor(method, levels = c("tn5", "lam")),
                      type = factor(type, levels = c("wt", "repl", "repl_rev",
                                                     "del", "inv", "plasmid")))

ggplot(pie %>% dplyr::filter(method == "tn5"))+
  geom_hline(yintercept = seq(5,35,5), linetype = "dashed", color = "grey")+
  geom_bar(aes(x = type, fill = method, y = fraction), stat = "identity",
           position = position_dodge2(), color = "black")+
  theme_classic()+
  scale_y_continuous(expand = c(0,0), breaks = seq(0,35,5))+
  coord_flip()+
  scale_fill_manual(values = c("tn5" = "grey50", "lam" = "white"),
                    labels = c("tn5" = "Tn5", "lam" = "LAM-HTGTS"),
                    name = "Method")+
  scale_x_discrete(labels = c("wt" = "WT", "repl" = "Replacement",
                              "repl_rev" = "Reversed\nreplacement",
                              "del" = "Deletion",
                              "inv" = "Inversion",
                              "plasmid" = "Plasmid\nintegration"),
                   limits = c("plasmid", "inv", "del", "repl_rev", "repl", "wt"))+
  ylab("Fraction, %")+
  xlab("Rearrangement type")

ggsave(paste('barplot_',
             format(Sys.time(), "%d%m%y_%H%M"),
             '.svg', sep = ''),
       width = 6, height = 5,
       device = "svg")


# Cake chart --------------------------------------------------------------
pie <- pie %>% mutate(method = factor(method, levels = c("tn5", "lam")),
                      type = factor(type, levels = c("plasmid", "inv", "del",
                                                     "repl_rev", "repl", "wt")))


ggplot(pie %>% filter(method == "lam"))+
  geom_bar(aes(x = method, y = fraction, fill = type), stat = "identity",
           position = position_stack(), color = "black")+
  theme_classic()+
  scale_y_continuous(expand = c(0,0), breaks = seq(0,100,10))+
  coord_polar(theta = "y")+
  theme(axis.line = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.ticks.y = element_blank(), axis.text.y = element_blank())+
  scale_fill_manual(values = c("wt" = rgb(200,200,200,maxColorValue = 255),
                               "repl" = rgb(255,0,0,maxColorValue = 255),
                               "repl_rev" = rgb(200,0,0,maxColorValue = 255),
                               "del" = rgb(100,100,100,maxColorValue = 255),
                               "inv" = rgb(150,150,150,maxColorValue = 255),
                               "plasmid" = rgb(0,0,0,maxColorValue = 255)),
                    labels = c("wt" = "WT", "repl" = "Replacement",
                               "repl_rev" = "Reversed\nreplacement",
                               "del" = "Deletion",
                               "inv" = "Inversion",
                               "plasmid" = "Plasmid\nintegration"),
                    limits = c("wt", "repl", "repl_rev", "del", "inv", "plasmid"),
                    name = "Rearrangement type")


ggsave(paste('piechart_',
             format(Sys.time(), "%d%m%y_%H%M"),
             '.svg', sep = ''),
       width = 6, height = 5,
       device = "svg")



# Detailed  ---------------------------------------------------------------
pie <- read.table("pie_chart_detailed.txt", header = TRUE,
                  sep = "\t", stringsAsFactors = FALSE)

pie <- pie %>% bind_rows(data.frame(type = "plasmid",
                                    count = 0,
                                    detailed.type = "intact"))


pie <- pie %>% mutate(detailed.type = factor(detailed.type, levels = c("ins", "del", "intact")),
                      type = factor(type, levels = c("plasmid","inv", "del", "repl_rev",
                                                     "repl", "wt")),
                      fraction = count/sum(pie$count))



ggplot(pie)+
  geom_hline(yintercept = seq(0,1,0.2), linetype = "dashed", color = "grey")+
  geom_bar(aes(x = type, fill = detailed.type, y = count), stat = "identity",
           position = position_fill(), color = "black")+
  theme_classic()+
  scale_y_continuous(expand = c(0,0), breaks = seq(0,1,0.2), labels = scales::percent)+
  coord_flip()+
  scale_fill_manual(values = c("intact" = "white", "del" = "grey50", "ins" = "grey"),
                    name = "Type")+
  # scale_x_discrete(labels = c("wt" = "WT", "repl" = "Replacement",
  #                             "repl_rev" = "Reversed\nreplacement",
  #                             "del" = "Deletion",
  #                             "inv" = "Inversion",
  #                             "plasmid" = "Plasmid\nintegration"),
  #                  limits = c("plasmid", "inv", "del", "repl_rev", "repl", "wt"))+
  ylab("Fraction, %")+
  xlab("Rearrangement type")

ggsave(paste('barplot_',
             format(Sys.time(), "%d%m%y_%H%M"),
             '.svg', sep = ''),
       width = 6, height = 5,
       device = "svg")


# LAM pie chart -----------------------------------------------------------
pie <- read.table("pie_chart_lam.txt", header = TRUE,
                  sep = "\t", stringsAsFactors = FALSE)

pie <- pie %>% dplyr::filter(method != "tn5")

pie <- pie %>% mutate(type = factor(type, levels = c("concat", "transloc", "target")))

ggplot(pie %>% filter(method == "tn5_coll"  & primer == "mcherry_rev"))+
  geom_bar(aes(x = method, y = fraction, fill = type), stat = "identity",
           position = position_stack(), color = "black")+
  theme_classic()+
  scale_y_continuous(expand = c(0,0), breaks = seq(0,100,10))+
  coord_polar(theta = "y")+
  theme(axis.line = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.ticks.y = element_blank(), axis.text.y = element_blank())+
  scale_fill_manual(values = c("target" = rgb(0,0,255,maxColorValue = 255),
                               "concat" = rgb(200,200,200,maxColorValue = 255),
                               "transloc" = rgb(255,220,220,maxColorValue = 255)),
                    labels = c("target" = "Target site\nintegration",
                               "concat" = "Concatemerization",
                               "transloc" = "Off-target\nintegration"),
                    name = "Integration type")

ggsave(paste('piechart_',
             format(Sys.time(), "%d%m%y_%H%M"),
             '.svg', sep = ''),
       width = 6, height = 5,
       device = "svg")




