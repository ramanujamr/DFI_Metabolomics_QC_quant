library(shiny)
library(shinythemes)
library(shinyjs)
library(shinyBS)
library(shinyscreenshot)
library(shinyalert)
library(rhandsontable)

library(tidyverse)
library(reshape2)
library(DT)
library(data.table)

library(broom)
library(ggsci)
library(gridExtra)
library(ggpmisc)




# 0. Initializations ###################################################################################################

# quant_compounds <- c("cholic acid", "deoxycholic acid", "lithocholic acid", "glycocholic acid", "taurocholic acid", 
#                      "isodeoxycholic acid", "alloisolithocholic acid", "3-oxolithocholic acid")
# 
# num_quant_compounds <- length(quant_compounds)
zero_threshold <- 1000

# 1. set up directory and files ########################################################################################

#wddir <- "/Volumes/chaubard-lab/shiny_workspace/csvs/"
wddir <- "/Users/ramanujam/GitHub/dfimmfshiny_test/test_files"

#'* Read quant tables*
df_quant_BileAcids <- read.csv("quant_compounds_bile_acids.csv")
df_quant_PFBBr <- read.csv("quant_compounds_pfbbr.csv")
df_quant_TMS <- read.csv("quant_compounds_tms.csv")
df_quant_Tryptophan <- read.csv("quant_compounds_tryptophan.csv")



# 2. FUNCTIONS #########################################################################################################

## 2.1 Read and clean data from input file =============================================================================


Function_readin_csv_1 <- function(filename, zero_threshold, recursive=F)
{
  df_input <- read.csv(file = filename) %>%
    select(-CAS.ID) %>%
    reshape2::melt(id.vars=c("Compound.Name", "Formula", "Mass", "RT")) %>%
    replace_na(list(value = 0)) %>%
    mutate(itsd = str_extract(Compound.Name,pattern="ITSD"),
           com = gsub("\\_ITSD","",Compound.Name),
           Data.File=variable) %>%
    separate(com,into=c("num","compound_name","letter"),sep="\\_") %>%
    separate(variable,into=c("num2","date_run","batch","sampleid","conc"),sep="\\_\\_") %>%
    mutate(num2 = gsub("[Xx]", "", num2),
           sampleid = paste(num2, sampleid, sep = "_")) %>% 
    select(Data.File, sampleid, date_run, Compound.Name, compound_name,
           batch, letter, itsd, conc, peakarea=value) %>%
    mutate(filename = filename,
           compound_name = gsub("D[0-9]+\\-","",compound_name),
           compound_name = tolower(compound_name),
           conc = ifelse(grepl("dil",conc),"diluted","concentrated"),
           peakarea = as.numeric(peakarea),
           peakarea = ifelse(peakarea <= zero_threshold, 0, peakarea),
           cc = str_extract(sampleid, pattern=".*CC.*"),
           cc = sub(".*_", "", cc))
  return(df_input)
}


Function_readin_csv_2 <- function(filename, zero_threshold, recursive=F){
  
  df_input_raw <- read.csv(file=filename, check.names=F)
  colnames(df_input_raw)[2] <- "garbage"
  colnames(df_input_raw)[3] <- "sampleid"
  colnames(df_input_raw)[4] <- "Data.File"
  colnames(df_input_raw)[5] <- "Type"
  colnames(df_input_raw)[6] <- "Level"
  colnames(df_input_raw)[7] <- "Acq.Date.Time"
  colnames(df_input_raw) <- gsub(" Results","",colnames(df_input_raw))
  df_input_raw = df_input_raw[-1,]
  
  df_input <- df_input_raw %>% 
    select(-Sample, -garbage,-Type,-Level,-Acq.Date.Time, -Data.File) %>% 
    slice(-1) %>% 
    reshape2::melt(id.vars="sampleid") %>% 
    mutate(value = as.numeric(value)) %>% 
    replace_na(list(value = 0)) %>% 
    dplyr::rename(compound_name = variable) %>% 
    mutate(itsd = str_extract(compound_name, pattern="ITSD"),
           compound_name = tolower(gsub("\\_ITSD","",compound_name))) %>% 
    separate(sampleid, into=c("inj_num","date_run","batch","sampleid","conc"),sep="\\_\\_") %>% 
    mutate(sampleid = paste0(inj_num,"_",sampleid)) %>% 
    select(sampleid, date_run, batch, compound_name, itsd, conc, peakarea=value) %>% 
    mutate(conc = ifelse(grepl("dil",conc),"diluted","concentrated"),
           peakarea = ifelse(peakarea <= zero_threshold, 0, peakarea),
           cc = str_extract(sampleid, pattern=".*CC.*"),
           cc = sub(".*_", "", cc))
  
  return(df_input)
}


## 2.1 Plot ITSD =======================================================================================================

Function_plot_itsd <- function(df_itsd, df_itsd_stats)
{
  df_itsd %>%
    filter(!grepl("MB",sampleid, ignore.case = T),
           !grepl("plasma",sampleid, ignore.case = T)) %>% 
    left_join(df_itsd_stats, by = "compound_name") %>%
    mutate(flag = ifelse(peakarea > average + (1.5 * stdev) |
                           peakarea < average - (1.5 * stdev), paste(num,batch,sampleid, sep = "_"), NA)) %>%
    group_by(compound_name) %>%
    ggplot(., aes(x = num, y = peakarea, group = compound_name, label = flag)) +
    geom_point(aes(color = compound_name, shape = cc_shape), fill = "black", alpha = 0.6, size = 2) +
    geom_line(aes(y = average, color = compound_name)) +
    geom_line(aes(y = average + 1.5*stdev, color = compound_name), linetype = "dashed") +
    geom_line(aes(y = average - 1.5*stdev, color = compound_name), linetype = "dashed") +
    geom_ribbon(aes(ymin = average - stdev, ymax = average + stdev,
                    fill = compound_name), alpha=0.2) +
    ggrepel::geom_label_repel(size = 1.2, max.overlaps = Inf,
                              min.segment.length = 0.1, label.padding = 0.1, na.rm=T) +
    theme_bw() +
    theme(panel.grid.minor= element_blank(),
          panel.grid.major.x = element_blank(),
          legend.text = element_text(color = "black", size = 8),
          legend.title = element_text(color = "black", size = 10),
          legend.position = "top",
          strip.text=element_text(color = "black", size=5),
          axis.text =element_text(color = "black", size=8),
          axis.title = element_text(color = "black", size = 14),
          plot.margin = margin(1,0.5,0,0.6, unit = 'cm')) +
    ggsci::scale_color_ucscgb() +
    ggsci::scale_fill_ucscgb() +
    guides(color = guide_legend(title = "Internal Standard Compound",
                                override.aes = list(size = 2.5), nrow = 2,
                                title.position="top", title.hjust = 0.5,
                                label.position = "right"), fill = "none",
           shape = guide_legend(title = "",
                                override.aes = list(size = 2.5), nrow = 2,
                                title.position="top", title.hjust = 0.5,
                                label.position = "right")) +
    scale_shape_manual(values = c(24,16))+
    scale_y_continuous(label = scales::scientific) +
    scale_x_continuous(breaks = seq(0,150,25)) +
    ylab("Raw Peak Area\n") +
    xlab("\nInjection Number")+
    facet_grid(~compound_name, scales="free_x")
}




# 2. Execution #########################################################################################################

source('ui.R', local=TRUE)
source('server.R', local=TRUE)
shinyApp(ui=ui, server=server)