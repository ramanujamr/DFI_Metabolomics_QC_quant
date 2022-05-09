start5 <- 125
compounds5 <- c("Cholic Acid", "Deoxycholic Acid", "Lithocholic Acid", "Glycocholic Acid", "Taurocholic Acid", 
                "Isodeoxycholic Acid", "Alloisolithocholic Acid", "3-Oxolithocholic Acid")
quant_conc5 <- c("dil","dil","dil","dil","dil","dil","dil","dil")
series5 <- 10 # dilution #

conc_tbl5 <- get_all_conc(start5, compounds=tolower(compounds5),series=series5,fold=3)

maxcc5 <- c("125,125,125,125,125,125,125,125")
mincc5 <- c("0,0,0,0,0,0,0,0")

# Input file -------------------------------------------------
filename="bile_acid_test.csv"
na.value=0
na.replacement=1000



# Internal std graph
# Internal Std stats table
# Calibration curves with plasma qc variation



df_itsd <- df_input %>% 
  filter(itsd == "ITSD") %>% 
  filter(!grepl("MB",sampleid, ignore.case = T),
         !grepl("Pooled",sampleid, ignore.case = T),
         !grepl("BHIQC",sampleid, ignore.case = T),
         !grepl("Plasma",sampleid, ignore.case = T),
         !grepl("Hexanes",sampleid, ignore.case = T),
         !grepl("Standard",sampleid, ignore.case = T),
         !grepl("50%_MeOH",sampleid, ignore.case = T),
         !grepl("50%MeOH",sampleid, ignore.case = T),
         !grepl("EA_BLANK",sampleid, ignore.case = T)) %>%
  mutate(peakarea = ifelse(peakarea <= zero_threshold, 0, peakarea),
         num = str_extract(sampleid, "[0-9][0-9][0-9]"),
         num = as.numeric(num),
         cc_shape = ifelse(grepl("CC[0-9]+", sampleid), "Calibration Curve Sample", "Standard Sample"))

df_itsd_stats <- df_itsd %>% 
  group_by(compound_name) %>%
  summarise(stdev = sd(peakarea),
            average = mean(peakarea),
            middle = median(peakarea),
            cv = stdev / average,
            cv_med = stdev / median(peakarea)) %>% 
  ungroup()
  

  














## Save report ##
output$ba_qc_report_download <- downloadHandler(
  
  filename = function(){
    paste0(rvalues$panel, "_QC_Report_",unique(rvalues$df_input$batch),"_",Sys.Date(),".pdf")
    },
  
  content = function(file) {
    
    pdf(file, height = 11, width = 8.5)
    
    print(Function_plot_itsd(rvalues$df_itsd, rvalues$df_itsd_stats) +
          ggtitle(paste("Bile Acid Quantitative QC Report\n", unique(meta5_1()$batch))) +
          theme(plot.title = element_text(color = "black", hjust = 0.5, size = 20, face = "bold")))
    
    print(rvalues$plot_plasma_qc)
    
    print( gridExtra::grid.arrange(gridExtra::tableGrob(rvalues$df_itsd_stats, rows = NULL, theme = ba_tt)) )
    
    print(rvalues$plot_cc)
    
    dev.off()
  })

