# BILE QUANT SERVER

server <- function(input, output, session) {
  
  rvalues <- reactiveValues(df_quant_compounds=NULL, df_input=NULL, df_itsd=NULL, df_itsd_stats=NULL, df_cc=NULL, 
                            df_linear_models=NULL, df_plasma_qc=NULL, df_method_blanks=NULL, df_calibrated=NULL, 
                            df_plasma_qc_targets=NULL, 
                            plot_cc=NULL, plot_plasma_qc=NULL, plot_bar1=NULL, plot_bar1=NULL,
                            quant_compounds=NULL, keeprows=NULL, panel=NULL)
  
  ranges <- reactiveValues(x = NULL, y = NULL)
  
  shinyjs::hide("Button_itsd_stats")
  shinyjs::hide("Button_itsd_plot")
  
  
  # 1. INPUT DATA ######################################################################################################

  ## 1.1 Refresh CSV list when hit button ==============================================================================
  
  observeEvent(input$Button_refresh_csv, ignoreInit = T, ignoreNULL = T, {
    updateSelectInput(session, 'filename', choices = list.files(wddir, pattern = "csv$"))
  })

  ## 1.2 Get input data from csv =======================================================================================

  observeEvent(input$Button_upload_csv, ignoreInit = T, ignoreNULL = T, {
    
    shinyjs::show("Button_itsd_stats")
    shinyjs::show("Button_itsd_plot")
    
    filename <- file.path(wddir,input$filename)
    
    # Identify panel and read equivalent quant compounds csv
    rvalues$panel <- case_when(
      grepl("Bile", input$filename, ignore.case = T) ~ "BileAcids",
      grepl("PFBBr", input$filename, ignore.case = T) ~ "PFBBr",
      grepl("Tryptophan|Indole", input$filename, ignore.case = T) ~ "Tryptophan",
      grepl("TMS", input$filename, ignore.case = T) ~ "TMS",
      TRUE ~ "NOT IDENTIFIED... CHECK FILENAME")
    
    # Read and clean input data
    if (rvalues$panel == "BileAcids") {
      rvalues$df_input <- Function_readin_csv_1(filename, zero_threshold) 
    } else {
      rvalues$df_input <- Function_readin_csv_2(filename, zero_threshold) 
    }
  
    # Check data availability for dil/conc for each compound
    dil_conc_found <- rvalues$df_input %>% 
      select(compound_name, conc) %>% 
      group_by(compound_name) %>% 
      mutate(found_dil = ifelse(any(grepl("diluted", conc)), "Yes", "No"),
             found_conc = ifelse(any(grepl("concentrated", conc)), "Yes", "No")) %>% 
      distinct() %>% 
      select(-conc)
    
    rvalues$df_quant_compounds <- get(paste0("df_quant_", rvalues$panel)) %>% 
      mutate(compound_name = tolower(compound_name)) %>% 
      mutate(found = ifelse(compound_name %in% rvalues$df_input$compound_name, "Yes","No"), 
             .after = conc) %>% 
      mutate(conc = factor(conc, levels = sort(unique(conc)))) %>% 
      left_join(dil_conc_found, by="compound_name") %>% 
      relocate(found_dil, .after = found) %>% 
      relocate(found_conc, .after = found) %>% 
      distinct()
    
    
    # Find number of cc points for each quant compound
    num_cc <- rvalues$df_input %>% 
      select(compound_name, cc) %>% 
      filter(compound_name %in% rvalues$df_quant_compounds$compound_name) %>% 
      group_by(compound_name) %>% 
      summarize(num_cc = n_distinct(cc, na.rm = T))

    rvalues$df_quant_compounds <- rvalues$df_quant_compounds %>% 
      left_join(num_cc, by="compound_name") %>% 
      relocate(num_cc, .after = found_dil)
    
    
    # Quant compounds
    rvalues$quant_compounds <- rvalues$df_quant_compounds %>% filter(found=="Yes") %>% pull(compound_name)
    
    # ITSD dataframe
    rvalues$df_itsd <- rvalues$df_input %>% 
      filter(itsd == "ITSD") %>% 
      filter(!grepl("MB|Pooled|Plasma|Standard",sampleid, ignore.case = T)) %>% 
      mutate(peakarea = ifelse(peakarea <= zero_threshold, 0, peakarea),
             num = str_extract(sampleid, "[0-9][0-9][0-9]"),
             num = as.numeric(num),
             cc_shape = ifelse(grepl("CC[0-9]+", sampleid), "Calibration Curve Sample", "Standard Sample"))

    
    # ITSD stats dataframe
    rvalues$df_itsd_stats <- rvalues$df_itsd %>% 
      filter(!grepl("MB|Pooled|Plasma|Standard|CC",sampleid, ignore.case = T)) %>% 
      group_by(compound_name) %>%
      summarise(stdev = sd(peakarea),
                average = mean(peakarea),
                middle = median(peakarea),
                cv = stdev / average,
                cv_med = stdev / median(peakarea)) %>% 
      ungroup()
    
    # Outputs
    output$Textout_panel <- renderText({paste("Panel: ", rvalues$panel)}) # Display identified panel
    output$Textout_filename <- renderText({paste("Uploaded file: ", input$filename)}) # Display filename
    
  })
  
  # Quant compounds settings table (rHandsontableOutput)
  output$Table_calibration_settings <- renderRHandsontable({
    if (!is.null(rvalues$df_quant_compounds))
      rhandsontable(rvalues$df_quant_compounds, stretchH = "all") %>%
      hot_col("conc", allowInvalid = T) %>% 
      hot_col("compound_name", readOnly = T) %>%
      hot_col("found", readOnly = T, halign = "htCenter") %>% 
      hot_col("found_conc", readOnly = T, halign = "htCenter") %>% 
      hot_col("found_dil", readOnly = T, halign = "htCenter") %>% 
      hot_col("num_cc", readOnly = T, halign = "htCenter") %>% 
      hot_col("plasma_qc_target", halign = "htCenter")
  })
  
  
  # ITSD stats (bsModal)
  output$Table_ITSD_stats <- DT::renderDataTable({
    rvalues$df_itsd_stats %>% 
      datatable(options = list(columnDefs = list(list(className='dt-center', targets="_all"))))
  })
  
  # ITSD plot (bsModal)
  output$Plot_ITSD_stats <- renderPlot({
    Function_plot_itsd(rvalues$df_itsd, rvalues$df_itsd_stats)
  })
  
  # 2. CALIBRATE #######################################################################################################

  ## 2.1 Button - Calibration ==========================================================================================
  
  observeEvent(input$Button_calibration, ignoreInit = T, ignoreNULL = T, {
    
    # Update df_cc based on inputs from Table_calibration_settings (rHandsontableOutput)
    rvalues$df_quant_compounds <- hot_to_r(input$Table_calibration_settings)
    saveRDS(rvalues$df_cc, "df_cc.rds")
    df_cc <- readRDS( "df_cc.rds")

    # Normalized dataframe ---------------------------------------------------------------------------------------------
    
    if (rvalues$panel == "BileAcids") { # For bile acids panel
      rvalues$df_normalized <- rvalues$df_input %>% 
        filter(compound_name %in% rvalues$quant_compounds) %>%
        inner_join(rvalues$df_quant_compounds[, c("compound_name", "conc")], by = c("compound_name", "conc")) %>%
        mutate(compound_name=factor(compound_name,levels = rvalues$quant_compounds)) %>%
        replace_na(list(itsd="peak")) %>%
        reshape2::dcast(sampleid+compound_name+conc+letter ~ itsd, value.var="peakarea",
                        fun.aggregate = mean) %>%
        group_by(sampleid, letter) %>%
        mutate(ITSD = zoo::na.locf(ITSD),
               norm_peak = ifelse(ITSD==0, 0, peak / ITSD),
               curveLab = str_extract(sampleid,pattern="CC[1-9][0-9]+|CC[1-9]+"))
    
    }else {
      rvalues$df_normalized <- rvalues$df_input %>% 
        filter(compound_name %in% rvalues$quant_compounds) %>%
        inner_join(rvalues$df_quant_compounds[, c("compound_name", "conc")], by = c("compound_name", "conc")) %>%
        mutate(compound_name=factor(compound_name,levels = rvalues$quant_compounds)) %>%
        replace_na(list(itsd="peak")) %>%
        reshape2::dcast(sampleid+compound_name+conc ~ itsd, value.var="peakarea") %>% 
        mutate(norm_peak = ifelse(ITSD==0, 0, peak / ITSD),
               curveLab = str_extract(sampleid,pattern="CC[1-9][0-9]+|CC[1-9]+"))
      }

    
    
    # Quant dataframe --------------------------------------------------------------------------------------------------
    
    rvalues$df_cc <- rvalues$df_quant_compounds %>% 
      filter(found=="Yes") %>%
      select(-plasma_qc_target, -conc, -found, -found_conc, -found_dil, -num_cc) %>% 
      pivot_longer(!compound_name, names_to = "curveLab", values_to = "conc_val") %>% 
      right_join( rvalues$df_normalized[grepl("CC", rvalues$df_normalized$curveLab), ] , by=c("compound_name","curveLab"))
    
    rvalues$keeprows <- rep(TRUE, nrow(rvalues$df_cc))
    
    # Display CC data table --------------------------------------------------------------------------------------------
    
    output$Table_cc_data <- DT::renderDataTable({
      keep <- rvalues$df_cc[ rvalues$keeprows, , drop = FALSE]
      keep %>% 
        datatable(options = list(columnDefs = list(list(className='dt-center', targets=4:length(rvalues$df_cc)))), 
                  class = 'cell-border stripe') %>% 
        formatRound(columns=c('ITSD', 'peak', 'norm_peak'), digits=2)
    })
    
    
    # Display Linear models table --------------------------------------------------------------------------------------
    
    output$Table_linear_models <- DT::renderDataTable({
      
      keep <- rvalues$df_cc[ rvalues$keeprows, , drop = FALSE]
      
      rvalues$df_linear_models <- keep %>%
        filter(!is.na(norm_peak)) %>% 
        group_by(compound_name) %>%
        summarize(r = cor(norm_peak,conc_val),
                  model_list <- broom::tidy(lm(norm_peak ~ conc_val))) %>%
        reshape2::dcast(compound_name+r ~ term,value.var="estimate") %>%
        dplyr::rename(R=r,
                      Intercept=`(Intercept)`,
                      Slope=conc_val) %>% 
        mutate_if(is.numeric, round, 3)  %>% 
        filter(!is.na(R))

      
      rvalues$df_linear_models %>% 
        dplyr::rename(Compound=compound_name) %>% 
        datatable(options = list(columnDefs = list(list(className='dt-center', targets="_all"))), 
                  class = 'cell-border stripe')
      })
  

    # Plot calibration curves ------------------------------------------------------------------------------------------
    
    output$Plot_calibration_curves <- renderPlot({
      
      keep    <- rvalues$df_cc[ rvalues$keeprows, , drop = FALSE]
      exclude <- rvalues$df_cc[!rvalues$keeprows, , drop = FALSE]

      rvalues$plot_cc <- ggplot(keep, aes(x=conc_val, y=norm_peak)) +
        geom_point(size=3) +
        geom_smooth(method = "lm") +
        geom_point(data = exclude, shape = 21, fill = NA, color = "black", alpha = 0.25) +
        theme(strip.text=element_text(size=15),
              axis.text=element_text(size=13)) +
        theme_bw() +
        facet_wrap("compound_name",scales="free") +
        xlab("Concentration") + ylab("Normalized peak area") +
        coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)
      
      rvalues$plot_cc
      })
    
  })
  
  ## 2.2 Update points to keep from graph click ========================================================================
  observeEvent(input$Plot_calibration_curves_click, {
    res <- nearPoints(rvalues$df_cc, input$Plot_calibration_curves_click, allRows = TRUE)
    rvalues$keeprows <- xor(rvalues$keeprows, res$selected_)
    })

  observeEvent(input$Plot_calibration_curves_dblclick, {
    brush <- input$Plot_calibration_curves_brush
    cat(file=stderr(), paste0(brush))
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges$x <- NULL
      ranges$y <- NULL
    }
  })
  
  ## 2.3 Button - Save calibration curve metrics csv ===================================================================
  
  observeEvent(input$Button_download_cc_metrics, {
    write.csv(rvalues$df_linear_models, paste0("/Volumes/chaubard-lab/shiny_workspace/calibration_metrics/",
                                 gsub(".csv","",input$filename),"_CC_Metrics",".csv"))
    })
  
  
  # 3. CALIBRATED RESULTS ##############################################################################################
  
  observeEvent(input$Button_continue, ignoreInit = T, ignoreNULL = T, {

    # Calibrated table
    rvalues$df_calibrated <- rvalues$df_normalized %>%
      filter(!grepl("^CC[0-9]",sampleid)) %>%
      left_join(rvalues$df_linear_models, by="compound_name") %>%
      mutate(quant_val =  (norm_peak - (Intercept))/Slope*as.numeric(input$Textin_x_factor)) %>%
      arrange(compound_name) %>%
      mutate(quant_val = ifelse(quant_val < 0,0,quant_val),
             quant_val = round(quant_val,2)) %>% 
      filter(!is.na(norm_peak)) %>% 
      filter(!is.na(R))
    
    # Plasma QC dataframe
    rvalues$df_plasma_qc <- rvalues$df_calibrated %>%
      filter(grepl("plasma", sampleid, ignore.case = T)) %>% 
      arrange(compound_name) %>%
      mutate(quant_val = ifelse(quant_val < 0,0,quant_val),
             quant_val = round(quant_val,2))

    
    # Plasma QC range dataframe (for barplot errorbars)
    df_plasma_qc_range <- rvalues$df_plasma_qc %>%
      group_by(compound_name) %>% 
      mutate(max_val = max(quant_val, na.rm = T),
             min_val = min(quant_val, na.rm = T)) %>% 
      select(compound_name, max_val, min_val) %>% 
      distinct()
    
    
    # Plasma QC plot
    output$Plot_plasma_qc <- renderPlot({
      
      rvalues$plot_plasma_qc <- ggplot(rvalues$df_plasma_qc, aes(x=compound_name, y=quant_val)) +
        geom_jitter(position=position_jitter(0.05)) +
        geom_point(data=rvalues$df_quant_compounds, aes(x=compound_name, y=as.numeric(plasma_qc_target)), color="red") +
        theme_bw() +
        xlab("Compound") + ylab(paste0("Concentration (", input$Select_conc_unit,")"))
      
      rvalues$plot_plasma_qc
    })
    
    
    
    # Method blanks dataframe
    rvalues$df_method_blanks <- rvalues$df_calibrated %>%
      filter(grepl("MB",sampleid))
    
    # Method blanks plot
    output$Plot_method_blanks <- renderPlot({
      ggplot(rvalues$df_method_blanks, aes(x=compound_name, y=quant_val)) +
        geom_jitter(position=position_jitter(0.05)) +
        theme_bw() +
        xlab("Compound") + ylab(paste0("Concentration (", input$Select_conc_unit,")"))
    })
    

    
    # Subtract avg method blank values
    df_method_blanks_avg <- rvalues$df_method_blanks %>% 
      group_by(compound_name) %>% 
      summarize(mean_mb = mean(quant_val))
      
    rvalues$df_calibrated <- rvalues$df_calibrated %>% 
      left_join(df_method_blanks_avg, by="compound_name") %>% 
      mutate(quant_val = quant_val - mean_mb)
    
    # Display calibrated table
    output$Table_calibrated_data <- DT::renderDataTable({
      rvalues$df_calibrated %>% 
        datatable(options = list(columnDefs = list(list(className='dt-center', targets="_all"))), 
                  class = 'cell-border stripe')
    })
    
    # Bar plot 1 (separated by compounds)
    
    df_bar <- rvalues$df_calibrated %>%
      ungroup() %>% 
      select(sampleid, compound_name, quant_val) %>%
      filter(!grepl("MB",sampleid, ignore.case = T),
             !grepl("Pooled",sampleid, ignore.case = T),
             !grepl("Plasma",sampleid, ignore.case = T),
             !grepl("Standard",sampleid, ignore.case = T),
             !grepl("CC[0-9]+", sampleid))
    
    output$Plot_bar1 <- renderPlot({
      rvalues$plot_bar1 <- df_bar %>% 
        left_join(df_plasma_qc_range) %>% 
        ggplot(aes(x = sampleid, y = quant_val, fill = compound_name)) +
        geom_bar(stat="identity") +
        geom_errorbar(aes(ymin = quant_val-(min_val+max_val)/2, ymax = quant_val + (min_val+max_val)/2), width = 0.2) +
        theme_bw() +
        theme(plot.margin=unit(c(5.5, 15, 5.5, 10),"points"),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              legend.position = "top",
              legend.text = element_text(size = 14),
              legend.title = element_text(size = 16),
              strip.text=element_text(size=18),
              axis.text.y =element_text(size=11),
              axis.title = element_text(size = 15),
              axis.text.x =element_text(vjust = 0.5,
                                        hjust = 1,
                                        angle = 90,
                                        size = ifelse(
                                          (((nrow(df_bar)))+(nrow(df_bar)))*1000 / (nrow(df_bar) * nrow(df_bar)) >= 11,
                                          11, (((nrow(df_bar)))+(nrow(df_bar)))*1000 / (nrow(df_bar) * nrow(df_bar))))) +
        facet_wrap(~compound_name, scales = "free_y") +
        scale_fill_manual(values = c(pal_ucscgb("default", alpha = 0.7)(7), "bisque4"))+
        xlab("\nSampleID") +
        ylab(paste0("Concentration (", input$Select_conc_unit,")")) +
        guides(fill = guide_legend(title="Compound    "))
      
      rvalues$plot_bar1
    })
    
    output$Plot_bar2 <- renderPlot({
      rvalues$plot_bar2 <- df_bar %>% 
        ggplot(aes(x = sampleid, y = quant_val, fill = compound_name)) +
        geom_bar(stat="identity", position="stack") +
        theme_bw() +
        theme(plot.margin=unit(c(5.5, 15, 5.5, 10),"points"),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              legend.position = "top",
              legend.text = element_text(size = 14),
              legend.title = element_text(size = 16),
              strip.text=element_text(size=18),
              axis.text.y =element_text(size=11),
              axis.title = element_text(size = 15),
              axis.text.x =element_text(vjust = 0.5,
                                        hjust = 1,
                                        angle = 90,
                                        size = ifelse(
                                          (((nrow(df_bar)))+(nrow(df_bar)))*1000 / (nrow(df_bar) * nrow(df_bar)) >= 15,
                                          15, (((nrow(df_bar)))+(nrow(df_bar)))*1000 / (nrow(df_bar) * nrow(df_bar))))) +
        scale_fill_manual(values = c(pal_ucscgb("default", alpha = 0.7)(7), "bisque4"))+
        xlab("\nSampleID") +
        ylab(paste0("Concentration (", input$Select_conc_unit,")")) +
        guides(fill = guide_legend(title="Compound    "))
      
      rvalues$plot_bar2
    })

    
  })
  
  # 4. SAVE RESULTS ####################################################################################################
  
  ## 4.1 Normalized csv ================================================================================================
  observeEvent(input$Button_download_normalized_csv, ignoreInit = T, ignoreNULL = T, {
    
    rvalues$df_normalized %>% 
      select(sampleid, compound_name, norm_peak) %>% 
      mutate(sampleid = ifelse(grepl("MB|Pooled|Plasma|CC|Standard", sampleid, ignore.case = T),
                               sampleid,
                               gsub("^[0-9]{3}_", "", sampleid))) %>% 
      pivot_wider(names_from = compound_name, values_from = norm_peak, values_fill = NA) %>% 
      write.csv(paste0("~/Downloads/normalized_results_",
                       gsub("\\.csv","",input$filename),"_",
                       gsub("\\-","",Sys.Date()),".csv"), row.names=F,quote=F)
    
    shinyalert(title = "File saved to ~/Downloads", type = "success")
  })
  
  
  ## 4.2 Normalized csv (No QCs) =======================================================================================
  observeEvent(input$Button_download_normalized_csv_no_qc, ignoreInit = T, ignoreNULL = T, {
    
    rvalues$df_normalized %>% 
      select(sampleid, compound_name, norm_peak) %>% 
      filter(!grepl("MB|Pooled|Plasma|CC|Standard",sampleid, ignore.case = T)) %>% 
      mutate(sampleid = gsub("^[0-9]{3}_", "", sampleid)) %>% 
      pivot_wider(names_from = compound_name, values_from = norm_peak, values_fill = NA) %>% 
      write.csv(paste0("~/Downloads/removed_qcs_normalized_results_",
                       gsub("\\.csv","",input$filename),"_",
                       gsub("\\-","",Sys.Date()),".csv"), row.names=F,quote=F)
    
    shinyalert(title = "File saved to ~/Downloads", type = "success")
    
  })
  
  ## 4.3 Quant csv =====================================================================================================
  
  observeEvent(input$Button_download_quant_csv, ignoreInit = T, ignoreNULL = T, {
    
    rvalues$df_calibrated %>% 
      select(sampleid, compound_name, norm_peak) %>% 
      mutate(sampleid = ifelse(grepl("MB|Pooled|Plasma|CC|Standard", sampleid, ignore.case = T),
                               sampleid,
                               gsub("^[0-9]{3}_", "", sampleid))) %>% 
      pivot_wider(names_from = compound_name, values_from = norm_peak, values_fill = NA) %>% 
      write.csv(paste0("~/Downloads/quant_results_",
                       gsub("\\.csv","",input$filename),"_",
                       gsub("\\-","",Sys.Date()),".csv"), row.names=F,quote=F)
    
    shinyalert(title = "File saved to ~/Downloads", type = "success")
  })
  
  
  ## 4.4 Quant csv (No QCs) ============================================================================================
  
  observeEvent(input$Button_download_quant_csv_no_qc, ignoreInit = T, ignoreNULL = T, {
    
    rvalues$df_calibrated %>% 
      select(sampleid, compound_name, norm_peak) %>% 
      filter(!grepl("MB|Pooled|Plasma|CC|Standard",sampleid, ignore.case = T)) %>% 
      mutate(sampleid = gsub("^[0-9]{3}_", "", sampleid)) %>% 
      pivot_wider(names_from = compound_name, values_from = norm_peak, values_fill = NA) %>% 
      write.csv(paste0("~/Downloads/removed_qcs_quant_results_",
                       gsub("\\.csv","",input$filename),"_",
                       gsub("\\-","",Sys.Date()),".csv"), row.names=F,quote=F)
    
    shinyalert(title = "File saved to ~/Downloads", type = "success")
    
  })

  ## 4.5 Barplots ======================================================================================================

  output$Button_download_barplots <- downloadHandler(

    filename = function(){
      paste0(rvalues$panel, "_Barplots_",unique(rvalues$df_input$batch),"_",Sys.Date(),".pdf")
    },

    content = function(file) {
      pdf(file, height = 11, width = 8.5)
      print(rvalues$plot_bar1)
      print(rvalues$plot_bar2)
      dev.off()
    })

  ## 4.6 QC report =====================================================================================================

  output$Button_download_qc_report <- downloadHandler(

    filename = function(){
      paste0(rvalues$panel, "_QC_Report_",unique(rvalues$df_input$batch),"_",Sys.Date(),".pdf")
    },

    content = function(file) {

      pdf(file, height = 11, width = 8.5)
      
      print(Function_plot_itsd(rvalues$df_itsd, rvalues$df_itsd_stats) +
              ggtitle(paste("Bile Acid Quantitative QC Report\n", unique(rvalues$df_input$batch))) +
              theme(plot.title = element_text(color = "black", hjust = 0.5, size = 20, face = "bold")))
      print(rvalues$plot_plasma_qc)
      
      ba_tt <- gridExtra::ttheme_default(
        colhead=list(fg_params=list(col="black", fontface=2L)),
        padding = unit(c(0.5,0.75), "cm"))
      print( gridExtra::grid.arrange(gridExtra::tableGrob(rvalues$df_itsd_stats, rows = NULL, theme = ba_tt)) )
      
      print(rvalues$plot_cc)

      dev.off()
    })
   
}