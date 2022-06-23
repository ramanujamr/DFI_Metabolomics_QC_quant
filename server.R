# QUANT SERVER

server <- function(input, output, session) {
  
  rvalues <- reactiveValues(df_quant_compounds=NULL, df_input=NULL, df_itsd=NULL, df_itsd_stats=NULL, df_cc=NULL, 
                            df_linear_models=NULL, df_plasma_qc=NULL, df_method_blanks=NULL, df_calibrated=NULL, 
                            df_plasma_qc_targets=NULL, 
                            plot_cc=NULL, plot_plasma_qc=NULL, plot_bar1=NULL, plot_bar1=NULL,
                            quant_compounds=NULL, keeprows=NULL, panel=NULL)
  
  ranges <- reactiveValues(x = NULL, y = NULL)
  
  shinyjs::hide("Button_itsd_stats")
  shinyjs::hide("Button_new_file")
  
  
  # 1. INPUT DATA ######################################################################################################

  ## 1.1 Refresh CSV list when hit button ==============================================================================
  observeEvent(input$Button_refresh_csv, ignoreInit = T, ignoreNULL = T, {
    updateSelectInput(session, 'filename', choices = rev(list.files(wddir, 
                                                                    pattern = ".*bile.*csv|.*PFBBr.*csv|.*Indole.*csv|.*Tryptophan.*csv", 
                                                                    ignore.case = T)))
  })
  
  ## 1.2 Upload new file ===============================================================================================
  observeEvent(input$Button_new_file, ignoreInit = T, ignoreNULL = T, {
    session$reload()
  })
  


  ## 1.2 Button Upload CSV #############################################################################################

  observeEvent(input$Button_upload_csv, ignoreInit = T, ignoreNULL = T, {
    
    shinyjs::hide("filename")
    shinyjs::hide("Button_refresh_csv")
    shinyjs::hide("Button_upload_csv")
    
    shinyjs::show("Button_new_file")
    shinyjs::show("Button_itsd_stats")

    
    filename <- file.path(wddir,input$filename)
    
    # Identify panel and read equivalent quant compounds csv --------------------------------
    rvalues$panel <- case_when(
      grepl("Bile", input$filename, ignore.case = T) ~ "BileAcids",
      grepl("PFBBr", input$filename, ignore.case = T) ~ "PFBBr",
      grepl("Tryptophan|Indole", input$filename, ignore.case = T) ~ "Tryptophan",
      grepl("TMS", input$filename, ignore.case = T) ~ "TMS",
      TRUE ~ "NOT IDENTIFIED... CHECK FILENAME")
    
    # Set zero threshold
    zero_threshold <- ifelse(rvalues$panel=="Tryptophan", 100, 1000)
    
    # Read and clean input data --------------------------------
    if (rvalues$panel == "BileAcids") {
      rvalues$df_input <- Function_readin_csv_1(filename, zero_threshold) 
    } else {
      rvalues$df_input <- Function_readin_csv_2(filename, zero_threshold) 
    }
  
    # Check data availability for dil/conc for each compound --------------------------------
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
    
    
    # Find number of cc points for each quant compound --------------------------------
    num_cc <- rvalues$df_input %>% 
      select(compound_name, cc) %>% 
      filter(compound_name %in% rvalues$df_quant_compounds$compound_name) %>% 
      group_by(compound_name) %>% 
      summarize(num_cc = n_distinct(cc, na.rm = T))

    rvalues$df_quant_compounds <- rvalues$df_quant_compounds %>% 
      left_join(num_cc, by="compound_name") %>% 
      relocate(num_cc, .after = found_dil)
    
    
    # Get quant compounds --------------------------------
    rvalues$quant_compounds <- rvalues$df_quant_compounds %>% filter(found=="Yes") %>% pull(compound_name)
    
    # ITSD dataframe (per sample for normalization) --------------------------------
    rvalues$df_itsd <- rvalues$df_input %>% 
      filter(itsd == "ITSD") %>% 
      filter(!grepl("MB|Pooled|Plasma|Standard|Spiked",sampleid, ignore.case = T)) %>% 
      mutate(peakarea = ifelse(peakarea <= zero_threshold, 0, peakarea),
             num = str_extract(sampleid, "[0-9][0-9][0-9]"),
             num = as.numeric(num),
             cc_shape = ifelse(grepl("CC[0-9]+", sampleid), "Calibration Curve Sample", "Standard Sample"))
    
    if (rvalues$panel == "BileAcids") { # For bile acids panel
      rvalues$df_itsd_stats <- rvalues$df_input %>%
        filter(itsd=="ITSD") %>%
        group_by(sampleid, letter) %>%
        summarize(avg = mean(peakarea),
                  med = median(peakarea)) %>% 
        ungroup()
    } else {  # For non bile acids panel
      
      rvalues$df_itsd_stats <- rvalues$df_input %>%
        filter(itsd=="ITSD") %>%
        group_by(sampleid) %>%
        summarize(avg = mean(peakarea),
                  med = median(peakarea)) %>% 
        ungroup()
    }
    
    
    # ITSD stats per compound (for QC) --------------------------------
    
    rvalues$df_itsd_stats_compounds <- rvalues$df_input %>%
      filter(itsd=="ITSD") %>%
      filter(!grepl("MB|Pooled|Plasma|CC|Standard",sampleid, ignore.case = T)) %>% 
      mutate(peakarea = ifelse(peakarea <= zero_threshold, 0, peakarea)) %>% 
      group_by(batch, compound_name) %>%
      summarise(stdev = sd(peakarea),
                average = mean(peakarea),
                middle = median(peakarea),
                cv = stdev / average,
                cv_med = stdev / median(peakarea)) # Don't turn into % here since it will be applied in the y-axis scale
    
    rvalues$df_itsd_compounds <- rvalues$df_input %>% 
      filter(itsd == "ITSD") %>% 
      filter(!grepl("MB|Pooled|Plasma|Standard|Spiked",sampleid, ignore.case = T)) %>% 
      mutate(peakarea = ifelse(peakarea <= zero_threshold, 0, peakarea),
             num = str_extract(sampleid, "[0-9][0-9][0-9]"),
             num = as.numeric(num),
             cc_shape = ifelse(grepl("CC[0-9]+", sampleid), "CC Sample", "ITSD")) %>% 
      left_join(rvalues$df_itsd_stats_compounds) %>% 
      mutate(num = as.numeric(num),
             flag = ifelse(peakarea > average + (1.5 * stdev) | peakarea < average - (1.5 * stdev), 
                           paste(num,batch,sampleid, sep = "_"), NA))

    
    
    # Outputs -------------------
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
    rvalues$df_itsd_stats_compounds %>% 
      dplyr::rename(Batch = batch,
                    `Internal Standard` = compound_name,
                    StDev = stdev,
                    Mean = average,
                    Median = middle,
                    CV = cv,
                    `CV Median` = cv_med) %>%
      mutate(StDev = round(StDev, digits = 0),
             Mean = round(Mean, digits = 0),
             Median = round(Median, digits = 0),
             `CV (%)` = round(CV * 100, digits = 1),
             `CV Median (%)` = round(`CV Median` * 100, digits = 1)) %>%
      select(Batch,`Internal Standard`,StDev,Mean,Median,`CV (%)`,`CV Median (%)`) %>% 
      datatable(options = list(columnDefs = list(list(className='dt-center', targets="_all"))))
      
  })
  
  
  # 2. CALIBRATE #######################################################################################################

  ## 2.1 Button - Calibration ==========================================================================================
  
  observeEvent(input$Button_calibration, ignoreInit = T, ignoreNULL = T, {
    
    # Update df_cc based on inputs from Table_calibration_settings (rHandsontableOutput)
    rvalues$df_quant_compounds <- hot_to_r(input$Table_calibration_settings)
    
    #saveRDS(rvalues$df_quant_compounds, "df_quant_compounds.rds")
    #df_quant_compounds <- readRDS("df_quant_compounds.rds")

    # Normalized dataframe ---------------------------------------------------------------------------------------------
    
    if (rvalues$panel == "BileAcids") { # For bile acids panel
      
      rvalues$df_normalized <- rvalues$df_input %>% 
        filter(compound_name %in% rvalues$quant_compounds) %>%
        filter(is.na(itsd)) %>% 
        inner_join(rvalues$df_quant_compounds[, c("compound_name", "conc")], by = c("compound_name", "conc")) %>%
        mutate(compound_name=factor(compound_name,levels = rvalues$quant_compounds)) %>%
        left_join(rvalues$df_itsd_stats[, c("sampleid", "letter", "avg")], by= c("sampleid", "letter")) %>% 
        rename(ITSD=avg, peak=peakarea) %>% 
        mutate(norm_peak = ifelse(ITSD==0, 0, peak / ITSD),
               curveLab = str_extract(sampleid,pattern="CC[1-9][0-9]+|CC[1-9]+")) %>% 
        ungroup() %>% 
        mutate(norm_peak = ifelse(is.na(norm_peak), 0, norm_peak))
    
    }else {
      rvalues$df_normalized <- rvalues$df_input %>% 
        filter(compound_name %in% rvalues$quant_compounds) %>%
        inner_join(rvalues$df_quant_compounds[, c("compound_name", "conc")], by = c("compound_name", "conc")) %>%
        mutate(compound_name=factor(compound_name,levels = rvalues$quant_compounds)) %>%
        replace_na(list(itsd="peak")) %>%
        reshape2::dcast(sampleid+compound_name+conc ~ itsd, value.var="peakarea", fun.aggregate = mean) %>% 
        mutate(norm_peak = ifelse(ITSD==0, 0, peak / ITSD),
               curveLab = str_extract(sampleid,pattern="CC[1-9][0-9]+|CC[1-9]+")) %>% 
        ungroup() %>% 
        mutate(norm_peak = ifelse(is.na(norm_peak), 0, norm_peak))
      }

    
    
    # Quant dataframe --------------------------------------------------------------------------------------------------
    
    rvalues$df_cc <- rvalues$df_quant_compounds %>% 
      filter(found=="Yes") %>%
      select(-plasma_qc_target, -conc, -found, -found_conc, -found_dil, -num_cc) %>% 
      pivot_longer(!compound_name, names_to = "curveLab", values_to = "conc_val") %>% 
      right_join( rvalues$df_normalized[grepl("CC", rvalues$df_normalized$curveLab), ] , by=c("compound_name","curveLab"))
    
    rvalues$df_cc_original <- rvalues$df_cc
    rvalues$keeprows <- rep(TRUE, nrow(rvalues$df_cc))
    
    # Display CC data table --------------------------------------------------------------------------------------------
    
    output$Table_cc_data <- DT::renderDataTable({
      
      rvalues$df_cc %>% 
        select(compound_name, curveLab, conc_val, peak, ITSD, norm_peak) %>% 
        mutate(status = rvalues$keeprows) %>% 
        datatable(options = list(columnDefs = list(list(className='dt-center', targets=2:7))), 
                  class = 'cell-border stripe') %>% 
        formatRound(columns=c('ITSD', 'peak', 'norm_peak'), digits=2)
    })
    
    
    # Update and display Linear models table ---------------------------------------------------------------------------
    
    output$Table_linear_models <- DT::renderDataTable({
      
      keep <- rvalues$df_cc[ rvalues$keeprows, , drop = FALSE]
      
      if( rvalues$panel == "Tryptophan" & input$Checkbox_intercept_zero==T ) { # For tryptophan panel only
        
        rvalues$df_linear_models <- keep %>%
          filter(!is.na(norm_peak)) %>% 
          group_by(compound_name) %>%
          summarize(r = cor(norm_peak,conc_val),
                    model_list <- broom::tidy(lm(norm_peak ~ conc_val + 0))) %>%
          reshape2::dcast(compound_name+r ~ term,value.var="estimate") %>%
          dplyr::rename(R=r,
                        Slope=conc_val) %>% 
          mutate_if(is.numeric, round, 3)  %>% 
          filter(!is.na(R)) %>% 
          mutate(Intercept = 0)
        
        } else {
          
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
        
        }
      
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
  
  ## Update CC range (max and min) =====================================================================================
  
  observeEvent(input$Button_cc_range, {

    rvalues$df_cc <- rvalues$df_cc_original %>% 
      filter(conc_val <= input$Numericinput_max_cc,
             conc_val >= input$Numericinput_min_cc)
    
    rvalues$keeprows <- rep(TRUE, nrow(rvalues$df_cc))
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
  
  observeEvent(input$Button_save_cc_metrics, {
    rvalues$df_linear_models %>% 
      left_join(rvalues$df_input %>% dplyr::select(compound_name, batch), by="compound_name") %>% 
      left_join(rvalues$df_quant_compounds %>% dplyr::select(compound_name, conc), by="compound_name") %>% 
      select(Batch = batch, `Internal Standard`=compound_name, `Concentration`=conc, r=R, Intercept, Slope) %>% 
      distinct() %>% 
      write.csv(paste0("/Volumes/chaubard-lab/shiny_workspace/calibration_metrics/",
                       gsub(".csv","",input$filename),"_CC_Metrics",".csv"))
    
    shinyalert(title = "File saved to /Volumes/chaubard-lab/shiny_workspace/calibration_metrics/", type = "success")
    
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
             quant_val = ifelse(is.na(quant_val), 0, quant_val),
             quant_val = round(quant_val,2)) %>% 
      filter(!is.na(norm_peak)) %>% 
      filter(!is.na(R)) %>% 
      ungroup()
    
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
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
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
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        xlab("Compound") + ylab(paste0("Concentration (", input$Select_conc_unit,")"))
    })
    

    
    # Subtract avg method blank values
    df_method_blanks_avg <- rvalues$df_method_blanks %>% 
      group_by(compound_name) %>% 
      summarize(mean_mb = mean(quant_val))
    
    rvalues$df_calibrated <- rvalues$df_calibrated %>% 
      left_join(df_method_blanks_avg, by="compound_name") %>% 
      rowwise() %>%
      mutate(quant_val = ifelse(input$Checkbox_subtract_MB==T, quant_val - mean_mb, quant_val))
    
    # Display calibrated table
    output$Table_calibrated_data <- DT::renderDataTable({
      if(rvalues$panel=="BileAcids"){
        rvalues$df_calibrated %>% 
          select(-Data.File, -date_run, -Compound.Name, - filename, -cc) %>% 
          datatable(options = list(columnDefs = list(list(className='dt-center', targets="_all"))), 
                    class = 'cell-border stripe')
      }else{
        rvalues$df_calibrated %>% 
          datatable(options = list(columnDefs = list(list(className='dt-center', targets="_all"))), 
                    class = 'cell-border stripe')
      }
    })
    
    # Bar plot 1 (separated by compounds)
    
    rvalues$df_bar <- rvalues$df_calibrated %>%
      ungroup() %>% 
      select(sampleid, compound_name, quant_val) %>%
      filter(!grepl("MB|Pooled|Plasma|Standard|CC[0-9]+",sampleid, ignore.case = T))
    
    output$Plot_bar1 <- renderPlot({
      rvalues$plot_bar1 <- rvalues$df_bar %>% 
        left_join(df_plasma_qc_range) %>% 
        ggplot(aes(x = sampleid, y = quant_val, fill = compound_name)) +
        geom_bar(stat="identity") +
        #geom_errorbar(aes(ymin = ifelse(quant_val - (max_val-min_val)/2 < 0, 0, quant_val - (max_val-min_val)/2), 
        #                  ymax = quant_val + (max_val-min_val)/2), width = 0.2) +
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
                                        size = 7)) +
        
        facet_wrap(~compound_name, scales = "free_y") +
        scale_fill_manual(values = c(paletteer::paletteer_d("ggsci::default_igv", length(rvalues$quant_compounds)))) +
        xlab("\nSampleID") +
        ylab(paste0("Concentration (", input$Select_conc_unit,")")) +
        guides(fill = guide_legend(title="Compound    ")) +
        ylim(c(0, NA))

      rvalues$plot_bar1
    })
    
    output$Plot_bar2 <- renderPlot({
      rvalues$plot_bar2 <- rvalues$df_bar %>% 
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
              axis.text.x =element_text(vjust = 0.5, hjust = 1, angle = 90,
                                        size = 10)) +
        scale_fill_manual(values = c(paletteer::paletteer_d("ggsci::default_igv", length(rvalues$quant_compounds))))+
        xlab("\nSampleID") +
        ylab(paste0("Concentration (", input$Select_conc_unit,")")) +
        guides(fill = guide_legend(title="Compound    "))
      
      rvalues$plot_bar2
    })

    
  })
  
  # 4. SAVE RESULTS ####################################################################################################

  
  ## 4.3 Quant csv =====================================================================================================
  
  output$Button_download_quant_csv <- downloadHandler(
    
    filename = function(){
      paste0("quant_results_",
             gsub("\\.csv","",input$filename),"_",
             gsub("\\-","",Sys.Date()),".csv")
    },
    
    content = function(file) {
      
      rvalues$df_calibrated %>% 
        select(sampleid, compound_name, quant_val) %>% 
        mutate(sampleid = ifelse(grepl("MB|Pooled|Plasma|CC|Standard", sampleid, ignore.case = T),
                                 sampleid,
                                 gsub("^[0-9]{3}_", "", sampleid))) %>% 
        pivot_wider(names_from = compound_name, values_from = quant_val, values_fill = NA) %>% 
        write.csv(file=file, row.names=F,quote=F)
    })
  
  
  ## 4.4 Quant csv (No QCs) ============================================================================================
  
  output$Button_download_quant_csv_no_qc <- downloadHandler(
    
    filename = function(){
      paste0("removed_qcs_quant_results_",
             gsub("\\.csv","",input$filename),"_",
             gsub("\\-","",Sys.Date()),".csv")
    },
    
    content = function(file) {
      
      rvalues$df_calibrated %>% 
        select(sampleid, compound_name, quant_val) %>% 
        filter(!grepl("MB|Pooled|Plasma|CC|Standard",sampleid, ignore.case = T)) %>% 
        mutate(sampleid = gsub("^[0-9]{3}_", "", sampleid)) %>% 
        pivot_wider(names_from = compound_name, values_from = quant_val, values_fill = NA) %>% 
        write.csv(file=file, row.names = F, quote=F)
    })
  

  ## 4.5 Barplots ======================================================================================================

  output$Button_download_barplots <- downloadHandler(

    filename = function(){
      paste0(rvalues$panel, "_Barplots_",unique(rvalues$df_input$batch),"_",Sys.Date(),".pdf")
    },

    content = function(file) {
      pdf(file, height = ncol(rvalues$df_bar)*7.5 / (ncol(rvalues$df_bar)), 
          width =  nrow(rvalues$df_bar)*33 / (nrow(rvalues$df_bar))+1.5)
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
      
      ### Plot ---------------------------------------------------------------------------------------------------------
      plot_itsd <- rvalues$df_itsd_compounds %>%
        ggplot(., aes(x = num, y = peakarea, group = compound_name, label = flag)) +
        geom_point(aes(color = compound_name, shape = cc_shape), fill = "black", alpha = 0.6, size = 2) +
        geom_line(aes(y = average, color = compound_name)) +
        geom_line(aes(y = average + 1.5*stdev, color = compound_name), linetype = "dashed") +
        geom_line(aes(y = average - 1.5*stdev, color = compound_name), linetype = "dashed") +
        geom_ribbon(aes(ymin = average - stdev, ymax = average + stdev,
                        fill = compound_name), alpha=0.2) +
        ggrepel::geom_label_repel(size = 1.2, max.overlaps = Inf,
                                  min.segment.length = 0.1, label.padding = 0.1, na.rm=T) +
        theme_bw()+
        theme(panel.grid.minor= element_blank(),
              panel.grid.major.x = element_blank(),
              legend.text = element_text(color = "black", size = 8),
              legend.title = element_text(color = "black", size = 10),
              legend.position = "top",
              strip.text=element_text(color = "black", size=5),
              axis.text =element_text(color = "black", size=8),
              axis.title = element_text(color = "black", size = 10),
              plot.margin = margin(1,0.5,0,0.6, unit = 'cm')) +
        scale_fill_manual(values = c(paletteer::paletteer_d("ggsci::default_igv", length(rvalues$quant_compounds)))) +
        guides(color = guide_legend(title = "Internal Standard Compounds",
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
        facet_wrap_paginate(~ compound_name, ncol = 2, nrow = 3) +
        ggtitle(paste0(rvalues$panel, " Quantitative QC Report - ", unique(rvalues$df_input$batch)[1]))
      
      
      
      
      
      pdf(file, onefile=T, height = 11, width = 12)
      
      for(i in 1:n_pages(plot_itsd)){
        print(plot_itsd + facet_wrap_paginate(~ compound_name, ncol = 2, nrow = 3, page = i))
      }
      
      ### Summary table ------------------------------------------------------------------------------------------------
      temp <- rvalues$df_itsd_stats_compounds %>%
        dplyr::rename(Batch = batch,
                      `Internal Standard` = compound_name,
                      StDev = stdev,
                      Mean = average,
                      Median = middle,
                      CV = cv,
                      `CV Median` = cv_med) %>%
        mutate(StDev = round(StDev, digits = 0),
               Mean = round(Mean, digits = 0),
               Median = round(Median, digits = 0),
               `CV (%)` = round(CV * 100, digits = 1),
               `CV Median (%)` = round(`CV Median` * 100, digits = 1)) %>%
        select(Batch,`Internal Standard`,StDev,Mean,Median,`CV (%)`,`CV Median (%)`)


      print( gridExtra::grid.arrange(gridExtra::tableGrob(temp, rows = NULL)) )
      
      invisible(dev.off())
    })
  
  
  
}