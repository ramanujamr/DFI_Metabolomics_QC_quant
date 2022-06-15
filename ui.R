# QUANT UI

ui <- fluidPage(shinytheme("journal"),
                useShinyjs(),
                titlePanel("DFI Metabolomics QC - Quant V2"),
                
                fluidRow(
                  h3("1. UPLOAD DATA"),
                  hr(),
                  column(width=3,
                  selectInput("filename", label="Select csv file to upload",
                              rev(list.files(wddir, 
                                             pattern = ".*bile.*csv|.*PFBBr.*csv|.*Indole.*csv|.*Tryptophan.*csv",
                                             ignore.case = T))),
                  actionButton("Button_refresh_csv", "Refresh", icon("sync"), width="100px"),
                  actionButton("Button_upload_csv", "Upload", icon("upload"), width="100px",
                                               style="color: #fff; background-color: #2346b0; border-color: #2e6da4"),
                  textOutput("Textout_filename"),
                  textOutput("Textout_panel"),
                  actionButton("Button_new_file", "Upload New File", icon("upload"), width="150px",
                               style="color: #fff; background-color: #2346b0; border-color: #2e6da4")),

                  column(width=3, offset=1,
                         br(),
                         fluidRow(actionButton("Button_itsd_stats", "ITSD stats", icon("tv"), width="150px")),
                         bsModal("Modal_ITSD_stats", "ITSD stats", "Button_itsd_stats", size = "large", DT::dataTableOutput("Table_ITSD_stats"))),
                         
                  column(width=1, shinythemes::themeSelector())
                ),
                
                hr(),
                fluidRow(h3("2. SETTINGS")),
                br(),
                fluidRow(
                  column(width=9, rHandsontableOutput("Table_calibration_settings")),
                  column(width=2, offset=1, align="center",
                         fluidRow(numericInput(inputId = "Textin_x_factor", value = 1, label="Multiplication factor"),
                                  selectInput("Select_conc_unit","Select concentration unit", choices = c("ug/mL", "mg/mL", "uM/mL")),
                                  checkboxInput("Checkbox_intercept_zero","Force intercept to zero (Tryptophan only)", value=F))
                         )),
                
                br(), br(),
                fluidRow(column(width=3, offset=9, align="center", 
                                shiny::actionButton("Button_calibration", "Calibration", icon("chevron-right"), width="200px",
                                                      style="color: #fff; background-color: #2346b0; border-color: #2e6da4"))
                         ),
                
                br(), hr(), br(),
                fluidRow(h3("3. CALIBRATION")),
                br(),
                fluidRow(
                  h4("3.1 Calibration plots"),
                  column(width=12,
                  plotOutput("Plot_calibration_curves", 
                             click = "Plot_calibration_curves_click",
                             dblclick = "Plot_calibration_curves_dblclick",
                             brush = brushOpts(id = "Plot_calibration_curves_brush", resetOnNew = TRUE)))
                  ),
                fluidRow(
                  column(width=2, numericInput("Numericinput_min_cc", "Min CC", value=0)),
                  column(width=2, numericInput("Numericinput_max_cc", "Max CC", value=125))
                ),
                fluidRow( column(width=1, shiny::actionButton("Button_cc_range", "Update"))),
                br(), br(), hr(),
                
                fluidRow(
                  column(width=6, 
                         h4("3.2 Calibration data"),
                         DT::dataTableOutput("Table_cc_data")),
                  column(width=4, offset=1, 
                         h4("3.3 Linear models"),
                         DT::dataTableOutput("Table_linear_models"))
                ),
                br(),
                fluidRow(
      
                  column(width=6, align="center", shiny::actionButton("Button_save_cc_metrics", "Save CC metrics", 
                                                                      icon("save"), width="200px",
                                                                      style="color: #fff; background-color: #00ab66; border-color: #2e6da4")),
                  column(width=3, offset=3, align="center", 
                         checkboxInput("Checkbox_subtract_MB","Subtract Method Blanks", value=T),
                         shiny::actionButton("Button_continue", "Continue", icon("chevron-right"), width="200px",
                                                                      style="color: #fff; background-color: #2346b0; border-color: #2e6da4"))

                  ),
                
                br(), br(), hr(),
                
                fluidRow(h3("4. CALIBRATED RESULTS")),
                br(),
                
                fluidRow(column(width=6, align="center", 
                                h4("4.1 Plasma QC"),
                         plotOutput("Plot_plasma_qc")),
                         column(width=6, align="center",
                                h4("4.2 Method blanks"),
                         plotOutput("Plot_method_blanks"))),
                
                br(), br(), hr(),
                
                
                fluidRow(
                  h4("4.3 Calibrated data table"),
                  column(width=12, DT::dataTableOutput("Table_calibrated_data"))),
                  
                  br(), br(), hr(),
                
                fluidRow(  
                  h4("4.4 Barplots (by compounds)"),
                  h5("Error bars are calculated by the range of plasma QCs (max-min) for each compound"),
                  plotOutput("Plot_bar1", height = "800px")),
                  
                  br(), br(), hr(),
                fluidRow(
                  h4("4.5 Barplots (overall)"),
                  plotOutput("Plot_bar2", height = "800px")),
                
                br(), br(),

                
                fluidRow(
                  column(width=3, align="center", 
                         fluidRow(
                           shiny::downloadButton("Button_download_quant_csv", "Quant Results", icon("file-csv"), width="200px",
                                               style="color: #fff; background-color: #00ab66; border-color: #2e6da4"),
                           shiny::downloadButton("Button_download_quant_csv_no_qc", "Quant Results (No QC)", icon("file-csv"), width="200px",
                                               style="color: #fff; background-color: #00ab66; border-color: #2e6da4")
                         )),
                  
                  column(width=3, align="center", shiny::downloadButton("Button_download_barplots", "Barplots",
                                                                      icon("file-pdf"), width="200px",
                                                                      style="color: #fff; background-color: #00ab66; border-color: #2e6da4")),
                  
                  column(width=3, align="center", shiny::downloadButton("Button_download_qc_report", "QC Report",
                                                                      icon("file-pdf"), width="200px",
                                                                      style="color: #fff; background-color: #00ab66; border-color: #2e6da4")),
                
                br(), br(), br()
                )
                
)