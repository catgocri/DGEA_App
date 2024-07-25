library(shiny)
library(shinyFiles)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(tidyverse)
library(limma)
library(EnhancedVolcano)
library(clipr)
library(ggfortify)
library(plotly)
library(processx)
library(RColorBrewer)
library(DT)
library(grid)

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      .container-fluid, .row, .col-md-4, .col-md-8 {
        height: 100%;
      }
      .navbar {
        margin-bottom: 0;
      }
      .tab-content {
        height: calc(100vh - 56px); /* Adjust height for the navbar */
        overflow-y: auto;
      }
      .plot-container {
        height: calc(100vh - 56px); /* Adjust height for the navbar */
        display: flex;
        align-items: center;
        justify-content: center;
        overflow: auto;
      }
      .plot-container .plot {
        max-width: 100%;
        max-height: 100%;
      }
    "))
  ),
  titlePanel("Pipeline Starter & DGE Analysis"),
  
  # Define a navbarPage with tab panels
  navbarPage("Pipeline & Plotting",
             tabPanel("Pipeline Settings",
                      sidebarLayout(
                        sidebarPanel(
                          tags$style(type='text/css', '#nfcore_command {white-space: pre-wrap; word-break: keep-all;}'),
                          selectInput("terminal", "Select Terminal (For viewing of nextflow processes)", 
                                      choices = list("Default for OS" = "", "cmd (Windows)" = "cmd", "PowerShell (Windows)" = "powershell", 
                                                     "Terminal (macOS)" = "terminal", "iTerm (macOS)" = "iterm", 
                                                     "GNOME Terminal (Linux)" = "gnome-terminal", "xterm (Linux)" = "xterm", 
                                                     "Konsole (Linux)" = "konsole")),
                          shinyDirButton("samples", "Select sample directory...", "Select sample directory for sheetMaker.sh"),
                          verbatimTextOutput("samples", placeholder = T),
                          shinyDirButton("outdir", "Select pipeline output directory...", "Select output directory for pipeline"),
                          verbatimTextOutput("outdir", placeholder = T),
                          fileInput("samplesheet", "Select sample annotations file", accept = ".tsv"),
                          actionButton("load_sheets", "Generate projSheet.csv Data"),
                          actionButton("load_existing_sheets", "Load existing projSheet.csv Data"),
                          actionButton("clear_term", "Clear Log"),
                          textOutput("divider"),
                          verbatimTextOutput("nfcore_command", placeholder = T),
                          selectInput("genome", "Select reference genome for transcript alignment.", choices=list("GRCh37" = "GRCh37", "GRCh38" = "GRCh38")),
                          checkboxInput("save_trimmed", "Save Trimmed", T),
                          checkboxInput("save_non_ribo_reads", "Save Non Ribosomal Reads", T),
                          checkboxInput("remove_ribo_rna", "Remove rRNA (RiboDetector)", T),
                          checkboxInput("save_unaligned", "Save unaligned", T),
                          checkboxInput("skip_gtf_filter", "Skip GTF filter", F),
                          checkboxInput("skip_gtf_transcript_filter", "Skip GTF transcript filter", F),
                          checkboxInput("skip_bbsplit", "Skip BBSplit", F),
                          checkboxInput("skip_umi_extract", "Skip UMI extraction", F),
                          checkboxInput("skip_trimming", "Skip trimming", F),
                          checkboxInput("skip_alignment", "Skip alignment", F),
                          checkboxInput("skip_pseudo_alignment", "Skip pseudo alignment", F),
                          checkboxInput("skip_markduplicates", "Skip MarkDuplicates", F),
                          checkboxInput("skip_bigwig", "Skip bigWig creation", F),
                          checkboxInput("skip_stringtie", "Skip StringTie", F),
                          checkboxInput("skip_fastqc", "Skip FastQC", F),
                          checkboxInput("skip_preseq", "Skip Preseq", F),
                          checkboxInput("skip_dupradar", "Skip dupRadar", F),
                          checkboxInput("skip_qualimap", "Skip Qualimap", F),
                          checkboxInput("skip_rseqc", "Skip RSeQC", F),
                          checkboxInput("skip_biotype_qc", "Skip biotype QC", F),
                          checkboxInput("skip_deseq2_qc", "Skip DESeq2 QC", F),
                          checkboxInput("skip_multiqc", "Skip MultiQC", F),
                          checkboxInput("skip_qc", "Skip all QC except MultiQC", F),
                          actionButton("run_pipeline", "Run nf-core/RNA-Seq Pipeline"),
                        ),
                        mainPanel(
                          tabsetPanel(
                            tabPanel("sheetMaker.sh CSV Preview", tableOutput("sheetMakerCsv")),
                            tabPanel("sheetMaker2.sh CSV Preview", tableOutput("sheetMaker2Csv")),
                            tabPanel("Samplesheet TSV Preview", tableOutput("samplesheetTsv")),
                            tabPanel("MultiQC Finder Pre-Alignment", uiOutput("multiqcViewerPre")),
                            tabPanel("MultiQC Finder Post-Alignment", uiOutput("multiqcViewer")),
                            tabPanel("Count Data Preview", tableOutput("count_data_pre")),
                            tabPanel("Log Messages", div(style="overflow-y:scroll;", verbatimTextOutput("logging")))
                          )
                        )
                      )
             ),
             tabPanel("Plotting Settings",
                      sidebarLayout(
                        sidebarPanel(
                          fileInput("load_dds", "Load saved DDS object"),
                          selectInput("intGroupSelect", "Select Comparison Group", choices = NULL),
                          selectInput("subGroupSelect", "Select Subgroup", choices = NULL),
                          sliderInput("plot_width", "Plot Width (px):", min = 100, max = 20000, value = 800),
                          sliderInput("plot_height", "Plot Height (px):", min = 100, max = 20000, value = 600),
                          numericInput("pvalue_cutoff", "(Volcano) P-value Cutoff:", value = 0.05, min = 0, max = 0.1, step = 0.0001, width = '200px'),
                          numericInput("fc_cutoff", "(Volcano) Log2 Fold Change Cutoff:", value = 1, min = 0, max = 20, step = 0.1, width = '200px'),
                          checkboxInput("full_data", "(Volcano) Show Full Data", value = FALSE),
                          actionButton("save_dds", "Save currently loaded DDS object"),
                        ),
                        mainPanel(
                          tabsetPanel(
                            tabPanel("MA Plot", div(class = "plot-container", plotOutput("maPlot", height = "100%"))),
                            tabPanel("PCA Plot", div(class = "plot-container", plotlyOutput("pcaPlot", height = "100%"))),
                            tabPanel("Cook's Distance Plot", div(class = "plot-container", plotOutput("cooksPlot", height = "100%"))),
                            tabPanel("MDS Plot", div(class = "plot-container", plotlyOutput("mdsPlot", height = "100%"))),
                            tabPanel("Volcano Plot", div(class = "plot-container", plotlyOutput("volcanoPlot", height = "100%"), )),
                            tabPanel("HC Heatmap Plot Top 50 DEG's", div(class = "plot-container", plotOutput("hcHeatmapPlot", height = "100%"))),
                            tabPanel("HC Heatmap Plot Sample Distances", div(class = "plot-container", plotOutput("hcHeatmapPlot2", height = "100%"))),
                            tabPanel("DESeq2 Analysis Summary", verbatimTextOutput("summary"), verbatimTextOutput("resname"), verbatimTextOutput("ddsview"), verbatimTextOutput("vsdview")),
                          ),
                          DTOutput("plotMeta"),
                        )
                      )
             )
  )
)

server <- function(input, output, session) {
  
  # Max input size 30MB
  options(shiny.maxRequestSize=100000*1024^2)
  
  ### START OF PIPELINE INPUT SELECTION ###
  # Define starting directories for both directory inputs
  shinyDirChoose(input, 'samples', roots = c(root = '/'))
  shinyDirChoose(input, 'outdir', roots = c(root = '/'))
  
  # Where to store the directories
  
  global <- reactiveValues(
    madeDds = F, 
    samples_dir = NULL, 
    outdir = NULL, 
    sheetMaker_file_dir = NULL, 
    sheetMaker2_file_dir = NULL, 
    nfcore_command_text = NULL, 
    counts_filepath = NULL,
    dds = NULL,
    vsd = NULL,
    result_names = NULL,
    result_summary = NULL,
    hover_metadata = NULL,
    pca_results = NULL,
    volcano_df = NULL
  )
  
  save_deseq_object <- function(deseq_object) {
    # Get current UNIX timestamp
    timestamp <- as.numeric(Sys.time())
    
    filename <- paste0("deseq_object_", timestamp, ".rds")
    
    saveRDS(deseq_object, file = file.path(global$outdir, filename))
    return(filename)
  }
  
  log_messages <- reactiveValues(value = "")
  append_to_log <- function(input) {
    log_messages$value <- paste(log_messages$value, paste(input, collapse = "\n"), sep = "\n")
  }
  
  
  
  # Get the directories from the input
  samples <- reactive(input$samples)
  outdir <- reactive(input$outdir)
  samplesheet <- reactive({
    if (is.null(input$samplesheet)) {
      return("")
    }
    read_tsv(file = input$samplesheet$datapath)
  })
  
  # What to do with the output
  output$samples <- renderText({
    global$samples_dir
  })
  
  output$outdir <- renderText({
    global$outdir
  })
  
  # Wait for inputs and render the path 
  observeEvent(input$samples, {
    if (!"path" %in% names(samples())) return()
    global$samples_dir <- file.path("/", paste(unlist(samples()$path[-1]), collapse = .Platform$file.sep))
    output$samples <- renderText({ global$samples_dir })
  })
  
  observeEvent(input$outdir, {
    if (!"path" %in% names(outdir())) return()
    global$outdir <- file.path("/", paste(unlist(outdir()$path[-1]), collapse = .Platform$file.sep))
    global$counts_filepath <- paste0(global$outdir, "/star_salmon/salmon.merged.gene_tpm.tsv")
    output$outdir <- renderText({ global$outdir })
  })
  
  ### START AUTOMATIC DESEQ ANALYSIS ###
  observe({
    req(global$counts_filepath, file.exists(global$counts_filepath), input$samplesheet)
    if (isFALSE(global$madeDds)) {
      showModal(modalDialog("Count data found, loading count preview and starting DESeq2 analysis. Please wait if the program freezes as this may take a while."))
      global$madeDds <- T
      counts <- read_tsv(global$counts_filepath)
      
      # Check for duplicates in gene names
      if (any(duplicated(counts$gene_name))) {
        append_to_log("Duplicate gene names found. Using Ensembl IDs as row indices.")
        counts <- counts %>% 
          select(-gene_name) %>%
          column_to_rownames("gene_id")
      } else {
        counts <- counts %>% 
          select(-gene_id) %>%
          column_to_rownames("gene_name")
      }
      
      counts <- counts %>% 
        mutate(across(everything(), ~ round(as.numeric(as.character(.)), 0)))
      
        append_to_log("Generating DESeqDataSet")
        #coldata <- samplesheet() %>% column_to_rownames("Delivery_name")
        dds <- DESeqDataSetFromMatrix(countData = counts, colData = samplesheet(), design = ~ cases)
        dds <- DESeq(dds)
        global$dds <- dds
        global$vsd <- vst(dds, blind=F)
        global$results <- capture.output(results(dds))
        global$result_names <- resultsNames(dds)
        global$result_summary <- capture.output(summary(results(dds)))
        append_to_log("DESeqDataSet created successfully.")
    }
  })
  ### END AUTOMATIC DESEQ ANALYSIS ###
  ### END OF PIPELINE INPUT SELECTION ###
  
  
  
  
  
  
  
  ### START PIPELINE CODE AND LOGGING FUNCTIONALITY ###
  
  temp_log_file <- reactiveVal(NULL)
  
  run_command <- function(input, temp_file) {
    system(paste(input, "1>", shQuote(temp_file), "2>&1"))
  }
  
  observeEvent(input$load_sheets, {
    if (!is.null(global$samples_dir) && !is.null(global$outdir)) {

      command <- paste0("./sheetMaker.sh ", global$samples_dir, " DGEA_projSheet ", global$outdir)
      
      # Create a temporary file to store the output
      temp_file <- tempfile(fileext = ".txt")
      temp_log_file(temp_file)
      
      # Execute the command and redirect output to the temporary file
      run_command(command, temp_file)
      
      # Read the contents of the temporary file
      if (file.exists(temp_file)) {
        log_content <- readLines(temp_file)
        append_to_log(log_content)
      } else {
        append_to_log("Error: Failed to read log file.")
      }
      
      global$sheetMaker_file_dir <- file.path(global$outdir, "sheetMaker", "DGEA_projSheet.csv")

      if (file.exists(global$sheetMaker_file_dir)) {
        append_to_log("Read sheetMaker csv file.")
      }
      else {
        append_to_log("Error: Failed to render sheetMaker csv file.")
      }
      
      # Remove temporary file
      unlink(temp_file)
    } else {
      if (is.null(global$samples_dir) || is.null(global$outdir)) {
        append_to_log("Error: Directories are not selected.")
      }
      if (is.null(input$samplesheet)) {
        append_to_log("Error: Sample annotations file not selected.")
      }
    }
  })
  
  observeEvent(input$load_existing_sheets, {
    if (!is.null(global$outdir)) {
      if (file.exists(file.path(global$outdir, "sheetMaker", "DGEA_projSheet.csv"))) {
        global$sheetMaker_file_dir <- file.path(global$outdir, "sheetMaker", "DGEA_projSheet.csv")
        
        if (file.exists(global$sheetMaker_file_dir)) {
          append_to_log("Read sheetMaker csv file.")
        }
        else {
          append_to_log("Error: Failed to render sheetMaker csv file.")
        }
      }
      if (file.exists(file.path(global$outdir, "sheetMaker", "DGEA_projSheet2.csv"))) {
        global$sheetMaker2_file_dir <- file.path(global$outdir, "sheetMaker", "DGEA_projSheet2.csv")
        
        if (file.exists(global$sheetMaker2_file_dir)) {
          append_to_log("Read sheetMaker2 csv file.")
        }
        else {
          append_to_log("Error: Failed to render sheetMaker2 csv file.")
        }
      }
    } else {
      append_to_log("Error: Output directory not selected to read existing sheetMaker csv files.")
    }
  })
  
  observeEvent(input$clear_term, {
    log_messages$value <- ""
  })
  
  output$logging <- renderText({log_messages$value})
    
  output$sheetMakerCsv <- renderTable({
    if (!is.null(global$sheetMaker_file_dir) && file.exists(global$sheetMaker_file_dir)) {
      read.csv(global$sheetMaker_file_dir)
    } else {
      data.frame("N/A")
    }
  })
  
  output$count_data_pre <- renderTable({
    if (!is.null(global$counts_filepath) && file.exists(global$counts_filepath)) {
      read_tsv(global$counts_filepath)
    } else {
      data.frame("N/A")
    }
  })
  
  output$sheetMaker2Csv <- renderTable({
    if (!is.null(global$sheetMaker2_file_dir) && file.exists(global$sheetMaker2_file_dir)) {
      read.csv(global$sheetMaker2_file_dir)
    } else {
      data.frame("N/A")
    }
  })
  
  observeEvent(input$samplesheet, {
    append_to_log("Sample annotations uploaded successfully.")
  })
  
  output$samplesheetTsv <- renderTable({
    if (!is.null(input$samplesheet)) {
      samplesheet()
    } else {
      data.frame("N/A")
    }
  })
  
  output$divider <- renderText({
    "Command Builder:"
  })
  
  ## START PIPELINE
  
  open_terminal <- function(terminal_app) {
    os_type <- .Platform$OS.type
    sys_name <- Sys.info()["sysname"]
    
    if (os_type == "windows") {
      if (tolower(terminal_app) == "cmd" || terminal_app == "") {
        system("cmd.exe")
      } else if (tolower(terminal_app) == "powershell") {
        system("powershell.exe")
      } else {
        append_to_log("Unsupported terminal application for Windows.")
      }
    } else if (sys_name == "Darwin") { # macOS
      if (tolower(terminal_app) == "terminal" || terminal_app == "") {
        system("open -a Terminal")
      } else if (tolower(terminal_app) == "iterm") {
        system(paste0("open -a iTerm"))
      } else {
        append_to_log("Unsupported terminal application for macOS.")
      }
    } else if (sys_name == "Linux") { # Linux
      if (tolower(terminal_app) == "gnome-terminal" || terminal_app == "") {
        system("gnome-terminal")
      } else if (tolower(terminal_app) == "xterm") {
        system("xterm")
      } else if (tolower(terminal_app) == "konsole") {
        system("konsole")
      } else {
        append_to_log("Unsupported terminal application for Linux.")
      }
    } else {
      append_to_log("Unsupported operating system.")
    }
  }
  
  output$nfcore_command <- renderText({
    if (input$remove_ribo_rna) {
      global$nfcore_command_text <- paste0("./pipelineStarter.sh ", global$outdir, " ", shQuote(paste0("nextflow run nf-core/rnaseq -r 3.14.0 -profile singularity --input ", global$outdir, "/sheetMaker/DGEA_projSheet.csv --outdir ", global$outdir, " -work-dir ", global$outdir, "/work ", "--genome ", input$genome, " --save_trimmed ", input$save_trimmed, " --save_non_ribo_reads ", input$save_non_ribo_reads, " --remove_ribo_rna FALSE", " --skip_gtf_filter ", input$skip_gtf_filter, " --skip_gtf_transcript_filter ", input$skip_gtf_transcript_filter, " --skip_bbsplit ", input$skip_bbsplit, " --skip_umi_extract ", input$skip_umi_extract, " --skip_trimming ", input$skip_trimming, " --skip_alignment TRUE", " --skip_pseudo_alignment TRUE", " --skip_markduplicates TRUE", " --skip_bigwig TRUE", " --skip_stringtie TRUE", " --skip_fastqc ", input$skip_fastqc, " --skip_preseq ", input$skip_preseq, " --skip_dupradar TRUE", " --skip_qualimap TRUE", " --skip_rseqc TRUE", " --skip_biotype_qc TRUE", " --skip_deseq2_qc TRUE", " --skip_multiqc ", input$skip_multiqc, " --skip_qc ", input$skip_qc)), " ", shQuote(paste0("nextflow run nf-core/rnaseq -r 3.14.0 -resume -profile singularity --input ", global$outdir, "/sheetMaker/DGEA_projSheet2.csv --outdir ", global$outdir, " -work-dir ", global$outdir, "/work ", "--genome ", input$genome, " --save_trimmed ", input$save_trimmed, " --save_non_ribo_reads ", input$save_non_ribo_reads, " --save_unaligned ", input$save_unaligned, " --skip_gtf_filter TRUE", " --skip_gtf_transcript_filter TRUE", " --skip_bbsplit TRUE", " --skip_umi_extract TRUE", " --skip_trimming TRUE", " --skip_alignment ", input$skip_alignment, " --skip_pseudo_alignment ", input$skip_pseudo_alignment, " --skip_markduplicates ", input$skip_markduplicates, " --skip_bigwig ", input$skip_bigwig, " --skip_stringtie ", input$skip_stringtie, " --skip_fastqc ", input$skip_fastqc, " --skip_preseq ", input$skip_preseq, " --skip_dupradar ", input$skip_dupradar, " --skip_qualimap ", input$skip_qualimap, " --skip_rseqc ", input$skip_rseqc, " --skip_biotype_qc ", input$skip_biotype_qc, " --skip_deseq2_qc ", input$skip_deseq2_qc, " --skip_multiqc ", input$skip_multiqc, " --skip_qc ", input$skip_qc)))
    } else {
      global$nfcore_command_text <- paste0("nextflow run nf-core/rnaseq -r 3.14.0 -profile singularity --input ", global$outdir, "/sheetMaker/DGEA_projSheet.csv --outdir ", global$outdir, " -work-dir ", global$outdir, "/work ", "--genome ", input$genome, " --save_trimmed ", input$save_trimmed, " --save_non_ribo_reads ", input$save_non_ribo_reads, " --remove_ribo_rna ", input$remove_ribo_rna, " --save_unaligned ", input$save_unaligned, " --skip_gtf_filter ", input$skip_gtf_filter, " --skip_gtf_transcript_filter ", input$skip_gtf_transcript_filter, " --skip_bbsplit ", input$skip_bbsplit, " --skip_umi_extract ", input$skip_umi_extract, " --skip_trimming ", input$skip_trimming, " --skip_alignment ", input$skip_alignment, " --skip_pseudo_alignment ", input$skip_pseudo_alignment, " --skip_markduplicates ", input$skip_markduplicates, " --skip_bigwig ", input$skip_bigwig, " --skip_stringtie ", input$skip_stringtie, " --skip_fastqc ", input$skip_fastqc, " --skip_preseq ", input$skip_preseq, " --skip_dupradar ", input$skip_dupradar, " --skip_qualimap ", input$skip_qualimap, " --skip_rseqc ", input$skip_rseqc, " --skip_biotype_qc ", input$skip_biotype_qc, " --skip_deseq2_qc ", input$skip_deseq2_qc, " --skip_multiqc ", input$skip_multiqc, " --skip_qc ", input$skip_qc)
    }
  })
  
  observeEvent(input$run_pipeline, {
    if (!is.null(input$samplesheet) && !is.null(global$sheetMaker_file_dir) && file.exists(global$sheetMaker_file_dir)) {
      append_to_log("Launching CLI for logging of nf-core/RNA-Seq. Remember to copy the command in the command builder into the CLI.")
      showNotification("Copied command to clipboard and opening CLI...", closeButton=T)
      Sys.sleep(1)
      open_terminal(input$terminal)
    } else {
      append_to_log("Error: Required files are not uploaded for start of pipeline.")
    }
  })
  
  output$summary <- renderText({
    if (!is.null(global$dds)) {
      paste(global$result_summary, collapse = "\n")
    } else {
      "No DESeqDataSet Found."
    }
  })
  
  output$resname <- renderText({
    if (!is.null(global$dds)) {
      paste(global$result_names, collapse="\n")
    } else {
      "No DESeqDataSet Found."
    }
  })
  
  output$ddsview <- renderText({
    if (!is.null(global$dds)) {
      paste(global$results, collapse = "\n")
    } else {
      "No DESeqDataSet Found."
    }
  })
  
  output$vsdview <- renderText({
    if (!is.null(global$vsd)) {
      paste(capture.output(global$vsd), collapse = "\n")
    } else {
      "No Variance Stabilized DESeqDataSet Found."
    }
  })
  
  output$multiqcViewerPre <- renderUI({
    if (!is.null(global$outdir) && file.exists(file.path(global$outdir, "multiqc/multiqc_report.html"))) {
      tags$div(
        h4("MultiQC report is ready."),
        p("Please open the following file in your web browser:"),
        p(file.path(global$outdir, "multiqc/multiqc_report.html"))
      )
    } else {
      h4("MultiQC report not found.")
    }
  })
  
  
  output$multiqcViewer <- renderUI({
    if (!is.null(global$outdir) && file.exists(file.path(global$outdir, "multiqc/star_salmon/multiqc_report.html"))) {
      tags$div(
        h4("MultiQC report is ready."),
        p("Please open the following file in your web browser:"),
        p(file.path(global$outdir, "multiqc/star_salmon/multiqc_report.html"))
      )
    } else {
      h4("MultiQC report not found.")
    }
  })
  
  ### PLOTS
  
  plotWidth <- reactive({
    input$plot_width
  })
  
  plotHeight <- reactive({
    input$plot_height
  })
  
  output$maPlot <- renderPlot({
    req(global$dds)
    DESeq2::plotMA(global$dds, ylim=c(-30,30), alpha=0.1)
  }, width = plotWidth, height = plotHeight)
  # 
  output$pcaPlot <- renderPlotly({
    req(global$vsd, input$intGroupSelect, input$subGroupSelect)
    
    pca_results <- plotPCA(global$vsd, intgroup = input$intGroupSelect, returnData = TRUE)
    
    # Add Delivery_name column to pca_results if not present
    pca_results$Delivery_name <- rownames(pca_results)
    
    global$pca_results <- pca_results
    
    p <- plot_ly(pca_results, x = ~PC1, y = ~PC2, color = ~get(input$intGroupSelect),
                 type = 'scatter', mode = 'markers',
                 marker = list(size = 10),
                 source = "pca_plot",
                 customdata = ~Delivery_name, 
                 text = ~paste('Sample:', Delivery_name, '<br>',
                               'PC1:', round(PC1, 2), '<br>',
                               'PC2:', round(PC2, 2)),
                 hoverinfo = 'text') %>%
      layout(title = 'PCA Plot',
             xaxis = list(title = 'PC1'),
             yaxis = list(title = 'PC2'))
    p <- event_register(p, "plotly_selected")
    p
  })
  
  observeEvent(event_data("plotly_selected", source = "pca_plot"), {
    selected_points <- event_data("plotly_selected", source = "pca_plot")
    
    # Debugging
    print("Selected Points Data:")
    print(selected_points)
    
    if (!is.null(selected_points) && length(selected_points$pointNumber) > 0) {
      # Extract indices of selected points
      selected_samples <- selected_points$customdata
      selected_indices <- selected_points$pointNumber
      
      # Debugging
      print("Selected Samples:")
      print(selected_samples)
      
      # Check if indices are valid
      if (all(selected_indices >= 0 & selected_indices < nrow(global$pca_results))) {
        # Map samples to PCA
        selected_pca_data <- global$pca_results[selected_samples, ]
        
        # Debugging
        print("Selected PCA Data:")
        print(selected_pca_data)
        
        # Extract Delivery_name for metadata lookup
        selected_names <- selected_pca_data$Delivery_name
        print("Selected Names Data:")
        print(selected_names)
        
        #  metadata from colData
        metadata <- as.data.frame(colData(global$dds))
        
        # Debugging
        print("Metadata Data:")
        print(metadata)
        
        # Ensure rownames of metadata match
        rownames(metadata) <- metadata$Delivery_name
        selected_metadata <- metadata[rownames(metadata) %in% selected_names, ]
        
        # Debugging
        print("Selected Metadata:")
        print(selected_metadata)
        
        # Display the selected metadata in the table
        output$plotMeta <- renderDT({
          datatable(selected_metadata)
        })
      } else {
        # Handle case out of bounds
        output$plotMeta <- renderDT({
          datatable(data.frame())
        })
      }
    } else {
      # Handle case when no points are selected
      output$plotMeta <- renderDT({
        datatable(data.frame())
      })
    }
  })
  
  output$cooksPlot <- renderPlot({
    req(global$dds)
    # Plot config.
    par(mar=c(5,5,1,1))
    # Plot cooks distances.
    boxplot(log10(assays(global$dds)[["cooks"]]), range=0, las=2)
  }, width = plotWidth, height = plotHeight)
  
  output$mdsPlot <- renderPlotly({
    req(global$vsd, input$intGroupSelect, input$subGroupSelect)
    
    # Perform MDS
    dist_matrix <- dist(t(assay(global$vsd)))
    mds_data <- cmdscale(dist_matrix, eig=T, k = 2)
    mds_df <- data.frame(
      MDS1 = mds_data$points[,1], 
      MDS2 = mds_data$points[,2], 
      group = colData(global$vsd)[[input$intGroupSelect]],
      Delivery_name = rownames(colData(global$vsd)) # Add Delivery_name column
    )

    global$mds_results <- mds_df
    
    p <- plot_ly(mds_df, x = ~MDS1, y = ~MDS2, color = ~group,
                 type = 'scatter', mode = 'markers',
                 marker = list(size = 10),
                 customdata = ~Delivery_name,
                 text = ~paste('Sample:', Delivery_name, '<br>',
                               'MDS1:', round(MDS1, 2), '<br>',
                               'MDS2:', round(MDS2, 2)),
                 hoverinfo = 'text',
                 source = "mds_plot") %>%
      layout(title = 'MDS Plot',
             xaxis = list(title = 'Dimension 1'),
             yaxis = list(title = 'Dimension 2'))
    p <- event_register(p, "plotly_selected")
    p
  })
  
  observeEvent(event_data("plotly_selected", source = "mds_plot"), {
    selected_points <- event_data("plotly_selected", source = "mds_plot")
    
    # Debugging
    print("Selected Points Data:")
    print(selected_points)
    
    if (!is.null(selected_points) && length(selected_points$pointNumber) > 0) {
      selected_samples <- selected_points$customdata
      selected_indices <- selected_points$curveNumber
      
      # Debuggin
      print("Selected Samples:")
      print(selected_samples)
      
      if (all(selected_indices >= 0 & selected_indices < nrow(global$mds_results))) {
        selected_mds_data <- global$mds_results[selected_samples, ]
        
        # Debugging
        print("Selected MDS Data:")
        print(selected_mds_data)
        
        # Extract Delivery_name for metadata lookup
        selected_names <- selected_mds_data$Delivery_name
        print("Selected Names Data:")
        print(selected_names)
        
        # metadata from colData againnn
        metadata <- as.data.frame(colData(global$dds))
        
        # Debugging
        print("Metadata Data:")
        print(metadata)
        
        rownames(metadata) <- metadata$Delivery_name
        selected_metadata <- metadata[rownames(metadata) %in% selected_names, ]
        
        # Debugging
        print("Selected Metadata:")
        print(selected_metadata)
        
        # Display the selected metadata in the table
        output$plotMeta <- renderDT({
          datatable(selected_metadata)
        })
      } else {
        # out of bounds
        output$plotMeta <- renderDT({
          datatable(data.frame())
        })
      }
    } else {
      # Handle case when no points are selected
      output$plotMeta <- renderDT({
        datatable(data.frame())
      })
    }
  })
  
  output$volcanoPlot <- renderPlotly({
    req(global$results)
    
    cb_palette <- brewer.pal(n = 4, name = "Dark2")
    
    # global$results is not the data but a logging of the output of it so we cant use that.
    res <- results(global$dds)
    
    res$gene <- rownames(res)
    
    # Classify genes based on input cutoffs
    res$Significance <- with(res, ifelse(pvalue < input$pvalue_cutoff & abs(log2FoldChange) > input$fc_cutoff, "Log2FC & pvalue",
                                       ifelse(pvalue < input$pvalue_cutoff, "pvalue",
                                              ifelse(abs(log2FoldChange) > input$fc_cutoff, "Log2FC", "Non-significant"))))
    
    # Separate significant and non-significant points
    significant_points <- res[res$Significance %in% c("Log2FC & pvalue", "pvalue", "Log2FC"), ]
    non_significant_points <- res[res$Significance == "Non-significant", ]
    
    # Downsampling
    if (!input$full_data) {
      set.seed(123)
      non_significant_points <- non_significant_points[sample(nrow(non_significant_points), min(5000, nrow(non_significant_points))), ]
    }
    
    # Combine significant and sampled non-significant points
    res <- rbind(significant_points, non_significant_points)
    
    global$volcano_df <- res
    
    p <- plot_ly(data = res, x = ~log2FoldChange, y = ~-log10(pvalue), text = ~paste("Gene:", gene), customdata = ~gene, 
                 color = ~Significance, colors = cb_palette, type = 'scattergl', mode = 'markers',
                 marker = list(size = 5, opacity = 0.7, line = list(width = 0.5, color = 'black')),
                 source = "volcano_plot") %>%
      layout(
        title = "Volcano Plot",
        xaxis = list(title = "log2 Fold Change"),
        yaxis = list(title = "-log10 p-value"),
        showlegend = TRUE
      )
    
    p <- event_register(p, "plotly_selected")
    p
  })
  
  observeEvent(event_data("plotly_selected", source = "volcano_plot"), {
    selected_points <- event_data("plotly_selected", source = "volcano_plot")
    
    # Debugging
    print("Selected Points Data:")
    print(selected_points)
    
    if (!is.null(selected_points) && length(selected_points$pointNumber) > 0) {
      selected_indices <- selected_points$pointNumber
      selected_genes <- selected_points$customdata
      
      # Debugging
      print("Selected Genes:")
      print(selected_genes)
      
      df <- global$volcano_df

      if (all(selected_indices >= 0 & selected_indices < nrow(res))) {
        selected_volcano_data <- res[selected_genes, ] 
        
        # Debugging
        print("Selected Volcano Data:")
        print(selected_volcano_data)
        
        output$plotMeta <- renderDT({
          datatable(selected_volcano_data)
        })
      } else {
        output$plotMeta <- renderDT({
          datatable(data.frame())
        })
      }
    } else {
      # Handle case when no points are selected
      output$plotMeta <- renderDT({
        datatable(data.frame())
      })
    }
  })
  
  hcHeatmapPlot2 <- reactive({
    req(global$vsd, input$intGroupSelect)

    comparisonGroup <- colData(global$vsd)[[input$intGroupSelect]]
    
    filtered_value <- input$subGroupSelect
    selected_rows <- comparisonGroup == filtered_value
    
    if (!any(selected_rows)) {
      stop("No rows selected. Check if comparisonGroup contains the value.")
    }
    
    # Calculate distance matrix and subset
    dist_matrix_2 <- as.matrix(dist(t(assay(global$vsd))))
    dist_matrix_subset <- dist_matrix_2[selected_rows, selected_rows]

    # Create heatmap
    pheatmap(dist_matrix_subset, 
             main = paste(filtered_value), 
             clustering_distance_rows = "euclidean", 
             clustering_distance_cols = "euclidean", 
             col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255), 
             show_rownames = TRUE,
             show_colnames = TRUE,
             fontsize_row = 8,
             fontsize_col = 8,
             width = plotWidth,
             height = plotHeight,
             angle_col = 45,  
             cellwidth = 15,  
             cellheight = 10)
  })
  
  #top 50 DEGs
  hcHeatmapPlot <- reactive({
    req(global$vsd, input$intGroupSelect)
    
    comparisonGroup <- colData(global$vsd)[[input$intGroupSelect]]
    
    # Filter data based on selected subgroup
    filtered_value <- input$subGroupSelect
    selected_rows <- comparisonGroup == filtered_value
    
    # Check if there are selected rows
    if (!any(selected_rows)) {
      stop("No rows selected. Check if comparisonGroup contains the value.")
    }
    
    res <- results(global$dds)
    
    # Order results by adjusted p-value and get top 50 genes
    res_ordered <- res[order(res$padj),]
    top_genes <- head(res_ordered, 50)
    
    top_gene_names <- rownames(top_genes)
    
    # Extract expression matrix and subset for top genes
    expression_matrix <- assay(global$vsd)
    top_expression_data <- expression_matrix[top_gene_names,]
    
    expression_data_subset <- top_expression_data[, selected_rows]
    
    colors <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(255)
    
    # Check for unique breaks in color scale
    range_values <- range(expression_data_subset, na.rm = TRUE)
    if (length(unique(expression_data_subset)) <= 1) {
      stop("Data for heatmap contains insufficient unique values")
    }
    
    pheatmap(expression_data_subset, 
             main = paste("Top 50 DEGs -", filtered_value), 
             clustering_distance_rows = "euclidean", 
             clustering_distance_cols = "euclidean", 
             col = colors, 
             show_rownames = TRUE,
             show_colnames = TRUE,
             fontsize_row = 8,
             fontsize_col = 8,
             angle_col = 45,  
             cellwidth = 15,  
             cellheight = 10)
  })
  
  output$hcHeatmapPlot2 <- renderPlot({
    hcHeatmapPlot2()
  })
  
  # Render heatmap of top 50 DEGs
  output$hcHeatmapPlot <- renderPlot({
    hcHeatmapPlot()
  })
  
  observeEvent(input$save_dds, {
    req(global$dds)
    req(global$outdir)
    save_deseq_object(global$dds)
    showModal(modalDialog("DDS object saved successfully!"))
  })
  
  observeEvent(input$load_dds, {
    if (file.exists(input$load_dds$datapath)) {
      global$madeDds <- T
      global$dds <- readRDS(input$load_dds$datapath)
      global$vsd <- vst(global$dds, blind=F)
      global$results <- capture.output(results(global$dds))
      global$result_names <- resultsNames(global$dds)
      global$result_summary <- capture.output(summary(results(global$dds)))
      append_to_log("DESeqDataSet loaded successfully.")
      showModal(modalDialog("DDS object loaded successfully!"))
    } else {
      append_to_log("Could not load DDS object.")
    }

  })
  
  observe({
    req(global$dds)
    
    # Extract column names from colData(global$dds)
    comparisonGroups <- colnames(colData(global$dds))
    
    # Update the selectInput choices based on the column names
    updateSelectInput(session, "intGroupSelect", 
                      choices = comparisonGroups)
  })
  
  subgroups <- reactive({
    req(input$intGroupSelect, global$dds)
    
    selectedGroup <- input$intGroupSelect
    colData(global$dds)[[selectedGroup]]
  })
  
  observe({
    req(subgroups())
  
    uniqueSubgroups <- unique(subgroups())
    
    updateSelectInput(session, "subGroupSelect", 
                      choices = uniqueSubgroups)
  })
}
shinyApp(ui, server)

