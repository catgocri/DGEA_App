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
library(DT)

# Fetch data from the pipeline
get_pipeline_data <- function() {
  samplesheet <- read_tsv("samplesheet_P28103.tsv")
  counts <- read_delim("P28103_expression_data_2.csv", delim=";") %>% data.frame()
  list(samplesheet = samplesheet, counts = counts)
}

ui <- fluidPage(
  titlePanel("Pipeline Starter & DGE Analysis"),
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
      actionButton("run_pipeline", "Run nf-core/RNA-Seq Pipeline")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("MA Plot", plotOutput("maPlot")),
        tabPanel("PCA Plot", plotOutput("pcaPlot")),
        tabPanel("Cook's Distance Plot", plotOutput("cooksPlot")),
        tabPanel("MDS Plot", plotOutput("mdsPlot")),
        tabPanel("Volcano Plot", plotOutput("volcanoPlot")),
        tabPanel("Heatmap Plot", plotOutput("heatmapPlot")),
        tabPanel("HC Heatmap Plot", plotOutput("hcHeatmapPlot")),
        tabPanel("Analysis Summary", verbatimTextOutput("summary")),
        tabPanel("Log Messages", div(style="overflow-y:scroll;", verbatimTextOutput("logging")))
      ),
      tabsetPanel(
        tabPanel("sheetMaker.sh CSV Preview", tableOutput("sheetMakerCsv")),
        tabPanel("sheetMaker2.sh CSV Preview", tableOutput("sheetMaker2Csv")),
        tabPanel("Samplesheet TSV Preview", tableOutput("samplesheetTsv")),
        tabPanel("MultiQC Finder Pre-Alignment", uiOutput("multiqcViewerPre")),
        tabPanel("MultiQC Finder Post-Alignment", uiOutput("multiqcViewer"))
      ),
    )
  )
)

server <- function(input, output, session) {
  
  # Max input size 30MB
  options(shiny.maxRequestSize=30*1024^2)
  
  ### START OF PIPELINE INPUT SELECTION ###
  # Define starting directories for both directory inputs
  shinyDirChoose(input, 'samples', roots = c(root = '/'))
  shinyDirChoose(input, 'outdir', roots = c(root = '/'))
  
  # Where to store the directories
  
  global <- reactiveValues(samples_dir = NULL, outdir = NULL, sheetMaker_file_dir = NULL, sheetMaker2_file_dir = NULL, nfcore_command_text = NULL)
  
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
    output$outdir <- renderText({ global$outdir })
  })
  ### END OF PIPELINE INPUT SELECTION ###
  
  
  
  
  
  
  
  ### START PIPELINE CODE AND LOGGING FUNCTIONALITY ###
  
  log_messages <- reactiveValues(value = "")
  
  temp_log_file <- reactiveVal(NULL)
  
  append_to_log <- function(input) {
    log_messages$value <- paste(log_messages$value, paste(input, collapse = "\n"), sep = "\n")
  }
  
  run_command <- function(input, temp_file) {
    system(paste(input, "1>", shQuote(temp_file), "2>&1"))
  }
  
  observeEvent(input$load_sheets, {
    if (!is.null(global$samples_dir) && !is.null(global$outdir)) {
      # Generate command to execute
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
  
}
shinyApp(ui, server)

