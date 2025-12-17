library(DT)
library(shiny)
library(shinythemes)

ui <- fluidPage(
  theme = shinythemes::shinytheme("flatly"),
  tags$head(tags$style(
    HTML(
      "
      #shiny-notification-panel {
        position: fixed;
        top: calc(50%);
        left: calc(50%);
        width: 500px;
        height: 200px;
        transform: translate(-50%, -50%);
        z-index: 9999;
      }
    "
    )
  )),
  sidebarLayout(
    sidebarPanel(fluidRow(
      column(
        12,
        textInput("ensemble_id_input", "Enter Ensemble id:", value = "ENSG00000100284"),
        checkboxInput("polymorphism_input", "Polymorphism analysis", value =
                        TRUE),
        checkboxInput("Conserved_input", "Conserved & Orthology", value =
                        TRUE),
        checkboxInput("linux_input", "Running on Linux-OS", value = TRUE),
        sliderInput(
          "oligo_length_range",
          "Oligo length:",
          min = 0,
          max = 50,
          value = c(14, 20)
        ),
        fluidRow(
          column(
            12,
            h5(HTML("<b>Amount of perfect matches</b>")),
            fluidRow(
              column(
                4,
                selectInput(
                  "dropdown_input_a",
                  "",
                  selected = "==",
                  choices = c("==", "!=", "<", ">", "<=", ">=")
                )
              ),
              column(5, numericInput("numeric_input_a", "", value = 1)),
              checkboxInput("perfect_input", "Enable", value = TRUE),
            ),
            
            h5(HTML("<b>Amount of 1 mismatch</b>")),
            fluidRow(
              column(
                4,
                selectInput(
                  "dropdown_input_b",
                  "",
                  selected = "<=",
                  choices = c("==", "!=", "<", ">", "<=", ">=")
                )
              ),
              column(5, numericInput("numeric_input_b", "", value = 50)),
              checkboxInput("mismatch_input", "Enable", value = TRUE),
            ),
            
            h5(HTML("<b>Accessibility</b>")),
            fluidRow(
              column(
                4,
                selectInput(
                  "dropdown_input_c",
                  "",
                  selected = ">",
                  choices = c("==", "!=", "<", ">", "<=", ">=")
                )
              ),
              column(5, numericInput("numeric_input_c", "", value = 1E-6)),
              checkboxInput("Accessibility_input", "Enable", value =
                              TRUE),
            ),
            
            h5(HTML("<b>Polymorphism frequence</b>")),
            fluidRow(
              column(
                4,
                selectInput(
                  "dropdown_input_d",
                  "",
                  selected = "<",
                  choices = c("==", "!=", "<", ">", "<=", ">=")
                )
              ),
              column(5, numericInput("numeric_input_d", "", value = 0.05)),
              checkboxInput("Poly_input", "Enable", value = TRUE),
            ),
            
            h5(HTML("<b>Toxicity score</b>")),
            fluidRow(
              column(
                4,
                selectInput(
                  "dropdown_input_e",
                  "",
                  selected = ">=",
                  choices = c("==", "!=", "<", ">", "<=", ">=")
                )
              ),
              column(5, numericInput("numeric_input_e", "", value = 70)),
              checkboxInput("tox_input", "Enable", value = TRUE),
            ),
            actionButton("run_button", "Run")
          )
        )
      )
    ), width = 3),
    
    # Updated main panel for RNase H script.
    mainPanel(tabsetPanel(
      id = "tabs_main",
      tabPanel(
        "Sequence results",
        downloadButton("Download_input", "Download Results", style = "margin-top: 20px"),
        hr(),
        DT::dataTableOutput('results1'),
        DT::dataTableOutput('results2')
      ),
      
      tabPanel(
        "RNase H cleavage results",
        h3(textOutput("rnaseh_title")),
        div(uiOutput("rnaseh_info"), style = "margin-bottom: 15px;"),
        hr(),
        
        fluidRow(column(
          6,
          numericInput(
            "mod_5prime",
            label = tagList(
              "Amount of modified nucleotides at the 5' end  ",
              tags$span(
                tags$img(src = "questionmark.png", height = "20px"),
                title = "Test on hover",
                `data-toggle` = "tooltip",
                style = "cursor: pointer;"
              )
            ),
            value = 0,
            min = 0,
            max = 10,
            width = "60%"
          ),
        ), column(
          6,
          numericInput(
            "mod_3prime",
            label = tagList(
              "Amount of modified nucleotides at the 3' end  ",
              tags$span(
                tags$img(src = "questionmark.png", height = "20px"),
                title = "Test on hover",
                `data-toggle` = "tooltip",
                style = "cursor: pointer;"
              )
            ),
            value = 0,
            min = 0,
            max = 10,
            width = "60%"
          )
        )),
        actionButton("add_mods", "Apply end modifications", class = "btn-primary"),
        hr(),
        
        downloadButton("download_rnaseh", "Download results", style = "margin-bottom: 15px;"),
        dataTableOutput("rnaseh_results"),
        hr(),
        
        h3("Visualised cleavage site: "),
        uiOutput("cleavage_visual"),
      ),
      
      tabPanel(
        "Off target results",
        h3(textOutput("offtarget_title")),
        textOutput("numb_offtargets"),
        textOutput("aso_seq"),
        hr(),
        selectInput(
          "user_mismatch",
          label = tagList(
            "Select number of mismatches allowed  ",
            tags$span(
              tags$img(src = "questionmark.png", height = "20px"),
              title = "Test on hover",
              `data-toggle` = "tooltip",
              style = "cursor: pointer;"
            )
          ),
          choices = list(
            "0" = 0,
            "1" = 1,
            "2" = 2,
            "3" = 3
          )
        ),
        actionButton("apply_mismatch", "Apply", class = "btn-primary"),
        hr(),
        downloadButton("download_offtarget", "Download results", style = "margin-bottom: 15px;"),
        DTOutput("offtarget_results")
      )
    ),
    tabPanel(
      "Off target results", 
      h3(textOutput("offtarget_title")),
      textOutput("aso_seq"),
      textOutput("numb_offtargets"),
      hr(),
      fluidRow(
        column(6, 
               selectInput("user_mismatch", "Select number of mismatches allowed", choices = list("0" = 0, "1" = 1, "2" = 2, "3" = 3)),
               actionButton("apply_mismatch", "Apply")
               ),
        column(6,
               fluidRow("Run off-target tissue expression and OMIM disease search (may take some time)"),
               
               fluidRow(
                 selectInput("target_tissue", "Select target tissue", choices = c("Brain",
                                                                                  "Eye",
                                                                                  "Endrocrine tissue",
                                                                                  "Respiratory system",
                                                                                  "Proximal digestive tract",
                                                                                  "Gastrointestinal tract",
                                                                                  "Liver & galbladder",
                                                                                  "Pancreas",
                                                                                  "Kidney & urinary bladder",
                                                                                  "Male tissues",
                                                                                  "Female tissues",
                                                                                  "Muscle tissues",
                                                                                  "Connective & soft tissues",
                                                                                  "Skin",
                                                                                  "Bone marrow & lymphoid tissues")),
                 actionButton("PAtlas_OMIM_search", "Run"))
               )
      ),
      hr(),
      downloadButton("download_offtarget", "Download results", style = "margin-bottom: 15px;"),
      DTOutput("offtarget_results"),
    ),
    width = 9)
  ),
  tags$script(
    HTML(
      '$(function () { $("[data-toggle=\'tooltip\']").tooltip(); });'
    )
  )
)