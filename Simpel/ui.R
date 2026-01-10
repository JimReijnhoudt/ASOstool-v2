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

      .navbar-brand {
        display: flex !important;
        align-items: center !important;
        height: 60px !important;
      }
    "
    )
  )),
  tags$nav(
    class = "navbar navbar-default",
    style = "margin-bottom: 5px;",
    tags$div(class = "container-fluid", tags$div(
      class = "navbar-header",
      tags$span(class = "navbar-brand", style = "font-weight: bold; font-size: 30px;", HTML("ER<u>ASO</u>R"))
    ))
  ),
  sidebarLayout(
    sidebarPanel(fluidRow(
      column(
        12,
        textInput(
          "ensemble_id_input",
          label = tagList(
            "Enter Ensembl id:  ",
            tags$span(
              tags$img(
                src = "questionmark.png",
                height = "20px",
                style = "margin-bottom: 3px;"
              ),
              title = "Please enter a valid Ensembl Gene ID (ENSG…) for the gene of interest.",
              `data-placement` = "right",
              `data-toggle` = "tooltip",
              style = "cursor: pointer;"
            )
          ),
          value = "ENSG00000100284"
        ),
        checkboxInput(
          "polymorphism_input",
          label = tagList(
            "Polymorphism analysis ",
            tags$span(
              tags$img(
                src = "questionmark.png",
                height = "20px",
                style = "margin-bottom: 3px;"
              ),
              title = "This feature identifies single nucleotide polymorphisms (SNPs) for the selected Ensembl gene. SNPs are genetic variations in which one nucleotide has changed, for example, a C to a T.",
              `data-placement` = "right",
              `data-toggle` = "tooltip",
              style = "cursor: pointer;"
            )
          ),
          value = TRUE
        ),
        checkboxInput(
          "ASO_ending_G",
          label = tagList(
            "Filter ASO ending in G ",
            tags$span(
              tags$img(
                src = "questionmark.png",
                height = "20px",
                style = "margin-bottom: 3px;"
              ),
              title = "Filters ASO sequences ending with G. These sequences could lead to unwanted toxic effects. ",
              `data-placement` = "right",
              `data-toggle` = "tooltip",
              style = "cursor: pointer;"
            )
          ),
          value = FALSE
        ),
        checkboxInput(
          "Conserved_input",
          label = tagList(
            "Conserved & Orthology ",
            tags$span(
              tags$img(
                src = "questionmark.png",
                height = "20px",
                style = "margin-bottom: 3px;"
              ),
              title = "This feature compares an Ensembl gene to the mouse genome to identify conserved regions and orthologous genes. Because of the strong genetic similarity between mice and humans, this helps determine whether the gene can be studied in vivo in mice.",
              `data-placement` = "right",
              `data-toggle` = "tooltip",
              style = "cursor: pointer;"
            )
          ),
          value = FALSE
        ),
        checkboxInput(
          "linux_input",
          label = tagList(
            "Running on Linux-OS ",
            tags$span(
              tags$img(
                src = "questionmark.png",
                height = "20px",
                style = "margin-bottom: 3px;"
              ),
              title = "Enable this setting only when running on a Linux operating system. Some features, such as the ViennaRNA analysis, require Linux. Enabling this option on other systems, such as Windows or Mac, may cause the program to fail.",
              `data-placement` = "right",
              `data-toggle` = "tooltip",
              style = "cursor: pointer;"
            )
          ),
          value = TRUE
        ),
        sliderInput(
          "oligo_length_range",
          label = tagList(
            "Oligo length: ",
            tags$span(
              tags$img(
                src = "questionmark.png",
                height = "20px",
                style = "margin-bottom: 3px;"
              ),
              title = "This option allows the user to input a length/range for the target mRNA sequences and complementary ASO sequences. Larger lengths/ranges will generally have a shorter runtime, while shorter lengths/ranges will take longer to process. ",
              `data-placement` = "right",
              `data-toggle` = "tooltip",
              style = "cursor: pointer;"
            )
          ),
          min = 0,
          max = 50,
          value = c(19, 20)
        ),
        fluidRow(
          column(
            12,
            
            h5(tagList(
              HTML("<b>Amount of perfect matches</b> "),
              tags$span(
                tags$img(
                  src = "questionmark.png",
                  height = "20px",
                  style = "margin-bottom: 3px;"
                ),
                title = "This option sets the amount of perfect matches allowed. The value is used as a quality control threshold of the target mRNA sequences. ",
                `data-toggle` = "tooltip",
                `data-placement` = "right",
                style = "cursor: pointer;"
              )
            )),
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
              column(5, numericInput("numeric_input_a", "", value = 5)),
              checkboxInput("perfect_input", "Enable", value = TRUE)
            ),
            
            h5(tagList(
              HTML("<b>Amount of 1 mismatch </b>"),
              tags$span(
                tags$img(
                  src = "questionmark.png",
                  height = "20px",
                  style = "margin-bottom: 3px;"
                ),
                title = "This option sets the amount of mismatches allowed. The value is used as a quality control threshold of the target mRNA sequences.",
                `data-toggle` = "tooltip",
                `data-placement` = "right",
                style = "cursor: pointer;"
              )
            )),
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
              column(5, numericInput("numeric_input_b", "", value = 100)),
              checkboxInput("mismatch_input", "Enable", value = TRUE)
            ),
            
            h5(tagList(
              HTML("<b>Accessibility </b>"),
              tags$span(
                tags$img(
                  src = "questionmark.png",
                  height = "20px",
                  style = "margin-bottom: 3px;"
                ),
                title = "This option is a quality control score for the accessibility of the target mRNA sequences. The score estimates how accessible the target mRNA sequences are, which is an important factor in finding potentially effective ASOs. ",
                `data-toggle` = "tooltip",
                `data-placement` = "right",
                style = "cursor: pointer;"
              )
            )),
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
                              TRUE)
            ),
            
            h5(tagList(
              HTML("<b>Polymorphism frequence </b>"),
              tags$span(
                tags$img(
                  src = "questionmark.png",
                  height = "20px",
                  style = "margin-bottom: 3px;"
                ),
                title = "This quality control score estimates the probability that a polymorphism (SNP) is present in the target mRNA sequences. Lower values are preferred, as they indicate a reduced likelihood of polymorphic variation.",
                `data-toggle` = "tooltip",
                `data-placement` = "right",
                style = "cursor: pointer;"
              )
            )),
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
              checkboxInput("Poly_input", "Enable", value = TRUE)
            ),
            
            h5(tagList(
              HTML("<b>Toxicity score </b>"),
              tags$span(
                tags$img(
                  src = "questionmark.png",
                  height = "20px",
                  style = "margin-bottom: 3px;"
                ),
                title = "This quality control score estimates the potential toxicity of the complementary ASO targeting the mRNA sequence, providing insight into the safety risk of a chosen ASO. Higher values correspond to lower toxicity, which is favourable.",
                `data-toggle` = "tooltip",
                `data-placement` = "right",
                style = "cursor: pointer;"
              )
            )),
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
              checkboxInput("tox_input", "Enable", value = TRUE)
            ),
            actionButton("run_button", "Run")
          )
        )
      )
    ),width = 3),
    
    # Updated main panel for RNase H script.
    mainPanel(tabsetPanel(
      id = "tabs_main",
      tabPanel(
        "Sequence results",
        div(
          "This is the main page of the application, which displays the primary output table. The table lists all potential target mRNA sequences for the provided Ensembl ID and includes information to help the user select suitable ASO targets. After selecting an ASO, the application automatically redirects to the RNase H tab, and the chosen ASO becomes available for further analysis on all other tabs.",
          style = "margin-top: 20px; font-size: 18px;"
        ),
        hr(),
        fluidRow(
          column(6,
                 p("This table shows the number of ASOs that didn't meet filtering criterea."),
                 DTOutput("unfiltered_results_table")
                ),
          column(2,
                 p("Download unfiltered results without off-target search here: "),
                 downloadButton("Download_unfiltered", "Download Unfiltered Results")
               )),
        hr(),
        DT::dataTableOutput('results1'),
        hr(),
        downloadButton("Download_filtered", "Download Filtered Results"),
      ),
      
      tabPanel(
        "RNase H cleavage results",
        div(
          "This page predicts the optimal binding site for RNase H on the chosen target mRNA sequence. It uses a matrix of dinucleotide values and a sliding window approach to evaluate all possible binding sites along the mRNA sequence. The script accounts for chemical modifications on the ASO’s 5’ and 3’ ends, which can affect possible binding sites. Each potential site is scored based on the average from the dinucleotide matrix, and the results are ranked to highlight the most promising locations for RNase H cleavage.",
          style = "margin-top: 20px; font-size: 18px;"
        ),
        hr(),
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
                title = "This setting adds modifications to the 5' end of the ASO sequence, which affect the accessible binding sites on the target mRNA. An overlap of up to 2 modified nucleotides on the 3' end of the target mRNA is allowed.",
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
                title = "This setting adds modifications to the 3' end of the ASO sequence, which affect the accessible binding sites on the target mRNA. An overlap of up to 4 modified nucleotides on the 5' end of the target mRNA is allowed.",
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
        div("This page displays the off-targets for the selected target mRNA sequence based on the allowed number of mismatches. For the complementary ASO sequence, it shows mismatches, deletions, insertions, protein name and additional information from GGGenome. The protein name is then given to Protein Atlas to retrieve the tissue expression data. This page also calculates the accessibility potential of each off-target.", style="margin-top: 20px; font-size: 18px;"),
        hr(),
        uiOutput("gggenome_status"),
        h3(textOutput("offtarget_title")),
        textOutput("aso_seq"),
        textOutput("numb_offtargets"),
        hr(),
        fluidRow(
          column(6, 
                 selectInput(
                   "user_mismatch",
                   label = tagList(
                     "Select number of mismatches allowed  ",
                     tags$span(
                       tags$img(src = "questionmark.png", height = "20px"),
                       title = "This setting identifies off-target sequences according to the number of mismatches permitted in the target mRNA, with a default value of one mismatch. Allowing more mismatches results in more off-targets being detected but also increases the processing time.",
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
          ),
          column(6,
                 fluidRow("Run off-target tissue expression and OMIM disease search (may take some time)"),
                 
                 fluidRow(
                   selectInput("target_tissue", label = tagList(
                     "Select target tissue  ",
                     tags$span(
                       tags$img(src = "questionmark.png", height = "20px"),
                       title = "For each off-target, tissue-specific expression data is retrieved from the Protein Atlas. This data, combined with the off-target information, is then used to query OMIM for associated diseases. Currently, the analysis supports the brain, eyes, and liver, with additional tissues planned as further developments become available.",
                       `data-toggle` = "tooltip",
                       style = "cursor: pointer;"
                     )
                   ),
                   choices = c(
                     "Brain",
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
                     "Bone marrow & lymphoid tissues"
                   )
                   ),
                   actionButton("PAtlas_OMIM_search", "Run"))
          )
        ),
        hr(),
        downloadButton("download_offtarget", "Download results", style = "margin-bottom: 15px;"),
        DTOutput("offtarget_results"),
      )
    ),
    width = 9
  )
  ),
  tags$script(
    HTML(
      '$(function () { $("[data-toggle=\'tooltip\']").tooltip(); });'
    )
  )
)