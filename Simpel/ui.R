library(shiny)
library(shinythemes)

fluidPage(theme = shinytheme("slate"),
  sidebarLayout(
    sidebarPanel(
      fluidRow(
        column(10,
      textInput("ensemble_id_input", "Enter Ensemble id:",value = "ENSG00000100284"),
      checkboxInput("polymorphism_input", "polymorphism analysis", value=TRUE),
      checkboxInput("Mouse orthology", "Mouse orthology", value=TRUE),
      checkboxInput("Conserved_input", "Conserved?", value=TRUE),
      sliderInput("oligo_length_range", "Oligo length:",min = 0, max = 50, value = c(14,20)),
      fluidRow(
        column(10,
               h5(HTML("<b>Amount of perfect matches</b>")),
               fluidRow(
                 column(3, selectInput("dropdown_input_a", "",selected = "==", choices = c("==", "!=","<",">","<=",">="))),
                 column(4, numericInput("numeric_input_a", "", value = 1)),
                 checkboxInput("perfect_input", "Enable", value=TRUE),
               ),
               
               h5(HTML("<b>Amount of 1 mismatch</b>")),
               fluidRow(
                 column(3, selectInput("dropdown_input_b", "",selected = "<=",  choices = c("==", "!=","<",">","<=",">="))),
                 column(4, numericInput("numeric_input_b", "", value = 50)),
                 checkboxInput("mismatch_input", "Enable", value=TRUE),
               ),
               
               h5(HTML("<b>Accessibility</b>")),
               fluidRow(
                 column(3, selectInput("dropdown_input_c", "",selected = ">",  choices = c("==", "!=","<",">","<=",">="))),
                 column(4, numericInput("numeric_input_c", "", value = 1E-6)),
                 checkboxInput("Accessibility_input", "Enable", value=TRUE),
               ),
               
               h5(HTML("<b>Polymorphism frequence</b>")),
               fluidRow(
                 column(3, selectInput("dropdown_input_d", "",selected = "<",  choices = c("==", "!=","<",">","<=",">="))),
                 column(4, numericInput("numeric_input_d", "", value = 0.05)),
                 checkboxInput("Poly_input", "Enable", value=TRUE),
               ),
               
               h5(HTML("<b>Toxicity score</b>")),
               fluidRow(
                 column(3, selectInput("dropdown_input_e", "",selected = ">=",  choices = c("==", "!=","<",">","<=",">="))),
                 column(4, numericInput("numeric_input_e", "", value = 70)),
                 checkboxInput("tox_input", "Enable", value=TRUE),
               ),
               actionButton("run_button", "Run")
        ),
        column(1,
               h5(HTML("<b>Equality operator: ==</b>")),
               h5(HTML("<b>Inequality operator: !=</b>")),
               h5(HTML("<b>Less than/greater than operator: < and ></b>")),
               h5(HTML("<b>Less than or equal to/greater than or equal to operator: <= and >=</b>"))),
      )
    )
  )
    
      
      
      
      
    ),
    mainPanel(
      dataTableOutput('results1'),
      dataTableOutput('results2')
    )
  )
)