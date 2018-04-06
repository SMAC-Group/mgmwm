#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyjs)
library(imudata)
library(mgmwm2)
library(gmwm)
library(wv)
library(simts)

data(KVH1750imu1kHzAcc)
data(KVH1750imu1kHzGyro)
data(MTIG710imu1kHz)
data(KVH1750imuAcc)
data(KVH1750imuGyro)
data(MTIG710imu50Hz)
data(ADIS16405imu100Hz)

const.RENDER_PLOT_WIDTH = 1000
const.RENDER_PLOT_HEIGHT = 600
const.RENDER_PLOT_RES = 100 # default is 72

const.FIGURE_PLOT_HEIGHT = "600px"
const.FIGURE_PLOT_HEIGHT_REDUCED = "400px"
const.FIGURE_PLOT_HEIGHT_LOGO = "100px"

const.nb_of_digits = 7

# convert degrees-per-second to radians-per-second
const.degps_2_radps = 1/360 * 2*pi

# constant default frequency for custom data
const.DEFAULT_FREQ = 1 # [Hz]
ui <- shinyUI(fluidPage(

  shinyjs::useShinyjs(),

  tags$style(HTML(".js-irs-0 .irs-single, .js-irs-0 .irs-bar-edge, .js-irs-0 .irs-bar {background: red}")),
  tags$style(HTML(".js-irs-1 .irs-single, .js-irs-1 .irs-bar-edge, .js-irs-1 .irs-bar {background: green}")),
  tags$style(type='text/css', '#summ {background-color: rgba(0,0,200,0.02); color: black; width: 500px; font-size: 14px;}'),


  title = "MGMWM GUI",
  tabsetPanel(id = "tabs",
              tabPanel("Model Data", plotOutput(outputId = "plot", height = const.FIGURE_PLOT_HEIGHT)),
              tabPanel("Selected Sensor", plotOutput(outputId = "plot2", height = const.FIGURE_PLOT_HEIGHT)),
              tabPanel("Summary", verbatimTextOutput(outputId = "summ", placeholder = FALSE)),
              tabPanel("Tutorial", htmlOutput("tuto")),
              tabPanel("Help",
                       # fluidPage("cluster"),
                       h4("Help Tab" ),
                       br(),
                       # actionButton("subClust", label = "Create Subcluster"),
                       #
                       uiOutput(outputId = "tabhelpurl"),
                       br(),br(),
                       fluidRow(
                         column(5,
                                plotOutput(outputId = "tabhelpplotlogo_pennstate", height = const.FIGURE_PLOT_HEIGHT_LOGO)
                         ),
                         column(5,
                                plotOutput(outputId = "tabhelpplotlogo_epfl", height = const.FIGURE_PLOT_HEIGHT_LOGO)
                         )
                       )
              )
  ),

  hr(),

  fluidRow(
    column(4,


             selectInput("imu_obj", "Select IMU file:",
                         c("KVH 1750 imu 1k Hz Accelerometers"="KVH1750imu1kHzAcc",
                           "KVH 1750 imu 1k Hz Gyroscopes"="KVH1750imu1kHzGyro",
                           "MTI-G-710 imu 1k Hz"="MTIG710imu1kHz",
                           "MTI-G-710 imu 50 Hz"="MTIG710imu50Hz",
                           "KVH 1750 imu 100 Hz Accelerometers"="KVH1750imuAcc",
                           "KVH 1750 imu 100 Hz Gyroscopes"="KVH1750imuGyro",
                           "ADIS 16405 imu 100Hz" = "ADIS16405imu100Hz"),
                         selected = 1),

            selectInput("sensors", "Select sensor", c("1"="1","2"="2", selected = 1)),



           actionButton("fit1", label = "Plot WV"),

           br(),

           uiOutput("choose_columns")
    ),

column(4,
       checkboxGroupInput("model", "Select Model",
                          c("Quantization Noise" = "QN",
                            "White Noise" = "WN",
                            "Random Walk" = "RW",
                            "Drift" = "DR",
                            "Auto-Regressive" = "AR"),
                          selected = "WN"),
       conditionalPanel(
         condition = "input.model.indexOf('AR')>-1",
         sliderInput("gm_nb", "Number of Gauss-Markov Processes", 1, 5, 2)
       ),

       checkboxInput("test", "Compute near-stationarity test", FALSE),

       actionButton("fit3", label = "Fit Model"), actionButton("ci", label = "Compute confidence Intervals"),

       br(),
       br(),
       br(),
       br(),

       actionButton("fit2", label = "Reduce Model Automatically")

),

column(4,

       checkboxInput("process_decomp", "Show latent processes", TRUE),
       checkboxInput("test_pval", "Paired test Selection", FALSE),

       br(),

       checkboxGroupInput("summary_plot", label = "Summary options:",
                          c("Show test result" = "test_result"),
                          selected = c("sum")),

       conditionalPanel(
         condition = "input.edit_intern == 1",
         numericInput("num", label = "Number of Simu. for Starting Values", value = 10^5),
         numericInput("seed", label = "Simulation seed", value = 2710)
       )

)
)
))

# Define server logic required to draw a histogram
server <- function(input, output, session) {

  # library or custom dataset
  v <- reactiveValues(plot = FALSE,
                      fit = FALSE,
                      gmwm = NULL,
                      form = NULL,
                      freq = 100,
                      first_gmwm = NULL,
                      n = NULL,
                      sensor_name = NULL,
                      sensor_column = NULL,
                      overlap_datasheet = FALSE,
                      y_label_with_dataunits = NA,


                      first_time_plotting_6_pack = TRUE,

                      custom_data = FALSE,
                      custom_data_name = NULL,
                      custom_data_type = NULL,
                      custom_data_size = NULL,
                      custom_data_tot_colums = NULL,
                      datasheet_noise_model = NULL,
                      datasheet_values_make_sense = FALSE)


  dsnames <- c()

  data_set <- reactive({
    inFile <- input$imu_obj

    if (is.null(inFile))
      return(kvh)

    data_set <- get(input$imu_obj)
  })

  observe({
    dsnames <- names(data_set())
    cb_options <- list()
    cb_options[ dsnames] <- dsnames
    data_set <- get(input$imu_obj)
    v$all = data_set
    updateSelectInput(session, "sensors",
                      label = "Selected sensor",
                      choices = cb_options,
                      selected = "")
  })


  # PUSHING ON BUTTON "Plot WV"
  observeEvent(input$fit1, {

    withProgress(message = 'Calculating empirical WV...', value = 0, {

      v$plot = TRUE
      v$fit = FALSE

      my_data = get(input$imu_obj)
      Xt = my_data[input$sensors][[1]]

      v$sensor_name = input$imu_obj
      v$sensor_column = input$sensors
      v$freq = attr(my_data, 'freq')
      v$custom_data = FALSE
      if (input$sensors == "Gyro. X" || input$sensors == "Gyro. Y" || input$sensors == "Gyro. Z"){
        v$y_label_with_dataunits = expression(paste("Wavelet Variance ", nu, " [", rad^2/s^2, "]"))
      } else if (input$sensors == "Acc.X" || input$sensors == "Acc.Y" || input$sensors == "Acc.Z"){
        v$y_label_with_dataunits = expression(paste("Wavelet Variance ", nu, " [", m^2/s^4, "]"))
      }


      v$form = Xt

      updateNavbarPage(session, "tabs", selected = "Selected Sensor")
    })
  })


  observeEvent(input$fit3, {

    withProgress(message = 'Fitting desired model...', value = 0, {

      if (is.null(v$first_gmwm)){
        v$first_gmwm = TRUE
      }
      v$fit = TRUE
      v$plot = FALSE

      my_data = get(input$imu_obj)
      Xt = my_data[input$sensors][[1]]

      first = TRUE
      counter_model_size = 0

      if ("AR" %in% input$model){
        for (i in 1:input$gm_nb){
          counter_model_size = counter_model_size + 1
          if (first == TRUE){
            model = AR1()
            first = FALSE
          }else{
            model = model + AR1()
          }
        }
      }

      if ("WN" %in% input$model){
        counter_model_size = counter_model_size + 1
        if (first == TRUE){
          model = WN()
          first = FALSE
        }else{
          model = model + WN()
        }
      }

      if ("QN" %in% input$model){
        counter_model_size = counter_model_size + 1
        if (first == TRUE){
          model = QN()
          first = FALSE
        }else{
          model = model + QN()
        }
      }


      if ("RW" %in% input$model){
        counter_model_size = counter_model_size + 1
        if (first == TRUE){
          model = RW()
          first = FALSE
        }else{
          model = model + RW()
        }
      }

      if ("DR" %in% input$model){
        counter_model_size = counter_model_size + 1
        if (first == TRUE){
          model = DR()
          first = FALSE
        }else{
          model = model + DR()
        }
      }

      if (is.null(model)){
        model = 3*AR1()
      }

      if (is.null(input$seed)){
        input$seed = 1982
      }

      if (is.null(input$num)){
        input$num = 10^5
      }
      v$gmwm = mgmwm(mimu = Xt, model = model, CI = FALSE, stationarity_test = FALSE, B_stationarity_test = NULL,
                     alpha_ci = NULL, alpha_near_test = NULL, seed = 2710, n_boot_ci_max = 300)
      v$form = v$gmwm
      v$first_gmwm = FALSE

      updateNavbarPage(session, "tabs", selected = "Selected Sensor")

    })

  })

  # BUTTON REDUCE MODEL WHICH WILL USE THE model_selection FUNCTION
  observeEvent(input$fit2, {

    withProgress(message = 'Reducing model automatically...', value = 0, {

      if (is.null(v$first_gmwm)){
        v$first_gmwm = TRUE
      }
      v$fit = TRUE
      v$plot = FALSE

      my_data = get(input$imu_obj)
      Xt = my_data[input$sensors][[1]]

      first = TRUE
      counter_model_size = 0

      if ("AR" %in% input$model){
        for (i in 1:input$gm_nb){
          counter_model_size = counter_model_size + 1
          if (first == TRUE){
            model = AR1()
            first = FALSE
          }else{
            model = model + AR1()
          }
        }
      }

      if ("WN" %in% input$model){
        counter_model_size = counter_model_size + 1
        if (first == TRUE){
          model = WN()
          first = FALSE
        }else{
          model = model + WN()
        }
      }

      if ("QN" %in% input$model){
        counter_model_size = counter_model_size + 1
        if (first == TRUE){
          model = QN()
          first = FALSE
        }else{
          model = model + QN()
        }
      }


      if ("RW" %in% input$model){
        counter_model_size = counter_model_size + 1
        if (first == TRUE){
          model = RW()
          first = FALSE
        }else{
          model = model + RW()
        }
      }

      if ("DR" %in% input$model){
        counter_model_size = counter_model_size + 1
        if (first == TRUE){
          model = DR()
          first = FALSE
        }else{
          model = model + DR()
        }
      }

      if (is.null(model)){
        model = 3*AR1()
      }

      a = model_selection(mimu, model, s_est = NULL,
                          alpha_paired_test = NULL, seed = 2710)
      v$form = a



      updateNavbarPage(session, "tabs", selected = "Selected Sensor")
    })
  })


  output$plot2 <- renderPlot({
    if (class(v$form) == "mgmwm"){
      plot(v$form, decomp = input$process_decomp)
    }else{
      plot(v$form)
    }
  })

  output$plot <- renderPlot({
    N = length(v$all)

    if (N > 3){
      par(mfrow = c(2,3))
    }else{
     par(mfrow = c(1,3))
    }


    for (i in 1:N){
      if (i == 1){
        plot(v$all[[i]])
      }else{
        plot(v$all[[i]], legend = FALSE)
      }

    }
  })

  output$summ <- renderPrint({
    if (v$fit){

      summmary_of_gmwm = list(v$form$estimate, v$form$obj.value, v$form$p_value, v$form$test_res)

      if("ci" %in% input$summary_plot){
        # summmary_of_gmwm
        cat("Objective Function: ", v$form$obj.value, "\n")
        cat("P-value: ", v$form$test_res, "\n\n")
        cat("Test Result: ", v$form$test.result, "\n\n\n")
        v$form$estimate
      } else {
        cat("Objective Function: ", v$form$obj.value, "\n")
        v$form$obj.value
        v$form$estimate
      }
    }
  })


  output$tuto <- renderUI({
    tags$iframe(src = "https://www.youtube.com/embed/HPPj6viIBmU", height=400, width=600)
  })
}



# Run the application
shinyApp(ui = ui, server = server)

