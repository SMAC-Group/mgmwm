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
library(mgmwm)
library(wv)
library(gmwm)
library(simts)
library(iterpc)

data(mtig1khrz)
data(adis16405)
data(kvh1750.gyro)
data(kvh1750.acc)


n1 = 10000
n2 = 50000
n3 = 100000
n4 = 500000

model1 = AR1(.995, sigma2 = 4e-6) + WN(.1) + RW (1e-8)

Wt =  gen_gts(n3, model1)
Xt =  gen_gts(n1, model1)
Yt =  gen_gts(n2, model1)
Zt =  gen_gts(n4, model1)



test_data1 = make_mimu(Wt ,Xt, Yt, Zt, freq = 100, unit = "s",
                       sensor.name = "Testing data", exp.name = c("K1", "K2", "K3", "K4"))

Wt =  gen_gts(n3, model1)
Xt =  gen_gts(n1, model1)
Yt =  gen_gts(n2, model1)
Zt =  gen_gts(n4, model1)

test_data2 = make_mimu(Wt ,Xt, Yt, Zt, freq = 100, unit = "s",
                       sensor.name = "Testing data", exp.name = c("K1", "K2", "K3", "K4"))

Wt =  gen_gts(n3, model1)
Xt =  gen_gts(n1, model1)
Yt =  gen_gts(n2, model1)
Zt =  gen_gts(n4, model1)

test_data3 = make_mimu(Wt ,Xt, Yt, Zt, freq = 100, unit = "s",
                       sensor.name = "Testing data", exp.name = c("K1", "K2", "K3", "K4"))

Wt =  gen_gts(n3, model1)
Xt =  gen_gts(n1, model1)
Yt =  gen_gts(n2, model1)
Zt =  gen_gts(n4, model1)

test_data4 = make_mimu(Wt ,Xt, Yt, Zt, freq = 100, unit = "s",
                       sensor.name = "Testing data", exp.name = c("K1", "K2", "K3", "K4"))

Wt =  gen_gts(n3, model1)
Xt =  gen_gts(n1, model1)
Yt =  gen_gts(n2, model1)
Zt =  gen_gts(n4, model1)

test_data5 = make_mimu(Wt ,Xt, Yt, Zt, freq = 100, unit = "s",
                       sensor.name = "Testing data", exp.name = c("K1", "K2", "K3", "K4"))

Wt =  gen_gts(n3, model1)
Xt =  gen_gts(n1, model1)
Yt =  gen_gts(n2, model1)
Zt =  gen_gts(n4, model1)

test_data6 = make_mimu(Wt ,Xt, Yt, Zt, freq = 100, unit = "s",
                       sensor.name = "Testing data", exp.name = c("K1", "K2", "K3", "K4"))


test_data = list(test_data1,test_data2,test_data3,test_data4,test_data5,test_data6)
names(test_data) = c("Gyro. X", "Gyro. Y","Gyro. Z", "Acc. X", "Acc. Y","Acc. Z")


const.RENDER_PLOT_WIDTH = 1000
const.RENDER_PLOT_HEIGHT = 800
const.RENDER_PLOT_RES = 100 # default is 72

const.FIGURE_PLOT_HEIGHT = "750px"
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
              tabPanel("All models", plotOutput(outputId = "plot3", height = const.FIGURE_PLOT_HEIGHT)),
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
                       c("ADIS 16405 imu 100Hz" = "adis16405",
                         "MTI-G-710 imu 1k Hz" = "mtig1khrz",
                         "KVH 1750 Accelerometer" = "kvh1750.acc",
                         "KVH 1750 Gyroscopes" = "kvh1750.gyro",
                         "Test data set" = "test_data"),
                       selected = 1),

           selectInput("sensors", "Select sensor", c("1"="1","2"="2", selected = 1)),



           actionButton("fit1", label = "Plot WV"),


           br(),

           uiOutput("choose_columns")

    ),
    column(4,
           checkboxGroupInput("model", "Select Latent Processes",
                              c("Quantization Noise" = "QN",
                                "White Noise" = "WN",
                                "Random Walk" = "RW",
                                "Drift" = "DR",
                                "Auto-Regressive" = "AR"),
                              selected = "WN"),
           conditionalPanel(
             condition = "input.model.indexOf('AR')>-1",
             sliderInput("gm_nb", "Number of Auto-Regressive Processes", 1, 5, 2)
           ),
           br(),
           checkboxInput("ci", "Compute confidence Intervals", FALSE),
           checkboxInput("test", "Compute near-stationarity test", FALSE),


           actionButton("fit3", label = "Fit Model")
           #checkboxInput("plot_ci", "Plot CI", FALSE)


    ),
    column(4,



           br(),
           br(),

           actionButton("fit2", label = "Reduce Model Automatically"),

           br(),
           br(),
           selectInput("sel_mod", "Select Estimated Model", c("1"="1","2"="2", selected = 1)),

           checkboxInput("rm", "Equivalent models", FALSE),
           checkboxInput("eq", "Plot equivalent models", FALSE),
           checkboxInput("compare", "Compare models", FALSE),
           checkboxInput("zoom", "Compare equivalent models", FALSE)


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
                      datasheet_values_make_sense = FALSE,
                      model_selection = FALSE,
                      model_selection_model_selected_name = NULL,
                      model_selection_model_selected = NULL,
                      user_selected_model = NULL)


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

  get_model = reactive({
    #print("22")
    #if (v$model_selection == TRUE){
    #  print(v$model_selection)


    #  print("3")
    #  print(input$sel_mod != "" && (is.null(v$model_selection_model_selected_name) || input$mod_sel != v$model_selection_model_selected_name))
    #  print(v$model_selection_model_selected_name)
    #  print(input$mod_sel)
    #  print(class(input$mod_sel))

    #  if (!is.null(input$mod_sel)){
    #  v$user_selected_model = input$mod_sel
    #    print("1111")
    #  }

    #print(v$user_selected_model)
    #if (v$user_selected_model != "" && (is.null(v$model_selection_model_selected_name) || v$user_selected_model != v$model_selection_model_selected_name)){
    #  updateNavbarPage(session, "tabs", selected = "Selected Model")
    #  print("4")
    #  v$model_selection_model_selected_name = v$form$model_name[which(v$form$model_name %in% input$sel_mod)]
    #  print("44")
    #  print(v$model_selection_model_selected_name)
    #  v$new_model_selection = TRUE
    v$form$model_nested[[which(v$form$model_name %in% input$sel_mod)]]
    #}
    #}
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


      v$gmwm = mgmwm(mimu = Xt, model = model, CI = input$ci, stationarity_test = input$test)
      v$form = v$gmwm
      v$first_gmwm = FALSE
      v$model_selection = FALSE
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

      a = model_selection(mimu = Xt, model = model, s_est = NULL,
                          alpha_ci_cvwvic = NULL,b_ci_wvic = 500,  seed = 2710)

      v$form = a
      v$model_selection = TRUE
      updateNavbarPage(session, "tabs", selected = "Selected Sensor")
      updateSelectInput(session, "sel_mod",
                        label = "Select Estimated Model",
                        choices = v$form$model_name,
                        selected = "")
    })
  })


  observeEvent(input$rm, {
    if (!is.null(v$form$model_name)){
      if (input$rm){
        updateSelectInput(session, "sel_mod",
                          label = "Select Estimated Model",
                          choices = v$form$model_name[v$form$selection_decision != "Model not appropriate"],
                          selected = "")
      }else{
        updateSelectInput(session, "sel_mod",
                          label = "Select Estimated Model",
                          choices = v$form$model_name,
                          selected = "")
      }

    }
  })

  output$model_selection <- renderUI({
    if(class(v$form == cvwvic)){
      name = v$form$model_name
    }
  })


  output$plot2 <- renderPlot({
    if (class(v$form) == "mgmwm"){
      # if(input$plot_ci == TRUE){
      #   plot_CI(v$form)
      # }else{
        plot.mgmwm(v$form, decomp = TRUE)
      #}
    }else if (class(v$form) == "cvwvic"){
      plot(v$form)
    }else{
      plot(v$form)
    }
  })

  output$plot3 <- renderPlot({
    if (class(v$form) == "cvwvic"){
      if (input$eq){
        plot(v$form, type = "equivalent")
      }else if (input$compare){
        plot(v$form, type = "wvic_all")
      }else if(input$zoom){
        plot(v$form, type = "wvic_equivalent")
      }else{
        if (input$sel_mod != ""){
          plot(v$form, decomp = TRUE,  model = get_model())
        }else{
          plot(v$form, decomp = TRUE)
        }
      }
    }else{
      plot(v$form)
    }
  })


  output$plot <- renderPlot({
    N = length(v$all)

    if (N > 3){
      par(mfrow = c(2,3), oma = c(2,2,0,0) + 0.1,
          mar = c(2.5,2.5,1,1) + 0.1)
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
      if(class(v$form) == "mgmwm"){
        summary.mgmwm(v$form)

      } else {
        summary.cvwvic(v$form)
      }
    }
  })

  output$tuto <- renderUI({
    tags$iframe(src = "https://www.youtube.com/embed/Rj-P-vu_7aI", height=400, width=600)
  })
}



# Run the application
shinyApp(ui = ui, server = server)

