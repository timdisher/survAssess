
ui <- dashboardPage(skin = "black",
  dashboardHeader(
    shiny::tags$li(
      class = "dropdown",
      shiny::tags$style(".main-header {max-height: 150px;}"),
      shiny::tags$style(".main-header .logo {min-height: 55px;}"),
      shiny::tags$style(".sidebar-toggle {height: 20px; padding-top: 1px;}"),
      shiny::tags$style(".navbar {min-height: 50px;}")
    ),

    title = shiny::tags$img(src = "Eversana_Logo_H_RGB.png", height = 50, width = 214.3, align = "middle")
    ),
  dashboardSidebar(
   disable = TRUE ),
  dashboardBody(
    shinyalert::useShinyalert(),
    shiny::tags$style('.nav-tabs-custom .nav-tabs li.active {
    border-top-color: #ED8B00;
}'),

    shiny::tags$head(shiny::tags$style(shiny::HTML('

                            .col-sm-12 {
                            padding: 0px;
                            }

                            .col-sm-6 {
                            padding: 0px;
                            }

                            .tab-content {
                            passing: 0px;
                            }'))),

    format_edits(),
    # Boxes need to be put in a row (or column)
    fluidPage(

    fluidRow(
      column(box(title = "Data Simulation",
        solidHeader = TRUE, shiny::numericInput("n.sim",
                              min = 100,
                              max = 10000,
                              value = 150,
                              step = 25,
                              label = "Number of patients (per arm)"),
          shiny::selectInput("mat", label = "Maturity of Data",
                             choices = c("Low" = "low",
                                         "Moderate" = "mod",
                                         "Max" = "max"),
                             selected = "max"),
          hr(),
          actionButton("go", "Simulate Data"),
        height = "20em",
        width = 11

      ), width = 6),
      column(box(title = "Guess Distribution",
          solidHeader = TRUE,
          height = "20em",
          width = 11,
        shiny::selectInput("guess.dist", label = "Guess for distribution",
                           choices = c( "exp","weibull", "weibullPH", "gompertz",
                                                  "gamma", "llogis", "lnorm"),
                           selected = "Mixed"),

        shiny::selectInput("guess.ph", label = "Guess for whether PH/AFT is satisfied",
                           choices = c("Yes", "No"),
                           selected = "Yes"),

        hr(),

        shiny::actionButton("guess", label = "Make Guess")
      ),width = 6),
      tabBox(
       width = 12,
        tabPanel("KM",plotOutput("plot", height = "700px")),
        tabPanel("Cumulative Hazard Plots", plotOutput("cumhaz", height = "700px")),
       tabPanel("QQ Plots", plotOutput("qq", height = "700px")),
       tabPanel("Smoothed hazards", plotOutput("smooth", height = "700px")),
       tabPanel("Schoenfeld residuals test and plot", plotOutput("zph", height = "700px")),
       tabPanel("Parametric Fits", DT::DTOutput("fit.tab")),
       tabPanel("Parametric Fits vs KM",  plotOutput("par.km", height = "700px")),
       tabPanel("Parametric Fits vs Smoothed Hazards", plotOutput("par.smooth", height = "700px")),
       tabPanel("About",


                shiny::h2("What is Assess Survival?"),
                shiny::p("Submissions to HTA agencies often involve time to event (eg, survival) analysis. The NICE technical
                                    support document 17 provides an overview of some of the ways to assess different survival models,
                                    but trianing can be difficult given the lack of materials where the correct distrubution is known."),
                shiny::p("That’s why we developed Assess Survival – to provide an opporutunity to work through assessment of survival data to develop skills and gain an
                                    understanding of the conditions where different approaches fail."),
                shiny::h3("How does Assess Survival Work?"),
                shiny::p("First, choose the number of patients and maturity of the data. Maturity refers to the proportion of the full
                         follow-up observed. Setting this to low will highlight the challenges with lack of data and eventual extrapolation.
                         Once the data have been simulated, work through each of the plots/diagnostics to try to understand which underlying
                         function may have generated the data you are seeing. Once you feel confident, guess the distribution and whether
                         or not proportional hazards/constant acceleration factor assumptions are satisfied"),



       )


       )
    )
  )
)
)
server <- function(input, output) {
  ## app.R ##
  library(shinydashboard)

  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(flexsurv)
  library(simsurv)
  library(survminer)
  library(shinyalert)
  library(muhaz)
  library(shinyvalidate)

  lapply(list.files("./R"), function(x) source(paste0("./R/", x)))

  iv <- InputValidator$new()
  iv$add_rule("n.sim", sv_between(50, 10000))
  iv$enable()

  res <- eventReactive(input$go, {

    if (!iv$is_valid()) {
      showNotification(
        "Please fix the errors in the form before continuing",
        type = "warning"
      )
    } else {


    shinyalert::shinyalert(text = "Simulating Data", closeOnEsc = FALSE, type = "info",
                           showConfirmButton = FALSE)

    dat <- .simsurv(n = input$n.sim, maturity = input$mat)


    tests <- try(.test.ph.aft(dat$dat))
  plot <- try(.km.plots(fits = tests$all.mods,
                                                  dat = dat$dat,
                                                  include = names(tests$all.mods),
                    maxt.mature = dat$maxt.mature))


  sm.haz.p <- try(.smooth.haz(
    dat = dat$dat,
    nbreaks = 10,
    fits = tests$all.mods,
    include = names(tests$all.mods)
  ))

  if(any(class(tests) == "try-error",
         class(plot) == "try-error",
         class(sm.haz.p) == "try-error")){
    shinyalert::closeAlert()
    shinyalert(title = "Not enough data for tests. Try larger number of patients or more mature data",
               type = "error"
               )

    return()
  }

  shinyalert::closeAlert()
    dplyr::lst(dat, tests, plot, sm.haz.p)
    }
  })


  output$plot <-  renderPlot({

    res()$plot$km

  })




  output$cumhaz <- renderPlot({

    (res()$plot$loglog_gomp + labs(title = "Log Cumulative Hazard vs Time")) | (res()$plot$loglog_weib + labs(title = "Log Cumulative Hazard vs Log(Time)"))

  })


  output$qq <- renderPlot({

    res()$plot$aft_plot
  })

  output$smooth <- renderPlot({

    res()$sm.haz.p$muhaz_pw_plot
  })


  output$zph <- renderPlot({
    res()$tests$zph.plot
  })

  output$fit.tab <-
    DT::renderDataTable({res()$tests$fits %>% dplyr::mutate_if(is.numeric, round, 2)},
                        #filter = "top",
                        # # allow navigating through table with arrows
                        # extensions = 'KeyTable', options = list(keys = TRUE),
                        # save option
                        extensions = 'Buttons', options = list(
                          paging = FALSE,
                          dom = 'Bfrtip',
                          #scrollX = TRUE,
                          #scrollY = "200px",
                          initComplete = DT::JS(
                            "function(settings, json) {",
                            "$(this.api().table().header()).css({'background-color': '#00549E', 'color': '#fff'});",
                            "}")
                        )
    )



  output$par.km <- renderPlot({
    res()$plot$overlay
  })

  output$par.smooth <- renderPlot({
    res()$sm.haz.p$muhaz_ol_p
  })


  observeEvent(input$guess,{
    msg = glue::glue("You picked {input$guess.dist} | Actual dist {res()$dat$dist}\nYou picked {input$guess.ph} | Actual PH/AFT {res()$dat$ph.aft}")


    shinyalert::shinyalert(text = msg)
  })


  }






shinyApp(ui, server)
