
ui <- dashboardPage(
  dashboardHeader(title = "Assess Survival"),
  dashboardSidebar(
   disable = TRUE ),
  dashboardBody(
    shinyalert::useShinyalert(),
    # Boxes need to be put in a row (or column)
    fluidPage(
    fluidRow(
      box(shiny::numericInput("n.sim",
                              min = 100,
                              max = 10000,
                              value = 150,
                              step = 25,
                              label = "Number of patients"),
          shiny::selectInput("mat", label = "Maturity of Data",
                             choices = c("Low" = "low",
                                         "Moderate" = "mod",
                                         "Max" = "max"),
                             selected = "max"),
          hr(),
          actionButton("go", "Simulate Data")

      ),
      box(
        shiny::selectInput("guess.dist", label = "Guess for distribution",
                           choices = c( "exp","weibull", "weibullPH", "gompertz",
                                                  "gamma", "llogis", "lnorm"),
                           selected = "Mixed"),

        shiny::selectInput("guess.ph", label = "Guess for whether PH/AFT is satisfied",
                           choices = c("Yes", "No"),
                           selected = "Yes"),

        hr(),

        shiny::actionButton("guess", label = "Make Guess")
      ),
      tabBox(
       width = 12,
        tabPanel("KM",plotOutput("plot", height = "700px")),
        tabPanel("Cumulative Hazard Plots", plotOutput("cumhaz", height = "700px")),
       tabPanel("QQ Plots", plotOutput("qq", height = "700px")),
       tabPanel("Smoothed hazards", plotOutput("smooth", height = "700px")),
       tabPanel("Schoenfeld residuals test and plot", plotOutput("zph", height = "700px")),
       tabPanel("Parametric Fits", tableOutput("fit.tab")),
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

  lapply(list.files("./R"), function(x) source(paste0("./R/", x)))

  res <- eventReactive(input$go, {

    shinyalert::shinyalert(text = "Simulating Data", closeOnEsc = FALSE, type = "info",
                           showConfirmButton = FALSE)

    dat <- .simsurv(n = input$n.sim, maturity = input$mat)


    tests <- .test.ph.aft(dat$dat)
  plot <- .km.plots(fits = tests$all.mods,
                                                  dat = dat$dat,
                                                  include = names(tests$all.mods))


  sm.haz.p <- .smooth.haz(
    dat = dat$dat,
    nbreaks = 10,
    fits = tests$all.mods,
    include = names(tests$all.mods)
  )
  shinyalert::closeAlert()
    dplyr::lst(dat, tests, plot, sm.haz.p)
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

  output$fit.tab <- renderTable({
    res()$tests$fits
  })

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
