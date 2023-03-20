library(shiny)
library(shinydashboard)
library(ggplot2)
library(doremi)

# Environment for global variables for the app
dec_env <- new.env()
dec_env$reactives <- reactiveValues(
    data = NULL,
    deconvoluted = NULL
)

deconvol <- function(time, reporter, halflife) {
    #' Deconvolutes reporter luminescence data to protein synthesis
    #' Adapted from Brown et al. Biotechnol. J. 2008
    #' Does NOT take into account volume dilution due to cell duplication
    #' Params:
    #' ------
    #' time -> time of measurement
    #' reporter -> reporter measurement (e.g. luminescence)
    #' halflife -> reporter half-life (in the same units as time)

    # Convert to rate
    halflife <- log(2) / halflife

    # Calculate derivative using the Generalized Orthogonal Local Derivative method
    # (GOLD) adapted from Deboeck 2010
    dBdT <- calculate.gold(time = time, signal = reporter, embedding = 9, n = 2)
    res <- dBdT$dsignal[, 2] + halflife * reporter

    res <- data.frame(
        Time = time,
        Luminescence = reporter,
        Synthesis = res
    )

    res
}

# Define UI
ui <- dashboardPage(
    dashboardHeader(),
    dashboardSidebar(
        # This avoids the download button to appear grayed out.
        # See https://stackoverflow.com/questions/36314780/
        tags$style(".skin-blue .sidebar a { color: #444; }"),
        # Instructions
        div(
            div(
                style = "margin-left: 15px;", h3("Instructions"),
                p("Upload a CSV file with a Time and a Luminescence column."),
                p("Select the reporter halflife (luciferase: ~3.7 hours)"),
                p("Click 'Deconvolute' to analyse the data, then 'Download' to save the results.")
            )
        ),
        fileInput("input_file", "Choose CSV File",
            accept = c(".csv")
        ),
        numericInput("halflife", "Reporter halflife (same unit as time)", 3.7, step = 0.1),
        actionButton("deconvolute", "Deconvolute"),
        div(style = "margin-left: 15px;", downloadButton("download", "Download"))
    ),
    dashboardBody(
        plotOutput("plot")
    )
)

server <- function(input, output) {
    # Read CSV file when "Read CSV" button is clicked
    observeEvent(input$input_file, {
        contents <- read.csv(input$input_file$datapath)
        dec_env$reactives$deconvoluted <- NULL # Reset deconvoluted data
        dec_env$reactives$data <- contents
    })

    # Deconvolute data when "Deconvolute" button is clicked
    observeEvent(input$deconvolute, {        
        req(dec_env$reactives$data)

        dec_env$reactives$deconvoluted <- reactive({
            deconvol(
                time = dec_env$reactives$data$Time,
                reporter = dec_env$reactives$data$Luminescence,
                halflife = input$halflife
            )
        })
    })

    # Plot data
    output$plot <- renderPlot({
        # If deconvoluted data is available, plot it
        # Otherwise, plot the raw data, if available
        print(dec_env$reactives$data)
        if (is.null(dec_env$reactives$deconvoluted)) {
            if (is.null(dec_env$reactives$data)) {
                g <- ggplot() +
                    theme_bw()
            } else {
                g <- ggplot(dec_env$reactives$data, aes(x = Time, y = Luminescence)) +
                    geom_line(col = "lightgray") +
                    labs(x = "Time (hours)", y = "Luminescence (a.u.)") +
                    theme_bw()
            }
        } else {
            g <- ggplot(dec_env$reactives$deconvoluted(), aes(x = Time, y = Luminescence)) +
                geom_line(col = "lightgray") +
                geom_line(aes(y = Synthesis), linewidth = 1.1) +
                labs(x = "Time (hours)", y = "Luminescence (a.u.)") +
                theme_bw()
        }

        print(g)
    })

    # Download deconvoluted data
    output$download <- downloadHandler(
        filename = function() {
            "deconvoluted.csv"
        },
        content = function(file) {
            write.csv(dec_env$reactives$deconvoluted(), file)
        }
    )
}

# Run the app
shinyApp(ui = ui, server = server)
