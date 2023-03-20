library(shiny)
library(shinydashboard)
library(ggplot2)
library(doremi)

# Environment for global variables for the app
dec_env <- new.env()
dec_env$reactives <- reactiveValues(
    data = NULL,
    deconvoluted = NULL,
    params = NULL
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

calculate_params <- function(time, reporter, points_initial_rate) {
    #' Calculates the requested parameters from the deconvoluted data
    #' Params:
    #' ------
    #' time -> time of measurement
    #' reporter -> deconvoluted reported levels
    #' points_initial_rate -> number of points to use for calculating the initial rate

    initial_rate <- NA
    initial_rate_model <- NA
    maximum <- NA
    maximum_time <- NA
    decay_rate <- NA
    decay_rate_model <- NA

    # Calculate initial rate    
    model <- lm(reporter ~ time,
        data = data.frame(
            time = time[1:points_initial_rate],
            reporter = reporter[1:points_initial_rate]
        )
    )
    initial_rate <- coef(model)[2]
    initial_rate_model <- model

    # Calculate maximum rate
    maximum <- max(reporter, na.rm = TRUE)
    maximum_time <- time[which.max(reporter)]

    # Calculate decay rate
    # Get the trace from the maximum onwards, then fit an exponential decay curve
    max_index <- which.max(reporter)

    data <- data.frame(
        time = time[max_index:length(time)],
        reporter = reporter[max_index:length(time)]
    )
    # Remove NAs
    data <- data[complete.cases(data), ]

    fit <- nls(reporter ~ SSasymp(time, yf, y0, log_alpha),
        data = data
    )

    decay_rate <- exp(coef(fit)["log_alpha"])
    decay_rate_model <- fit

    params <- list(
        initial_rate = initial_rate,
        maximum = maximum,
        maximum_time = maximum_time,
        decay_rate = decay_rate,
        initial_rate_model = initial_rate_model,
        decay_rate_model = decay_rate_model
    )

    return(params)
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
                p("Upload a CSV file with a Time (in hours) and a Luminescence column."),
                p("Select the reporter halflife (luciferase: ~3.7 hours)"),
                p("Click 'Deconvolute' to analyse the data, then 'Download' to save the results.")
            )
        ),
        fileInput("input_file", "Choose CSV File",
            accept = c(".csv")
        ),
        numericInput("halflife", "Reporter halflife (hours)", 3.7, step = 0.1),
        numericInput("points_initial_rate", "Number of points for initial rate", 10, step = 1),
        actionButton("deconvolute", "Deconvolute"),
        div(style = "margin-left: 15px;", downloadButton("download", "Download"))
    ),
    dashboardBody(
        plotOutput("plot"),
        tableOutput("params_table")
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

        # Calculate parameters
        dec_env$reactives$params <- calculate_params(
            time = dec_env$reactives$data$Time,
            reporter = dec_env$reactives$deconvoluted()$Synthesis,
            points_initial_rate = input$points_initial_rate
        )
    })

    # Plot data
    output$plot <- renderPlot({
        # If deconvoluted data is available, plot it
        # Otherwise, plot the raw data, if available
        if (is.null(dec_env$reactives$deconvoluted)) {
            if (is.null(dec_env$reactives$data)) {
                g <- ggplot() +
                    theme_bw()
            } else {
                g <- ggplot(dec_env$reactives$data, aes(x = Time, y = Luminescence)) +
                    geom_line(col = "lightgray") +
                    labs(x = "Time (hours)", y = "Luminescence (a.u.)") +
                    theme_bw() +
                    theme(
                        axis.text = element_text(size = 14),
                        axis.title = element_text(size = 16)
                    )
            }
        } else {
            g <- ggplot(dec_env$reactives$deconvoluted(), aes(x = Time, y = Luminescence)) +
                geom_line(col = "lightgray") +
                geom_line(aes(y = Synthesis), linewidth = 1.1) +
                labs(x = "Time (hours)", y = "Luminescence (a.u.)") +
                theme_bw() +
                theme(
                    axis.text = element_text(size = 14),
                    axis.title = element_text(size = 16)
                )

            # If we calculated parameters, add them to the plot
            if (length(dec_env$reactives$params)) {
                # Maximum
                maximum <- dec_env$reactives$params[["maximum"]]
                maximum_time <- dec_env$reactives$params[["maximum_time"]]
                g <- g +
                    geom_point(aes(x = maximum_time, y = maximum), size = 3, col = "red")

                # Initial rate
                initial_rate_model <- dec_env$reactives$params[["initial_rate_model"]]
                prediction <- predict(initial_rate_model, newdata = data.frame(time = dec_env$reactives$data$Time[1:input$points_initial_rate]))
                new_data <- data.frame(Time = dec_env$reactives$data$Time[1:input$points_initial_rate], Luminescence = prediction)
                g <- g +
                    geom_line(data = new_data, aes(x = Time, y = Luminescence), col = "orange", linewidth = 1.3)

                # Decay rate
                decay_rate_model <- dec_env$reactives$params[["decay_rate_model"]]

                maximum_index <- which(dec_env$reactives$data$Time == maximum_time)
                prediction <- predict(decay_rate_model,
                    newdata = data.frame(time = dec_env$reactives$data$Time[maximum_index:nrow(dec_env$reactives$data)])
                )
                new_data <- data.frame(
                    Time = dec_env$reactives$data$Time[maximum_index:nrow(dec_env$reactives$data)],
                    Luminescence = prediction
                )
                g <- g +
                    geom_line(data = new_data, aes(x = Time, y = Luminescence), col = "#81d100", lty = "dashed", linewidth = 1.1)
            }
        }

        print(g)
    })

    output$params_table <- renderTable({
        req(dec_env$reactives$params)

        params <- data.frame(
            Parameter = c("Initial rate (AU/hour)", "Maximum (AU)", "Time to maximum (hours)", "Decay rate (1/hour)"),
            Value = unlist(dec_env$reactives$params[1:4])
        )
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
