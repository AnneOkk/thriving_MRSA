library(shiny)
library(shinyWidgets)

ui <- fluidPage(
        title = "Box Title",
        box(
        dropdownMenu = dropdown(
            width = "200px",
            icon = icon("gear")),

            materialSwitch(inputId = "Id079", label = "Color:"),
            materialSwitch(inputId = "666", label = "Display Goal:")
        ),

        textOutput("text")
    )

server <- function(input, output, session) {
    output$text <- renderText("Hello World!")
}

shinyApp(ui, server)
