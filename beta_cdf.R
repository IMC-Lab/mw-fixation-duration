library(tidyverse)
library(shiny)

ui <- fluidPage(titlePanel("Beta CDF"),
                fluidRow(column(2,
                                sliderInput(inputId="N_timer", label="N_timer", min=1, max=30, value=15),
                                sliderInput(inputId="N_labile", label="N_labile", min=1, max=30, value=15),
                                sliderInput(inputId="mod", label="mod", min=0, max=1, value=.5)),
                         column(5, plotOutput(outputId='cdfPlot')),
                         column(5, plotOutput(outputId='differencePlot'))))

server <- function(input, output) {
    d <- reactive({
        expand_grid(alpha=input$N_timer, beta=input$N_labile, mod=c(1, input$mod),
                    x=seq(0, 1, length.out=101)) %>%
            mutate(y=pbeta(x, mod*alpha, mod*beta)) %>%
            mutate(mod=ifelse(mod == 1, 'baseline', 'modulation'))
    })
    
    output$cdfPlot <- renderPlot({
        ggplot(d(), aes(x=x, y=y, color=mod, group=mod)) +
            geom_vline(xintercept=((input$N_timer) / (input$N_timer + input$N_labile))) +
            geom_line(size=1) +
            theme_bw()
    })

    output$differencePlot <- renderPlot({
        if (input$mod == 1.0) {
            return(tibble(x=seq(0, 1, length.out=101), difference=0) %>%
                   ggplot(aes(x=x, y=difference)) +
                   geom_vline(xintercept=((input$N_timer) / (input$N_timer + input$N_labile))) +
                   geom_line(size=1) +
                   geom_area(alpha=0.25) +
                   theme_bw())
        }
        
        d() %>%
            pivot_wider(names_from=mod, values_from=y) %>%
            mutate(difference=modulation-baseline) %>%
            ggplot(aes(x=x, y=difference)) +
            geom_vline(xintercept=((input$N_timer) / (input$N_timer + input$N_labile))) +
            geom_line(size=1) +
            geom_area(alpha=0.25) +
            theme_bw()
    })
}

shinyApp(ui, server)


expand_grid(alpha=15, beta=15, mod=c(1, 1/3), x=seq(0, 1, length.out=101)) %>%
    mutate(y=pbeta(x, mod*alpha, mod*beta)) %>%
    mutate(mod=factor(ifelse(mod == 1, 'beta[N] == 1', 'beta[N] == 1/3'), levels=c('beta[N] == 1', 'beta[N] == 1/3'))) %>%
    ggplot(aes(x=x, y=y, color=mod)) +
    geom_line(size=1) +
    coord_fixed() + xlab(expression(r[timer] / r[timer] + r[labile])) + ylab('Cancellation Probability') +
    scale_color_discrete(name='', labels=parse_format()) +
    theme_bw()
ggsave('beta_cdf.pdf', width=6, height=4)
