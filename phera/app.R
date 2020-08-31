library(shiny)
library(tidyverse)
library(viridis)
library(viridisLite)

ui <- fluidPage(
  h1("roGFP2 responds curves from phera data"),
  hr(),
  p(strong("File formats:")),
  p(em("Phera data:"), "Export .csv file from MARS (ASCII). MARS settings for file export: Tick \"add header\" and \"transpose table\". "),
  p(em("Layout:"), "Save layout (table format) as a .csv file. File needs to include the following columns:"), 
  p(tags$div
    (tags$ul
      (tags$li("well (A1 to H12"))),
    (tags$ul
     (tags$li("strain (what protein is tagged with roGFP2?)"))),
    (tags$ul
     (tags$li("background (wt/deletion/BY/W ?etc.)"))),
    (tags$ul
     (tags$li("treatment"))),
    (tags$ul
     (tags$li("max (indicate max_ox eg. diamide / high H2O2, and max_red eg. DTT treated samples")))
    ),
  hr(),
  sidebarLayout(
    sidebarPanel(fileInput(inputId = "phera_data", 
                            label = "Phera data: choose CSV File", 
                            accept = ".csv"),
                  fileInput(inputId = "layout", 
                            label = "Plate layout: choose CSV File", 
                            accept = ".csv"),
                  
                  # #ratio or OxD
                  # selectInput(inputId = "y_axis",
                  #             label = "Plot roGFP2 ratio or OxD?", 
                  #             choices = c("OxD" = "OxD", "ratio" = "ratio"),
                  #             selected = NULL,
                  #             multiple = FALSE,
                  #             selectize = TRUE),
                  
                  #plot_title
                  textInput(inputId = "plot_title_ratio", 
                            label = "Plot title ratio", 
                            placeholder = "My wonderful phera data", 
                            value = "My wonderful phera data"),
                  
                  textInput(inputId = "plot_title_oxd", 
                            label = "Plot title OxD", 
                            placeholder = "My wonderful phera data", 
                            value = "My wonderful phera data"),
                  
                  #change y-axis values: 
                  numericInput(inputId = "y_min_ratio",
                              label = "y-axis min (ratio)",
                              value = 0, min = 0, max = 50),
                  numericInput(inputId = "y_max_ratio",
                              label = "y-axis max (ratio)",
                              value = 20, min = 0, max = 100),
                  
                  numericInput(inputId = "y_min_oxd",
                              label = "y-axis min (OxD)",
                              value = 0, min = -3, max = 3),
                  numericInput(inputId = "y_max_oxd",
                              label = "y-axis max (OxD)",
                              value = 1, min = -3, max = 3),
                  
                  #change time to be displayed
                  sliderInput(inputId = "time_max",
                              label = "time max [min]",
                              value = 140, min = 1, max = 500),
      ),
  mainPanel( plotOutput("plot_ratio"),
             plotOutput("plot_oxd"))
    )
)


server <- function(input, output) {
  output$plot_ratio <- renderPlot({
  #phera_data
    phera_data <- input$phera_data
    ext <- tools::file_ext(phera_data$datapath)
    req(phera_data)
    validate(need(ext == "csv", "Please upload a csv file"))
    
   #layout
    layout <- input$layout
    ext <- tools::file_ext(layout$datapath)
    req(layout)
    validate(need(ext == "csv", "Please upload a csv file"))
    
  #tidy data
    df <- read_csv(phera_data$datapath, skip = 6, col_names = FALSE)
    col_names <- t(df[1:4, 3:length(df)]) %>%  # transpose 
      as_tibble() %>%
      unite("well_group", V1, V2, V4) %>%
      pull(well_group)
    
    col_names <- c("data", "time", col_names)

    colnames(df) <- col_names
    df <- df[5:length(df$data), ]
    df_bc <- df %>%
      dplyr::filter(str_detect(data, "Blank corrected based on Raw Data"))%>%
      mutate(data = if_else(str_detect(data, "400"), "bc400", "bc485")) %>%
      pivot_longer(cols = 3:length(df), names_to = "sample")%>%
      pivot_wider(names_from = data, values_from = value) %>%
      separate(sample, into = c("well_row", "well_col", "group")) %>%
      unite("well", well_row, well_col, sep = "")

    #annotate
    layout <- read_csv(layout$datapath, col_names = TRUE)

    df_bc <- left_join(df_bc, layout, by = "well") %>%
      mutate(time = round(as.numeric(time), digits = 0))

    #remove below blank
    df_bc <- df_bc %>% filter(bc400 > 0 & bc485 > 0)%>%
      mutate(bc400 = as.numeric(bc400),
             bc485 = as.numeric(bc485),
             ratio = as.numeric(bc400)/as.numeric(bc485)) %>%
      group_by(well) %>%
      mutate(timepoints_n = n_distinct(time)) %>%
      ungroup()
    
    max_redox <- df_bc %>%
      dplyr::filter(treatment %in% c("DTT", "diamide")) %>%
      group_by(strain, background, treatment) %>%
      summarise(max(bc400), max(bc485))%>%
      pivot_wider(names_from = treatment, values_from = c("max(bc400)", "max(bc485)"))
    
    colnames(max_redox) <- c("strain", "background", "max_red_bc400", "max_ox_bc400", "max_red_bc485", "max_ox_bc485")
    
    df_bc <- df_bc %>%
      left_join(max_redox, by = c("strain", "background")) %>%
      mutate(OxD = (bc400 * max_red_bc485 - max_red_bc400 * bc485) / (bc400 * max_red_bc485 - bc400 * max_ox_bc485 + max_ox_bc400 * bc485 - max_red_bc400 * bc485))
    

    #plot
    df_bc %>%
      ggplot(aes(x = time, y = ratio))+
      geom_point(aes(colour = treatment), size = 0.4)+
      facet_wrap(vars(strain, background))+
      coord_cartesian(ylim = c(input$y_min_ratio, input$y_max_ratio), xlim = c(0, input$time_max))+
      scale_color_viridis_d()+
      labs(title = input$plot_title_ratio)+
      xlab("time [min]")+
      ylab("roGFP2 ratio")+
      theme_bw()
  })
  
  
  output$plot_oxd <- renderPlot({
    #phera_data
    phera_data <- input$phera_data
    ext <- tools::file_ext(phera_data$datapath)
    req(phera_data)
    validate(need(ext == "csv", "Please upload a csv file"))
    
    #layout
    layout <- input$layout
    ext <- tools::file_ext(layout$datapath)
    req(layout)
    validate(need(ext == "csv", "Please upload a csv file"))
    
    #tidy data
    df <- read_csv(phera_data$datapath, skip = 6, col_names = FALSE)
    col_names <- t(df[1:4, 3:length(df)]) %>%  # transpose 
      as_tibble() %>%
      unite("well_group", V1, V2, V4) %>%
      pull(well_group)
    
    col_names <- c("data", "time", col_names)
    
    colnames(df) <- col_names
    df <- df[5:length(df$data), ]
    df_bc <- df %>%
      dplyr::filter(str_detect(data, "Blank corrected based on Raw Data"))%>%
      mutate(data = if_else(str_detect(data, "400"), "bc400", "bc485")) %>%
      pivot_longer(cols = 3:length(df), names_to = "sample")%>%
      pivot_wider(names_from = data, values_from = value) %>%
      separate(sample, into = c("well_row", "well_col", "group")) %>%
      unite("well", well_row, well_col, sep = "")
    
    #annotate
    layout <- read_csv(layout$datapath, col_names = TRUE)
    
    df_bc <- left_join(df_bc, layout, by = "well") %>%
      mutate(time = round(as.numeric(time), digits = 0))
    
    #remove below blank
    df_bc <- df_bc %>% filter(bc400 > 0 & bc485 > 0)%>%
      mutate(bc400 = as.numeric(bc400),
             bc485 = as.numeric(bc485),
             ratio = as.numeric(bc400)/as.numeric(bc485)) %>%
      group_by(well) %>%
      mutate(timepoints_n = n_distinct(time)) %>%
      ungroup()
    
    max_redox <- df_bc %>%
      dplyr::filter(max %in% c("max_ox", "max_red")) %>%
      group_by(strain, background, treatment) %>%
      summarise(max(bc400), max(bc485))%>%
      pivot_wider(names_from = treatment, values_from = c("max(bc400)", "max(bc485)"))
    
    colnames(max_redox) <- c("strain", "background", "max_red_bc400", "max_ox_bc400", "max_red_bc485", "max_ox_bc485")
    
    df_bc <- df_bc %>%
      left_join(max_redox, by = c("strain", "background")) %>%
      mutate(OxD = (bc400 * max_red_bc485 - max_red_bc400 * bc485) / (bc400 * max_red_bc485 - bc400 * max_ox_bc485 + max_ox_bc400 * bc485 - max_red_bc400 * bc485))
    
    
    #plot
    df_bc %>%
      ggplot(aes(x = time, y = OxD))+
      geom_point(aes(colour = treatment), size = 0.4)+
      facet_wrap(vars(strain, background))+
      coord_cartesian(ylim = c(input$y_min_oxd, input$y_max_oxd), xlim = c(0, input$time_max))+
      scale_color_viridis_d()+
      labs(title = input$plot_title_oxd)+
      xlab("time [min]")+
      ylab("roGFP2 OxD")+
      theme_bw()
    
  })
}

shinyApp(ui = ui, server = server)