## Author: Nicholas Winn
## ntwinn23@gmail.com
## Visualization of Transcriptional Reversion of Cardiac Myocyte Fate During Mammalian Cardiac Regeneration

# Import libraries
library(shiny)
library(bslib)
library(ggplot2)
library(colourpicker)
library(DT)
library(bslib)
library(tidyr)
library(tidyverse)
library(stringr)
#install.packages("shinycssloaders")
library(shinycssloaders)
library(edgeR)
library(gplots)
#install.packages("ggbeeswarm")
library(ggbeeswarm)
library(patchwork)
#install.packages("DESeq2")
library(DESeq2)





                                ##### Define UI #####





ui <- fluidPage(
  theme = bs_theme(version = 5, bootswatch = "darkly"),
    titlePanel("BF591 RShiny Assignment"),
      sidebarLayout(
        sidebarPanel(
          fileInput(inputId = "upload", label = "Upload a CSV", accept = ".csv"),
          actionButton("submit_button", "Submit", icon = NULL),
          sliderInput("variance_slider", "Adjust the variance", 
                      min = 0, max = 100, step = 1, value = 50),
          sliderInput("non_zero_slider", "Genes that contain non-zero samples",
                      min = 0, max = 8, step = 1, value = 4)
        ),
        mainPanel(style='height: 90vh; overflow-y: auto;', 
          textOutput(outputId = "fileInfo"),
          actionButton("change_button", "Apply Changes", icon = NULL),
          tabsetPanel(
            tabPanel("Sample",
                      tableOutput("sample"),
              tabsetPanel(
                tabPanel("Metadata",
                      withSpinner(tableOutput("metadata"))
                ),
                tabPanel("Summary Table",
                      withSpinner(tableOutput("summary"))
                ),
                tabPanel("Histogram",
                         withSpinner(plotOutput("histogram"))
                )
              )
            ),
            tabPanel("Count",
              tabsetPanel(
                tabPanel("Statistics",
                  
                  textOutput("num_samples"),
                  textOutput("num_genes"),
                  textOutput("num_genes_passing_filter"),
                  textOutput("num_genes_not_passing_filter")
                ),
                 tabPanel("Scatter Plot",
                          withSpinner(plotOutput("scatter_variance")),
                          withSpinner(plotOutput("scatter_zeros"))
                 ),
                 tabPanel("Heatmap",
                          withSpinner(plotOutput("heatmap"))
                 ),
                 tabPanel("PCA",
                          withSpinner(plotOutput("pca")),
                          sliderInput("pca_slider", "Select desired principal components",
                                      min = 1, max = 6, step = 1, value = 1)
                 )
               )
            ), 
            tabPanel("DE",
              tabsetPanel(
                tabPanel("Table",
                     withSpinner(DT::dataTableOutput('de'))
                     ),
                tabPanel("Volcano Plot",
                     withSpinner(plotOutput("volcano")),
                     radioButtons("button_x", "Choose a column for the X-axis", 
                                  choices = c("baseMean", "log2FoldChange", "pvalue", "padj")),
                     radioButtons("button_y", "Choose a column for the Y-axis",
                                  choices = c("baseMean", "log2FoldChange", "pvalue", "padj")),
                     colourInput("base", "Base point color", value = "#D4A5D9"),
                     colourInput("highlight", "Highlight point color", value = "#50C7B5"),
                     sliderInput("scale_adjust", "Select the magnitude of adjustment", 
                                 min = -300, max = 0, step = 10, value = -100)
                )
              )
            ),
            tabPanel("Gene Counts",
                     withSpinner(plotOutput("gene_counts")),
                     selectInput("categorical_variable", "Choose Categorical Variable",
                                 choices = c("Timepoint", "Replicate")),
                     textInput("gene_search", "Search Gene"),
                     selectInput("plot_type", "Choose Plot Type",
                                 choices = c("Bar Plot", "Boxplot", "Violin Plot", "Beeswarm Plot")),
                     actionButton("visualize_button", "Visualize Gene Counts")
            )
          )
        )
      )
)





                            ##### Define server #####





server <- function(input, output, session) {

# Increase data size for bigger data sets
  options(shiny.maxRequestSize = 30*1024^2)

  ####
  thematic::thematic_shiny()
  ####
  
  
##### FUNCTIONS #####  
  
  
  
# Metadata table function
  meta_table <- function(data) {
    sample_names <- colnames(data)
    timepoints <- str_sub(sample_names, 2, 3)
    samples <- str_sub(sample_names, -1, -1)
  
    meta_tibble <- tibble(
      sample = sample_names,
      timepoint = timepoints,
      replicate = samples
    )
  
    return(meta_tibble)  
  } 
  
# Counts per million normalization function
  normalize_cpm <- function(datan) {
    count_data <- datan
    cpm_data <- count_data / colSums(count_data) * 1e6
    return(cpm_data)
  }
  
# DESEQ function
  deseq_function <- function(datad) {
    col_data <- meta_table(datad)
    
    dds_matrix <- DESeqDataSetFromMatrix(countData = datad,
                                         colData = col_data,
                                         design = ~timepoint)
    
    dds <- DESeq(dds_matrix)
    
    return(dds)
  }
  
# Volcano Plot function
  volcano_plot <-function(dataf, x_name, y_name, slider, color1, color2) {
    if (y_name=='padj'){
      dataf$padj<- log10(dataf$padj)
    }
    else {
      slider<-10^slider
    }
    p <- ggplot(dataf, aes(x = !!sym(x_name), y = !!sym(y_name), color = padj <= slider)) +
      geom_point() +
      scale_color_manual(values = c('FALSE'=color1, 'TRUE'=color2)) +
      labs(x = x_name, y = y_name, color=paste0('padj<10^',slider), title = "Volcano Plot")#+
    ifelse(y_name=='padj',return(p+scale_y_reverse()),return(p))
  }   
  
filter_data <- function(datal){
  data <-  datal %>% column_to_rownames(var='gene')
  
  row_variances <- apply(data, 1, var, na.rm = FALSE)
  non_zero_counts <- rowSums(data != 0, na.rm = FALSE)

  fd<-data %>% filter(row_variances>quantile(row_variances,na.rm=TRUE,probs=(input$variance_slider/100))
                      ,rowSums(data>0)>=input$non_zero_slider)

  return(fd)
}
 
  
  
### REACTIVES ###
  
  
  
# Sample data reactive
  load_data <- reactive({
    req(input$upload)
    data <- read_csv(input$upload$datapath)
    return(data)
  })
  
# Filtered data reactive
  filtered <- reactive({
    data_load <- read_csv(input$upload$datapath)
    f <- filter_data(data_load)
    return(f)
  })
  
# DESEQ results reactive  
  deseq_results <- reactive ({
    data_load <- read_csv(input$upload$datapath)

    data <- filter_data(data_load)

    deseq_dds <- deseq_function(data)
    
    d_results <- results(deseq_dds)
    
    return(as.data.frame((d_results)))
  })
  
  
  
##### OUTPUTS #####
  
  
  
## SUMMARY TAB OUTPUTS ## 
  
  
  # Metadata Table
  output$metadata <- renderTable({
    data <- load_data()

      sample_names <- colnames(data[-1])
      timepoints <- str_sub(sample_names, 2, 3)
      samples <- str_sub(sample_names, -1, -1)
      
      meta_tibble <- tibble(
        sample = sample_names,
        timepoint = timepoints,
        replicate = samples
      )
      
      return(meta_tibble)
  })
  
  # Summary Table
  output$summary <- renderTable({
    data <- load_data()
    select_data <- data %>%
      select(-1)
    
    col_names <- colnames(select_data)
    col_means <- colMeans(select_data)
    col_sd <- apply(select_data, 2, sd)
    col_max <- apply(select_data, 2, max)
    col_min <- apply(select_data, 2, min)
    
    
    summary_table <- tibble(
      Columns = col_names,
      Means = col_means,
      SD = col_sd,
      Max = col_max,
      Min = col_min
    )
    
    return(summary_table)
  })
  
  
  # Histogram Plot
  output$histogram <- renderPlot({
    data <- load_data()
      df_long <- data %>%
        select(-(1)) %>%
        pivot_longer(cols = everything(), names_to = "Timepoint", values_to = "Counts")
      
      p <- df_long %>%
        mutate(FirstLetter = substr(Timepoint, 1, 1)) %>%
        ggplot(aes(x = Timepoint, y = Counts, fill = Timepoint)) +
        geom_bar(stat = "identity", position = "stack", alpha = 0.7) +
        labs(title = "Bar Plot of Gene Expression Counts", x = "Condition", y = "Count") +
        #theme_minimal() +
        facet_wrap(~FirstLetter, scales = "free_x", ncol = 1)
      
      return(p)
  })   
  
  
## COUNTS TAB OUTPUTS ##
  
  
  # Number of Samples
  output$num_samples <- renderText({
    paste("Number of Samples:", ncol(filtered()))
  })
  
  # Number of Genes
  output$num_genes <- renderText({
    paste("Number of Genes:", nrow(load_data()))
  })
  
  # Number and % of Genes Passing Current Filter
  output$num_genes_passing_filter <- renderText({
    data <- normalize_cpm(filtered())
    total_genes <- nrow(load_data())
    genes_passing_filter <- nrow(data)
    
    percentage_passing_filter <- (genes_passing_filter / total_genes) * 100
    
    paste("Number of Genes Passing Current Filter:", genes_passing_filter,
          "(", round(percentage_passing_filter, 2), "%)")
  })  
  
  # Number and % of Genes Not Passing Current Filter
  output$num_genes_not_passing_filter <- renderText({
    data <- normalize_cpm(filtered())
    total_genes <- nrow(load_data())
    genes_passing_filter <- nrow(data)
    genes_not_passing_filter <- total_genes - genes_passing_filter
    
    percentage_not_passing_filter <- (genes_not_passing_filter / total_genes) * 100
    
    paste("Number of Genes Not Passing Current Filter:", genes_not_passing_filter,
          "(", round(percentage_not_passing_filter, 2), "%)")
  })
  
  # Variance slider text
  output$variance_slider_text <- renderText({
    paste("Variance Slider Value:", input$variance_slider)
  })
  
  # Non-zero slider text
  output$non_zero_slider_text <- renderText({
    paste("Non-zero Slider Value:", input$non_zero_slider)
  })
  
  # Scatter Plot for Variance
  output$scatter_variance <- renderPlot({
    data <- normalize_cpm(filtered())
    row_variances <- apply(data, 1, var, na.rm = FALSE)
    non_zero_counts <- rowSums(data != 0, na.rm = FALSE)
    
    plot_data <- data.frame(
      MedianCount = apply(data, 1, median, na.rm = TRUE),
      Variance = row_variances
    )
    
    ggplot(plot_data, aes(x = log10(Variance), y = log10(MedianCount))) +
      geom_point(color = 'lightblue', alpha = 0.7) +
      labs(
        title = "Median Count vs Variance",
        x = "Variance (log scale)",
        y = "Median Count (log scale)"
      ) 
  })
  
  # Scatter plot for number of zeros
  output$scatter_zeros <- renderPlot({
    data <- normalize_cpm(filtered())
    row_variances <- apply(data, 1, var, na.rm = FALSE)
    non_zero_counts <- rowSums(data != 0, na.rm = FALSE)
    
    plot_data <- data.frame(
      MedianCount = apply(data, 1, median, na.rm = TRUE),
      NumZeros = non_zero_counts
    )
    
    ggplot(plot_data, aes(x = NumZeros, y = log10(MedianCount))) +
      geom_point(color = 'lightblue', alpha = 0.7) +
      labs(
        title = "Median Count vs Number of Zeros",
        x = "Number of Zeros (log scale)",
        y = "Median Count (log scale)"
      ) 
  })
  
  # Heatmap
  output$heatmap <- renderPlot({
    data <- normalize_cpm(filtered())
    
    heatmap_data <- t(data)  # Transpose the data for proper heatmap orientation
    heatmap_data <- log2(heatmap_data + 1)  # Log-transform the counts
    
    heatmap(heatmap_data,
            col = rev(colorRampPalette(c("blue", "orange"))(512)),
            scale = "row",  # Scale rows (genes)
            key = TRUE,  # Include a color key
            keysize = 1.0,  # Set size of color key
            trace = "none",  # Remove trace lines
            density.info = "none",  # Remove density plot
            main = "Clustered Heatmap of Counts"
    )
  })
  
  # PCA
  output$pca <- renderPlot({
    filtered_data <- normalize_cpm(filtered())
    
    non_zero_genes <- filtered_data[rowSums(filtered_data != 0) > 0, ]
    
    transposed_data <- t(non_zero_genes)
    
    pca <- prcomp(transposed_data, scale. = TRUE)
    
    pca_df <- as_tibble(pca$x) %>%
      select(1:input$pca_slider) %>%
      pivot_longer(everything(), names_to = "PC", values_to = "projection") %>%
      mutate(PC=fct_relevel(PC,str_c("PC",1:input$pca_slider))) %>%
      ggplot(aes(x = PC, y = projection, color = PC)) +
      geom_beeswarm() + labs(title="PCA Projection Plot")
      labs(title = "PCA Projection Plot") +
      theme(axis.text.x=element_text(angle=90,hjust=0,vjust=0.5))
    
    print(pca_df)
  })  
  

## DESEQ TAB OUTPUTS ##
  
  
  # DESEQ Table
  output$de <- DT::renderDataTable({
    des<-data.frame(deseq_results())%>%rownames_to_column(var='gene')
    DT::datatable(des, class='table-light',style='bootstrap5', width='auto', rownames=FALSE,
                  options = list(paging = TRUE, searching = TRUE, scrollX=TRUE,
                                 fixedColumns = TRUE, autoWidth = TRUE,
                                 ordering = TRUE, dom = 'Bfrtip'))
  })

  # Volcano Plot
  output$volcano <- renderPlot({
    volcano_plot(deseq_results(), input$button_x, input$button_y, input$scale_adjust, 
                 input$base, input$highlight)
  })
  

## ADVENTURE TAB OUTPUTS ##


  # Gene count visualizer
  output$gene_counts <- renderPlot({
    data <- filtered()
    
    gene_name <- input$gene_search
    
    if (gene_name %in% rownames(data)) {
      data <- data[gene_name, , drop = FALSE]
    } else {
      return(paste("Gene not found in the row names"))
    }
    
    categorical_variable <- input$categorical_variable
    plot_type <- input$plot_type
    
    
    # Extract the timepoints or replicates based on the chosen categorical variable
    if (categorical_variable == "Timepoint") {
      variable_values <- unique(str_extract(colnames(data)[-1], "(?<=v)\\w+"))
      grouping_variable <- str_sub(colnames(data)[-1], 2, 3)
    } else if (categorical_variable == "Replicate") {
      variable_values <- unique(str_extract(colnames(data)[-1], "\\d+$"))
      grouping_variable <- str_sub(colnames(data)[-1], -1, -1)
    } else {
      variable_values <- NULL
      grouping_variable <- NULL
    }
    
    # Prepare data for plotting
    plot_data <- gather(data, key = "Sample", value = "Count", -1) %>%
      mutate(!!categorical_variable := ifelse(categorical_variable == "Timepoint",
                                              str_extract(Sample, "(?<=v)\\w+"),
                                              str_extract(Sample, "\\d+$")),
             Group = interaction(variable_values, grouping_variable, sep = "_")) %>%
      filter(!is.na(!!as.name(categorical_variable)))
    
    # Generate the plot based on the chosen plot type
    if (plot_type == "Bar Plot") {
      ggplot(plot_data, aes(x = Group, y = Count, fill = Sample, group = Sample)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(title = gene_name, x = categorical_variable, y = "Count") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    } else if (plot_type == "Boxplot") {
      ggplot(plot_data, aes(x = Group, y = Count, group = Sample)) +
        geom_boxplot() +
        labs(title = gene_name, x = categorical_variable, y = "Count") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    } else if (plot_type == "Violin Plot") {
      ggplot(plot_data, aes(x = Group, y = Count, fill = Group)) +
        geom_violin() +
        labs(title = gene_name, x = categorical_variable, y = "Count") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    } else if (plot_type == "Beeswarm Plot") {
      ggplot(plot_data, aes(x = Group, y = Count, group = Sample)) +
        geom_beeswarm() +
        labs(title = gene_name, x = categorical_variable, y = "Count") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }
  })
}




##### Run the application #####
shinyApp(ui = ui, server = server)