#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(GenomicRanges)
library(ggplot2)
library(data.table)
library(patchwork)

# Define UI for application that draws a histogram
fluidPage(

    # Application title
    titlePanel("Generate arc plots from Tiled-C results"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
          textInput("region.plot", "Enter the coordinates of the region to plot", value = "chr10:86652547-89647572"),
          p("The format of the region coordinates: chr:start-end"),
          p("By default, the minimum and maximum coordinates present in the file will be used."),
          splitLayout(actionButton("choose.region","Apply"),actionButton("reset.region","Reset")),
          textInput("distance", "Enter the minimal distance between interacting regions (Kb):", value = ""),
          actionButton("choose.dist","Apply")
        ),

        # Show a plot of the generated distribution
        mainPanel(
            h6("Upload the interaction data"),
            p("The file MUST contain at least seven columns in this order: "),
            p("- chr1, start1, end1 specifying the chromosome start and end position of the first region"),
            p("- chr2, start2, end2 specifying the chromosome start and end position of the second region"),
            p("- significance indicating the weight of the interactions"),
            splitLayout(textInput("sampleID", "Give a name to your sample before uploading:", value = "samp1"),
                        fileInput("upload.bedpe", "Upload interation coordinates in bedpe format", buttonLabel = "Upload")
                        ),
            checkboxInput("merge.ints","Merge interactions from several samples",value = TRUE),
            actionButton("reset.plot","Reset"),
            #plotOutput("significance.plot"),
            plotOutput("chromInt")
        )
    )
)
