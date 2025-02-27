#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define server logic required to draw a histogram
function(input, output, session) {
  rv <- reactiveValues(coords.regions = "chr10:86652547-89647572", interactions = NULL, distance = 2000)
  observeEvent(input$upload.bedpe, {
    file.coords <- input$upload.bedpe
    if (is.null(file.coords)) {
      return(NULL)
    }
    exp_name <- input$sampleID
    coords.regions <- fread(file.coords$datapath)[,1:7,with=F]
    colnames(coords.regions) <- c("chr1", "start1", "end1", "chr2", "start2", "end2", "significance")
    coords.regions$sample <- exp_name
    if (is.null(rv$interactions)){
      rv$interactions <- coords.regions
      rv$interactions$sample <- factor(as.vector(rv$interactions$sample),levels=unique(rv$interactions$sample))
    } else {
      rv$interactions <- rbind(rv$interactions,coords.regions)
      rv$interactions$sample <- factor(as.vector(rv$interactions$sample),levels=unique(rv$interactions$sample))
    }

  })

  observeEvent(input$choose.region,{
    rv$coords.regions <- input$region.plot
  })
  observeEvent(input$reset.region,{
    rv$coords.regions <- "chr10:86652547-89647572"
  })

  observeEvent(input$reset.plot,{
    rv$interactions <- NULL
  })

  observeEvent(input$choose.dist,{
    rv$distance <- as.numeric(input$distance)*1000
  })

  genes <- reactive({
    coords <- fread("data/genes_hg38.bed")
  })

  enhancers <- reactive({
    coords <- fread("data/ATAC_selected.bed")
  })

  significance.plot <- reactive({
    interactions <- rv$interactions
    dist.min <- ifelse(is.na(rv$distance), 2000, as.numeric(rv$distance))
    coords.plot <- as.numeric(unlist(strsplit(rv$coords.regions,":|-"))[2:3])
    genes <- genes()
    coords <- c(87861000,87867000)
    if (is.null(interactions)){
      g <- ggplot()
    } else if (length(unique(interactions$sample))==1) {
      interactions$distance <- abs(interactions$start2 - interactions$start1)
      interactions.arc <- interactions[distance > dist.min]
      interaction.sig <- interactions[(start1>coords[1] & end1<coords[2]) | (start2>coords[1] & end2<coords[2])]
      g <- ggplot(interaction.sig,aes(x = start2, y = significance)) +
        geom_point() +
        theme_light() +
        labs(x="",y="") +
        coord_cartesian(xlim = coords.plot,ylim=c(0,100))
    } else {
      palette.samp <- colours()[c(24,124,144,33,613,367,57,511,227)]
      interactions$distance <- abs(interactions$start2 - interactions$start1)
      interactions.arc <- interactions[distance > dist.min]
      interaction.sig <- interactions[(start1>coords[1] & end1<coords[2]) | (start2>coords[1] & end2<coords[2])]
      g <- ggplot(interaction.sig,aes(x = start2, y=sample,color = sample, alpha=significance)) +
        geom_jitter() +
        theme_light() +
        #scale_colour_gradient2(low = "yellow",mid ="red",midpoint = 100 , high = "black") +
        scale_colour_manual(values=palette.samp) +
        scale_alpha_continuous(range=c(0,50)) +
        labs(x="",y="") + theme(legend.position = "top") +
        coord_cartesian(xlim = coords.plot)
    }
    return(g)
  })

  interaction.plot <- reactive({
    interactions <- rv$interactions
    dist.min <- ifelse(is.na(rv$distance), 2000, as.numeric(rv$distance))
    coords.plot <- as.numeric(unlist(strsplit(rv$coords.regions,":|-"))[2:3])
    genes <- genes()
    coords <- c(87861000,87867000)
    print(length(unique(interactions$sample)))
    if (is.null(interactions)){
      g <- ggplot()
    } else {
      palette.samp <- colours()[c(24,124,144,33,613,367,57,511,227)]
      interactions$distance <- abs(interactions$start2 - interactions$start1)
      interactions.arc <- interactions[distance > dist.min]
      interaction.sig <- interactions[(start1>coords[1] & end1<coords[2]) | (start2>coords[1] & end2<coords[2])]
      if (input$merge.ints){
        g <- ggplot(interactions.arc,aes(x = start1, y = 0)) +
          geom_curve(aes(xend = start2, yend = 0, alpha=significance, color=sample)) +
          scale_colour_manual(values=palette.samp) +
          theme_light() +
          labs(x="",y="") + theme(legend.position = "none") +
          coord_cartesian(xlim = coords.plot,ylim = c(-50,0))
      } else {
        g <- ggplot(interactions.arc,aes(x = start1, xend= start2, y = as.numeric(sample),yend  = as.numeric(sample))) +
          geom_curve(aes(alpha=significance, color=sample)) +
          scale_colour_manual(values=palette.samp) +
          scale_y_continuous(breaks = 1:length(levels(interactions.arc$sample)),labels = levels(interactions.arc$sample)) +
          theme_light() +
          expand_limits(y = c(0,length(levels(interactions.arc$sample)) + 0.5)) +
          labs(x="",y="") + theme(legend.position = "none") +
          coord_cartesian(xlim = coords.plot)
    }

    }
    return(g)
  })

  gene.plot <- reactive({
    genes <- genes()
    coords.plot <- as.numeric(unlist(strsplit(rv$coords.regions,":|-"))[2:3])
    g <- ggplot(genes,aes(xmin=start1,xmax=end1,ymax=-2,ymin=-6)) +
      geom_rect(fill=c(colours()[124],"black",colours()[144]),alpha=0.5) +
      theme_void() +
      labs(x="",y="") + theme(legend.position = "none") +
      coord_cartesian(xlim = coords.plot,ylim=c(-6.1,-1.9))
    return(g)

  })

  enhancer.plot <- reactive({
    enh <- enhancers()
    coords.plot <- as.numeric(unlist(strsplit(rv$coords.regions,":|-"))[2:3])
    g <- ggplot(enh,aes(xmin=start1,xmax=end1,ymax=-2,ymin=-6)) +
      geom_rect(fill=rep(c(colours()[124],colours()[613]),ceiling(nrow(enh)/2)),alpha=0.5) +
      theme_void() +
      labs(x="",y="") + theme(legend.position = "none") +
      coord_cartesian(xlim = coords.plot,ylim=c(-6.1,-1.9))
    return(g)
  })

  output$chromInt <- renderPlot({
    p1 <- interaction.plot()
    p2 <- significance.plot()
    p3 <- gene.plot()
    p4 <- enhancer.plot()
    if (input$merge.ints | (!is.null(rv$interactions) & length(unique(rv$interactions$sample))<=2)){
      p2 + p3 + p4 + p1 + plot_layout(heights = c(6, 0.5,0.5, 4), ncol=1)
    } else {
      if (!is.null(rv$interactions) & length(unique(rv$interactions$sample))>2){
        l <- length(unique(rv$interactions$sample)) - 2
        p2 + p3 + p4 + p1 + plot_layout(heights = c(6, 0.5,0.5, 4+(0.5*l)), ncol=1)
      }
    }
  },height = 1000)

  #output$significance.plot <- renderPlot({
  #  significance.plot()
  #})

}
