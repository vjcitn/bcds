#' provide an app that develops gene set representations from TF target selections
#' @import shiny
#' @importFrom GSEABase setName geneIds
#' @export
tftFinder = function() {
 msigTFs = sort(sapply(TFutils::tftColl, setName))
 ui = fluidPage(
  sidebarLayout(
   sidebarPanel(width=3,
    helpText(h3("TF target finder.  Based on MSigDb c3/TFT gene sets.")),
    selectInput("build", "build", choices=c( "EnsDb.Hsapiens.v75", "EnsDb.Hsapiens.v79"),
                 selected = "EnsDb.Hsapiens.v75" ),
    selectInput("pickedTF", "TF id", choices=msigTFs, selected=msigTFs[1]),
    uiOutput("biotypes"),
    uiOutput("nschrom"),
    actionButton("btnSend", "stop app")
    ),  
   mainPanel(
    dataTableOutput("targets")
    )   
   )
  )
 server = function(input, output) {
  getGenes = reactive({
    library(input$build, character.only=TRUE)
    edb = get(input$build)
    genes(edb, filter=EntrezFilter(geneIds(TFutils::tftColl[[input$pickedTF]])))
    })
  output$targets = renderDataTable({
    gn = getGenes()
    std = as.character(c(1:22, "X", "Y"))
    chr2keep = std
    if (!is.null(input$nschr)) chr2keep = c(chr2keep, input$nschr)
    kp = which(gn$gene_biotype %in% input$biotype & as.character(seqnames(gn)) %in% chr2keep)
    gn <<- gn[kp]
    as.data.frame(gn[kp])
    })  
  output$biotypes = renderUI( {
    gn = getGenes()
    biotypes = sort(unique(gn$gene_biotype))
    selectInput("biotype", "biotypes to keep", choices=biotypes,
        selected=biotypes, multiple=TRUE)
    } )
  output$nschrom = renderUI( {   # nonstandard chromosomes
    gn = getGenes()
    std = as.character(c(1:22, "X", "Y"))
    allseq = unique(as.character(seqnames(gn)))
    nschr = sort(setdiff(allseq, std))
    selectInput("nschr", "non-std chr to keep", choices=nschr,
        multiple=TRUE)
    } )

   observe({
            if(input$btnSend > 0)
               isolate({
                 stopApp(returnValue=0)
                        })  
           })  

 }
 runApp(list(ui=ui, server=server))
 gn
}

