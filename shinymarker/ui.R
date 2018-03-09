source("setup.R")

shinyUI(
    fluidPage(
        titlePanel("FSA Visualizer"),
        sidebarLayout(
            sidebarPanel(
                checkboxInput("projManage","Project management",FALSE),

                conditionalPanel(
                    condition = "input.projManage == true",
                    textInput('project',"Project (no spaces)","default")
                    ),
                conditionalPanel(
                    condition = "input.projManage == true",
                    fileInput('filename', 'Choose zipped archive of fsa files to upload',
                              accept = c(
                                  '.zip',
                                  '.ZIP',
                                  '.Zip'
                                  )
                              )
                    ),
                conditionalPanel(
                    condition = "input.projManage == true",
                    actionButton('hoseProject',"Trash project (will have to upload files again)")
                    ),
                
                uiOutput("chooseind"),
                actionButton("updateIndividuals","Press to update individuals list"),
                checkboxGroupInput("dye", 
                                   label = ("Which dye (loci in paren)"), 
                                   choices = c("6-FAM "=1,
                                       "VIC"=2,
                                       "NED"=3,"PET"=4),
                                   selected = 1),
#                radioButtons('dye',"Which color?",c("6-FAM"=1,"VIC"=2,"NED"=3,"PET"=4)),
                checkboxInput("standard","Display ladder",FALSE),
                textOutput("size"),
#                checkboxInput("show.alleles","Show previous allele calls",FALSE),
                sliderInput("xzoom","X-window range (in base-pairs)",min=10,max=500,value=c(100)),
                sliderInput("pointer","Center of window",min=75,max=250,value=c(140)),
                radioButtons('grid',"Draw lines at every X base",c(100,50,10,5,1)),
                sliderInput("yrange","Y-range",min=0,max=1,value=c(0,1)),
                sliderInput("order","Order of polynomial for size fitting",min=1,max=5,value=c(3)),
                sliderInput("min.signal","Min signal for ladder peaks",min=100,max=10000,value=c(500)),
                sliderInput("min.time","Clip times less than this for peak finding",min=1,max=3000,value=c(1200))

                ),
            mainPanel(
                tabsetPanel(
                    tabPanel("1 individual per chromatogram",
                             plotOutput('aligned.dye',width="100%",
                                        hover = hoverOpts(id="plot_hover"))
                             ),
                    tabPanel("Overlapping chromatograms",
                             plotOutput('sameplot',width="100%",
                                        hover = hoverOpts(id="plot_hover"))
                             )
                    )
                )
            )
        ) 
    )
