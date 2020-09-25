#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)


# Define UI for application that draws a histogram
ui <- fluidPage(

  titlePanel("Selection (Evolution class, Jodie Wiggins)"),

  #fitnessSurface=density(rnorm(100000,mean=20,sd=1)),startingPopulation=200,mutationMean=0,mutationSD=0.1,ngen=25

  sidebarPanel(
    radioButtons("dist", "Distribution of fitness",
                 c("One peak" = 1,
                   "Two equal peaks" = 2,
                   "Major and minor peak" = 3,
                   "Flat with bounds" = 4)),
    sliderInput('startingPopulation', 'Starting population size', min=10, max=1000, value=200, round=1),
    sliderInput('mutationMean', 'Mean bias of mutations', min=-1.1, max=1.1, value=0, round=FALSE),
    sliderInput('mutationSD', 'Standard deviation of mutations', min=0, max=1, value=.1),
    sliderInput('ngen', 'Number of generations to simulate', min=1, max=100, value=25, round=1)


  ),

  mainPanel(
    plotOutput('plot')
  )
)


# Define server logic required to draw a histogram
server <- function(input, output) {

  output$plot <- renderPlot({

    #fitnessSurface <- input$fitnessSurface
    fitnessSurface <- density(rnorm(100000,mean=20,sd=1))
    if(input$dist == 1) {
      fitnessSurface <- density(rnorm(100000,mean=20,sd=1))
    }
    else if (input$dist == 2) {
      fitnessSurface <- density(c(rnorm(100000,20,1),rnorm(100000,25,1)))
    }
    else if (input$dist == 3) {
      fitnessSurface <- density(c(rnorm(100000,20,1),rnorm(80000,25,1.5)))
    }
    else if (input$dist == 4) {
      fitnessSurface <- list(x=seq(from=10, to=20, length.out=10000), y=rep(0.5, 10000))
    }
    startingPopulation <- input$startingPopulation
    mutationMean <- input$mutationMean
    mutationSD <- input$mutationSD
    ngen <- input$ngen

    fitnessSurface$y<-fitnessSurface$y/max(fitnessSurface$y)
    popSize<-length(startingPopulation)
    if(length(startingPopulation)==1) { #just have a single value
      popSize<-round(startingPopulation[1])
      startingPopulation<-runif(n=popSize,min=min(fitnessSurface$x),max=max(fitnessSurface$x))
    }
    currentPopulation<-startingPopulation
    plot(x=0.9*c(min(fitnessSurface$x),1.1*max(fitnessSurface$x)),y=c(1.1*ngen,0),type="n",bty="n",yaxt="n",ylab="",xlab="trait")
    for (i in 1:(length(fitnessSurface$y)-1)) {
      meanY=mean(fitnessSurface$y[i:i+1])
      rect(xleft=fitnessSurface$x[i],xright=fitnessSurface$x[i+1],ybottom=0,ytop=ngen+0.1*ngen*meanY/max(fitnessSurface$y),border=NA,col=rgb(1-meanY/4,1-meanY/4,1-meanY/4))
    }
    lines(x=fitnessSurface$x,y=ngen+0.1*ngen*fitnessSurface$y/max(fitnessSurface$y))
    lines(x=fitnessSurface$x,y=ngen+0.1*0*fitnessSurface$y)

    for (generation in 1:ngen) {
      populationSurface<-density(currentPopulation)
      fitness<-(approx(fitnessSurface$x,fitnessSurface$y,xout=currentPopulation,yleft=0,yright=0))$y
      if(max(fitness)==0) {
        print("Your population has gone extinct, probably mutation moved it too far from the peak -- everything had zero fitness")
        break
      }
      selectedParents<-as.vector(rmultinom(1,popSize,fitness))
      selectedParentsToPrune<-selectedParents
      newPopulation<-c()
      while (max(selectedParentsToPrune)>0) {
        particularParent<-which.max(selectedParentsToPrune)
        childState<-rnorm(1,mean=currentPopulation[particularParent]+mutationMean,sd=mutationSD)
        newPopulation<-append(newPopulation,childState)
        lines(x=c(currentPopulation[particularParent],childState),y=c(generation-1, generation),col=rgb(1,0,0,0.3))
        selectedParentsToPrune[particularParent]= selectedParentsToPrune[particularParent]-1
      }
      lines(x=c(median(currentPopulation),median(newPopulation)),y=c(generation-1, generation),col=rgb(0,0,1,0.9),lwd=2)
      lines(x=c(quantile(currentPopulation,0.05),quantile(newPopulation,0.05)),y=c(generation-1, generation),col=rgb(0,0,1,0.9),lwd=1)
      lines(x=c(quantile(currentPopulation,0.95),quantile(newPopulation,0.95)),y=c(generation-1, generation),col=rgb(0,0,1,0.9),lwd=1)

      currentPopulation<-newPopulation
    }
  })

}

# Run the application
shinyApp(ui = ui, server = server)
