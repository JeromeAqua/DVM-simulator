#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyWidgets)

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  h1('DVM simulator'),
  p('Uses game theory to derive the optimal vertical migration patterns for a population of zooplankton prey and a population of forage fish predators.')
  ,
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
   #   sliderInput("bmax",
    #              "Maximum clearance rate of fish (m^3 / day)",
     #             min = 10^-2,
      #            max = 100,
       #           step=0.01,
        #          value = 1)
      #,
      
      shinyWidgets::sliderTextInput("bmax","Maximum clearance rate of fish (m^3 / day)",
                                    choices=c(10^-2, 0.05, 0.1, 0.5, 1, 5, 10, 50, 100),
                                    selected=1, grid = T),
      shinyWidgets::sliderTextInput("gmax","Maximum growth rate of prey (m^3 / day)",
                                    choices=10^-2*(1:100),
                                    selected=0.1, grid = T)
      ,
      sliderInput("N0",
                  "Average concentration of zooplankton in the water column (m^-3)",
                  min = 100,
                  max = 10000,
                  step = 10,
                  value = 5000)
      ,
      sliderInput("P0",
                  "Average concentration of fish in the water column (m^-3)",
                  min = 0.1,
                  max = 10,
                  step = 0.1,
                  value = 1)
      ,
      sliderInput("M",
                   "Number of depth bins (-)",
                   min = 2,
                   max = 100,
                   value = 30,
                   step = 1),
      sliderInput("sigma",
                  "Fraction of daytime during a day (-)",
                  min = 0,
                  max = 1,
                  value = 0.6,
                  step = 0.1),
      sliderInput("H",
                  "Depth of the water column (m)",
                  min = 20,
                  max = 1000,
                  value = 300,
                  step = 10),
      sliderInput("Lmax",
                  "Light irradiance at the surface during daytime (W/m^2)",
                  min = 20,
                  max = 800,
                  value = 100,
                  step = 10),

      shinyWidgets::sliderTextInput("L0","Half saturation parameter for the fish visual range (W/m^2)",
                                    choices=10^-2*(1:1000),
                                    selected=1, grid = T),
      shinyWidgets::sliderTextInput("klight","Light attenuation coefficient (m^-1)",
                                    choices=10^-2*(1:100),
                                    selected=0.07, grid = T),
      shinyWidgets::sliderTextInput("rho","Fractional difference between night and day light levels",
                                    choices=c(10^-6,10^-5,10^-4,10^-3,10^-2),
                                    selected=10^-3, grid = T),
      
      shinyWidgets::sliderTextInput("mu","Density-dependent mortality rate for predators (m^3/day)",
                                    choices=c(10^-6,10^-5,10^-4,10^-3,10^-2,10^-1,1,10),
                                    selected=10^-3, grid = T),
      shinyWidgets::sliderTextInput("C","Unitary migration cost (m^-1 day^-1)",
                                    choices=c(10^-8,10^-7,10^-6,10^-5,10^-4,10^-3),
                                    selected=10^-5, grid = T),
    sliderInput("eta",
                "Predator growth efficiency",
                min = 10^-4,
                max = 10^-1,
                value = 10^-2,
                step = 10^-2),
    sliderInput("zo",
                "Depth of mixed layer (m)",
                min = 0,
                max = 200,
                value = 50,
                step = 10),
    sliderInput("zm",
                "Thickness of the mixed layer transition zone",
                min = 1,
                max = 100,
                value = 10,
                step = 1)
    
  ),
    
    # Show a plot of the generated positions
    mainPanel(
      plotOutput(outputId = "Positions", height = "700px")
    )))

# Define server logic required to compute the positions
server <- function(input, output) {
  
  #Create the output function to draw the plots
  output$Positions <- renderPlot({
 
  dH <- input$H / input$M;
  
  #Create necessary vectors for the calculations
  Vi <- 1:input$M; #water layers
  Ai <-  do.call("cbind", replicate(input$M, matrix(Vi, nrow = input$M, ncol = 1), simplify = FALSE));
  zi <- -Vi*dH + dH/2;
  
  #Create growth ,light and so on functions
  #l <- function(z) {    output <- exp(input$klight * z)  };
  #g <- function(z) {    output <- (1- tanh((-input$zo-z)/input$zm))/2  };
  #bday <- function(z) {    output <- input$Lmax * exp(input$klight * z) / (input$L0 + input$Lmax*exp(input$klight * z)) };
  #bnight <- function(z) {    output <- input$rho * input$Lmax * exp(input$klight * z) / (input$L0 + input$rho*input$Lmax*exp(input$klight * z)) };
  
  v <- input$bmax * input$Lmax * exp(input$klight * zi) / (input$L0 + input$Lmax*exp(input$klight * zi));
  o <- input$bmax * input$rho * input$Lmax * exp(input$klight * zi) / (input$L0 + input$rho*input$Lmax*exp(input$klight * zi));
  g <- input$gmax * (1- tanh((-input$zo-zi)/input$zm))/2;
  
  
  #initialization
  N <- matrix(1, nrow = input$M, ncol = input$M);
  Nsum <- sum(sum(N));
  N <- N / Nsum;
  P <- N;
  
  Nday <- apply(N,2,sum);
  Nnight <- matrix(apply(N,1,sum), nrow = input$M, ncol = 1);
  Pday <- apply(P,2,sum);
  Pnight <- matrix(apply(P,1,sum), nrow = input$M, ncol = 1); 
  
  CNP <- input$C * dH * abs(Ai-t(Ai)); #migration cost, assumed the same for prey and predator for now
  GN <- (1-input$sigma) * do.call("cbind", replicate(input$M, matrix(g, nrow = input$M, ncol = 1), simplify = FALSE)) + 
           input$sigma* t(do.call("cbind", replicate(input$M, matrix(g, nrow = input$M, ncol = 1), simplify = FALSE))) - CNP ;
  MN <- input$M * input$P0* ((1-input$sigma)* do.call("cbind", replicate(input$M, matrix(o, nrow = input$M, ncol = 1)*Pnight, simplify = FALSE)) +
                               input$sigma*do.call("rbind", replicate(input$M, Pday*v, simplify = FALSE)));
  GP <- input$M * input$eta * input$N0 * ((1-input$sigma)*do.call("cbind", replicate(input$M, matrix(o, nrow = input$M, ncol = 1)*Nnight, simplify = FALSE)) +
                                             input$sigma * do.call("rbind", replicate(input$M, Nday*v, simplify = FALSE))) - CNP;
  MP <- input$mu * input$M * input$P0 * ((1-input$sigma)*do.call("cbind", replicate(input$M, Pnight, simplify = FALSE)) +
                                            input$sigma *do.call("rbind", replicate(input$M, Pday, simplify = FALSE)));
  
  FN <- GN - MN;
  FP <- GP - MP;
  
  #Now the forloop for the replicator equation
  for (i in 1:10000){
    FNmax <- max(max(FN)); FNmin = min(min(FN));
    FPmax <- max(max(FP)); FPmin = min(min(FP));
    fact <- 0.5 / max(c(FNmax, FPmax, -FNmin, -FPmin)); #max 20% increase per time step
    
    N <- N + fact* FN*N;
    P <- P + fact* FP*P;
    
    N[N < 10^-11] <- 10^-11; #to prevent extinction in the population
    P[P < 10^-11] <- 10^-11;
    
    N <- N / sum(sum(N)); P <- P/ sum(sum(P));
    
    Nday <- apply(N,2,sum);
    Nnight <- matrix(apply(N,1,sum), nrow = input$M, ncol = 1);
    Pday <- apply(P,2,sum);
    Pnight <- matrix(apply(P,1,sum), nrow = input$M, ncol = 1); 
    
    MN <- input$M * input$P0* ((1-input$sigma)* do.call("cbind", replicate(input$M, matrix(o, nrow = input$M, ncol = 1)*Pnight, simplify = FALSE)) +
                                 input$sigma*do.call("rbind", replicate(input$M, Pday*v, simplify = FALSE)));
    GP <- input$M * input$eta * input$N0 * ((1-input$sigma)*do.call("cbind", replicate(input$M, matrix(o, nrow = input$M, ncol = 1)*Nnight, simplify = FALSE)) +
                                              input$sigma * do.call("rbind", replicate(input$M, Nday*v, simplify = FALSE))) - CNP;
    MP <- input$mu * input$M * input$P0 * ((1-input$sigma)*do.call("cbind", replicate(input$M, Pnight, simplify = FALSE)) +
                                             input$sigma *do.call("rbind", replicate(input$M, Pday, simplify = FALSE)));
    
    FN <- GN - MN;
    FP <- GP - MP;
  }
  
  
par(mfrow = c(1,3))
plot(g,zi, type = "l", col = "green")
title("Resource distribution")

plot(Nday,zi, type = "l", col = "red")
lines(Nnight,zi, col = "blue")
title("Prey distribution")

plot(Pday,zi, type = "l", col = "red")
lines(Pnight,zi, col = "blue")
title("Predator distribution")

})
  }



# Run the application 
shinyApp(ui = ui, server = server)

