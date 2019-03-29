##########################################################
# Starting script to the module 'SIR models of epidemics'
##########################################################

# implements the basic SIR model, and plots simulation results

## Here you only need to know basic things:
##  #	To make comments
##  <-	Assigning a value to a symbol
##  c	To create a vector c(v1, v2, ..., vn)
##  e	is *10^
##  seq(a,b,c)	To create a sequence from a to b with intervals c


###################################
# FUNCTION DEFINITIONS
###################################

###
# SIR <- function(t, x, parms)
# Use: calculates the derivatives for the SIR model
# Input: 
#      t: time (not used here, because there is no explicit time dependence)
#      x: vector of the current values of all variables (S, I, R)
#      parms: vector of all parameter values
# Output:
#      der: vector of derivatives

# To use the lsoda function, the model function has to be a function of t (time),
# x (the vector of the variables) and parms (the parameters of the model).


SIR <- function(t, x, parms){
  
  # Beta and r are not global variables. This means that if you type beta in R, the output will be 'beta',
  # and not its value. You have to specify that you want to use the value of beta from 'parms' to solve the ODEs.
  # Similarly, the variables of the model are taken from the vector x. This is done by the 'with' function.
  
  with(as.list(c(parms,x)),{
    dS <- - beta*S*I
    dI <- + beta*S*I - r*I
    dR <- r*I         # Note: because S+I+R=constant, this equation could actually be omitted,
    # and R at any time point could simply be calculated as N-S-I.
    der <- c(dS, dI,dR)
    list(der)  # the output must be returned    
  }) # end of 'with'
  
}  # end of function definition


###########################
# MAIN PROGRAM
###########################

### LOAD LIBRARIES
#load R library for ordinary differential equation solvers
library(deSolve)

#For loop to try different Beta and r values
index = 0
for(i in 1:5) {
  for(j in 1:5) {
    
    index = index + 1
    
    ### INITIALIZE PARAMETER SETTINGS
    
    parms <- c(beta=10^(-i), r=10^(-j+2))		# set the parameters of the model - per time unit!
    inits <- c(S=499, I=1, R=0)		# set the initial values
    dt    <- seq(0,100,0.1)			# set the time points for evaluation
    
    # Calculate and print R_0 on the screen
    N <- sum(inits)
    R_0 <- with(as.list(parms),{beta*N/r})
    
    
    
    ### SIMULATE THE MODEL
    
    ## Use lsoda to solve the differential equations numerically. The syntax should be
    ## lsoda(initial values, time points, function, parameters)
    
    
    
    simulation <- as.data.frame(lsoda(inits, dt, SIR, parms = parms)) # this way our set 'parms' will be used as default
    
    ### PLOT THE OUTPUT
    ##ZSÓFI's ggplot graphs
    library(ggplot2)
    
    #itt írom be, hogy kiírja nekem a beta-kat meg r-eket is a plotra
    title = paste("R_0 = ",R_0, sep = "")
    beta = parms[1]
    r = parms[2]
    beta2 = paste("beta = ",beta, sep = "")
    r2 = paste("r = ",r,sep = "")
    subtitle = paste(beta2,r2, sep = ", ")
    
    #plot script
    plot = ggplot(simulation, aes(x = time)) + 
      geom_line(aes(y = S, colour = "Susceptible"), size=3, alpha=0.45) + 
      geom_line(aes(y = I, colour = "Infected"), size=3.1, alpha=0.45) + 
      geom_line(aes(y = R, colour = "Recovered"), size=3.2, alpha=0.45) +
      ylab(label="Number of individuals") + 
      xlab("Time") +
      ggtitle(title, subtitle = subtitle) +
      scale_color_manual(name = " ",
                         values = c("Susceptible" = "turquoise4", "Infected" = "tomato3", "Recovered" = "goldenrod2")) +
      theme(plot.background = element_rect(fill = "antiquewhite4"),
            legend.background = element_rect(fill = "antiquewhite4"),
            plot.title = element_text(size = 18, lineheight=.8, hjust=0.5, face="bold", colour="antiquewhite"),
            plot.subtitle = element_text(size = 16, lineheight=.8, hjust=0.5, face="italic", colour="antiquewhite"),
            axis.title = element_text(size = 15, face="bold",colour="antiquewhite"),
            axis.text = element_text(size = 12,colour="antiquewhite"),
            legend.position="top",
            legend.key = element_rect(fill = "seashell"),
            legend.text = element_text(size = 12,colour="antiquewhite"),
            legend.box.background = element_rect(colour = "antiquewhite", size=1.5),
            panel.background = element_rect(fill = "seashell", colour = "seashell", size = 0.5, linetype = "solid"),
            panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "antiquewhite"), 
            panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "antiquewhite"))
    print(plot)
    
    attach(simulation) # this command allows you to refer to the columns of the data frame directly.
    
    #print(c(toString(5 * i + j), max(I)))
    
    print(paste("index: ", index))
    
    initI = unname( inits[2])
    if(initI < max(I))
    {
      print("Járvány")
    }
    else
    {
      print("Nem járvány")
    }
    
    print(paste("R_0 =",R_0),quote=FALSE)
    print(unname( inits[2]))
    print(max(I))
    print("==================================")
    
    detach(simulation) # clean up the search path
    
  }
}
