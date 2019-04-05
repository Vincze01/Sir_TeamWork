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

### INITIALIZE PARAMETER SETTINGS
for (v in 0:20) {
  f = v/20
  index = f*100
  S0 = 499
  parms <- c(beta=1e-3, r=1e-1)		# set the parameters of the model
  inits <- c(S=(1-f)*S0, I=1, R=f*S0)		# set the initial values
  dt    <- seq(0,100,0.1)			# set the time points for evaluation
  
  
  # Calculate and print R_0 on the screen
  N <- sum(inits)
  R_0 <- with(as.list(parms),{beta*(N-f*S0)/r})
  print(paste("hany szazalek oltott", index))
  print(paste("R_0 =",R_0),quote=FALSE)
  
  
  ### SIMULATE THE MODEL
  
  ## Use lsoda to solve the differential equations numerically. The syntax should be
  ## lsoda(initial values, time points, function, parameters)
  
  simulation <- as.data.frame(lsoda(inits, dt, SIR, parms=parms)) # this way our set 'parms' will be used as default
  
  ### PLOT THE OUTPUT
  
  # If you remove the # before pdf(...) and dev.off(), the output will be written in a pdf file,
  # in the working directory. If you don't, a window containing your graph will just pop up.
  
  #pdf("startingscript.pdf")
  #par(cex=1.7)
  # Plot S according to time, in blue, and add the graph I and R according to time,
  # in red and dark green respectively. Call help(plot) for further details.
  
  attach(simulation) # this command allows you to refer to the columns of the data frame directly.
  library(ggplot2)
  
  title = paste("R_0 = ",R_0, sep = "")
  beta = parms[1]
  r = parms[2]
  beta2 = paste("beta = ",beta, sep = "")
  r2 = paste("r = ",r,sep = "")
  percent = paste(index, "% oltott")
  subtitle = paste(percent,", ","\n",beta2,", ",r2)
  
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
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "antiquewhite"))+
    expand_limits(y=c(0,500))
  print(plot)
  
  detach(simulation) # clean up the search path
}