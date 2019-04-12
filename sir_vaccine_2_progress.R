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
R0_vektor = vector(mode="list")
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
  R0_vektor[v] = (R_0)
  names(R0_vektor[v]) = f*100


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

  plot(dt, S, type="l", col="blue", ylim=c(0,sum(inits)), main = R_0, xlab="time", ylab="number of individuals",lwd=3)
  lines(dt, I, type="l", col="red",lwd=3)
  lines(dt, R, type="l", col="darkgreen",lwd=3)

  # Add a legend to the graph
  legend(70,400, legend=c("S","I","R"), col=c("blue", "red", "darkgreen"), lty=1,lwd=2)
  #dev.off()

  detach(simulation) # clean up the search path
}
print(R0_vektor)