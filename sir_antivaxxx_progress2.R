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
    dSw <- - beta*(1-p)*Sw*Iw - beta*p*Sw*In #szuletes = halal, ezert a halalt parameterezem csak
    # d*sig+gi a total szuletesszam, ds a S-ek halala
    dSn <- - beta*(1-p)*Sn*In - beta*p*Sn*Iw
    # dS <- (d*(Sw+Iw+Rw)+g*Iw - beta*(1-p)*Sw*Iw - beta*p*Sw*In - d*Sw - fw*Sw) + 
    # (d*(Sn+In+Rn)+g*In - beta*(1-p)*Sn*In - beta*p*Sn*Iw - d*Sn - fn*Sn)
    dIw <- + beta*(1-p)*Sw*Iw + beta*p*Sw*In - r*Iw
    #dI a sima halalozasi rataval halo I-k, g*I amennyivel jobban meghalnak az I-k
    dIn <- + beta*(1-p)*Sn*In + beta*p*Sn*Iw - r*In
    # dI <- (+ beta*(1-p)*Sw*Iw + beta*p*Sw*In - r*Iw - d*Iw -g*Iw) +
    # (+ beta*(1-p)*Sn*In + beta*p*Sn*Iw - r*In - d*In -g*In)
    dRw <- + r*Iw        # Note: because S+I+R=constant, this equation could actually be omitted,
    # and R at any time point could simply be calculated as N-S-I.
    dRn <- + r*In
    # dR <- (r*Iw - d*Rw + fw*Sw) +
    # (r*In - d*Rn + fn*Sn)
    der <- c(dSn, dSw, 
             # dS , 
             dIn, dIw, 
             # dI, 
             dRn, dRw # , 
             # dR
    )
    list(der)  # the output must be returned    
  }) # end of 'with'

}  # end of function definition


###########################
# MAIN PROGRAM
###########################

### LOAD LIBRARIES
#load R library for ordinary differential equation solvers
library(deSolve)
R0_vektor = vector(mode="list") #ez lesz amiben az R0-okat tárolom
### INITIALIZE PARAMETER SETTINGS
for (w in 1:21) {        # most csak azt nezem, ahol az oltasellenesek aranya 0
  fn = 0.85 # 0.5 + (v-1)/10
  fw = 0.15 # (w-1)/10
  q = (w-1)/20 # Waldorfosok aránya a popuban
  p = 0.2 # a valószínűség, hogy akivel találkozol, nem olyan, mint te
  f = q*fw + (1-q)*fn
  index = f*100
  S0 = 499
  parms <- c(beta=1e-3, r=1e-1, fw=fw, fn=fn, q=q, p=p)		# set the parameters of the model
  inits <- c(Sn=(1-q)*(1-fn)*S0, 
             Sw=q*(1-fw)*S0, 
             # S=((1-q)*(1-fn)*S0) +  (q*(1-fw)*S0),
             In= (1-q),
             Iw= q,
             # I=1,
             Rn=(1-q)*fn*S0, 
             Rw=q*fw*S0 # ,
             # R=((1-q)*fn*S0) + (q*fw*S0)
    )
  print(q*(1-fw)*S0)
  # set the initial values
  dt    <- seq(0,500,0.1)			# set the time points for evaluation

  # Calculate and print R_0 on the screen
  N <- sum(inits)
  R_0 <- with(as.list(parms),{beta*(N-f*S0)/r})
#   print(paste("hany szazalek oltott", index))
#   print(paste("R_0 =",R_0),quote=FALSE)
#   R0_vektor[v] = (R_0) #a vektor/lista/kisfaszom v-edik eleme legyen az aktualis R0 ertek
#   names(R0_vektor)[v] = f*100 #a vektor/lista/kisfaszom v-edik elemenek a neve legyen az oltottsagi %


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

  par(mfrow=c(3,1))
  plot(dt, Sn, type="l", col="blue", 
       ylim=c(0,sum(inits)), 
       main = (1-q), xlab="time", ylab="number of individuals",lwd=3)
  lines(dt, In, type="l", col="red",lwd=3)
  lines(dt, Rn, type="l", col="darkgreen",lwd=3)
  
  # Add a legend to the graph
  legend(70,400, legend=c("Sn","In","Rn"), col=c("blue", "red", "darkgreen"), lty=1,lwd=2)
  #dev.off()
  
  plot(dt, Sw, type="l", col="blue", 
       ylim=c(0, sum(inits)), 
       main = q, xlab="time", ylab="number of individuals",lwd=3)
  lines(dt, Iw, type="l", col="red",lwd=3)
  lines(dt, Rw, type="l", col="darkgreen",lwd=3)
  
  # Add a legend to the graph
  legend(70,400, legend=c("Sw","Iw","Rw"), col=c("blue", "red", "darkgreen"), lty=1,lwd=2)
  #dev.off()
  
  plot(dt, (Sn+Sw), type="l", col="blue", 
       ylim=c(0, sum(inits)), 
       main = (q*fw+(1-q)*fn), xlab="time", ylab="number of individuals",lwd=3)
  lines(dt, (In+Iw), type="l", col="red",lwd=3)
  lines(dt, (Rn+Rw), type="l", col="darkgreen",lwd=3)

  # Add a legend to the graph
  legend(70,400, legend=c("S","I","R"), col=c("blue", "red", "darkgreen"), lty=1,lwd=2)
  #dev.off()

  detach(simulation) # clean up the search path
}
# print(R0_vektor)
# volte = FALSE #ebben tarolom, hogy van-e olyan oltottsag, ami megelozi a jarvanyt
# for (a in 1:length(R0_vektor)) {
#   if (R0_vektor[a] < 1) {
#     print(names(R0_vektor)[a])
#     volte = TRUE
#     break
#   }
# }
# if (volte == FALSE) {
#   print("Mind meg fogunk dogleni")
# }