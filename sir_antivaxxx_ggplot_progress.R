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
# oltottsági ráta = 1/időegység, teljes oltottsági ráta = oltottsági arány*1/időegység
# felezési idő

SIR <- function(t, x, parms){

# Beta and r are not global variables. This means that if you type beta in R, the output will be 'beta',
# and not its value. You have to specify that you want to use the value of beta from 'parms' to solve the ODEs.
# Similarly, the variables of the model are taken from the vector x. This is done by the 'with' function.

  with(as.list(c(parms,x)),{
    dSw <- d*(Sw+Iw+Rw)+g*Iw - beta*Sw*Iw - beta*p*Sw*In - d*Sw - ffw*Sw
   
    dSn <- d*(Sn+In+Rn)+g*In - beta*Sn*In - beta*p*Sn*Iw - d*Sn - ffn*Sn

    dIw <- + beta*Sw*Iw + beta*p*Sw*In - r*Iw - d*Iw -g*Iw

    dIn <- + beta*Sn*In + beta*p*Sn*Iw - r*In - d*In -g*In

    dRw <- r*Iw - d*Rw + ffw*Sw        # Note: because S+I+R=constant, this equation could actually be omitted,
                  # and R at any time point could simply be calculated as N-S-I.
    dRn <- r*In - d*Rn + ffn*Sn

    der <- c(dSn, dSw, 

             dIn, dIw, 

             dRn, dRw
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
library(ggplot2)
# pdf("sir_antivaxx_birthdeath_gg.pdf")
### INITIALIZE PARAMETER SETTINGS
for (w in 1:21) {
  for(j in 1:11){
    fn = 0.85
    fw = 0.85
    q = (w-1)/20 # Waldorfosok aránya a popuban
    p = (j-1)/10 # a valószínűség, hogy akivel találkozol, nem olyan, mint te
    f = q*fw + (1-q)*fn
    index = f*100
    S0 = 499
    parms <- c(beta=1e-2, r=1e-1, d=0.1, g=0.5, ffw=0, ffn=fn, q=q, p=p)		# set the parameters of the model
    inits <- c(Sn=(1-q)*(1-fn)*S0, 
               Sw=q*(1-fw)*S0, 
               In=(1-q),
               Iw=q,
               Rn=(1-q)*fn*S0, 
               Rw=q*fw*S0
               )
    
    # set the initial values
    dt    <- seq(0,50,0.1)			# set the time points for evaluation
  
  
    # Calculate and print R_0 on the screen
    N <- sum(inits)
    # R_0 <- with(as.list(parms),{beta*(N-f*S0)/(r+d+g)})
  
    ### SIMULATE THE MODEL
  
    ## Use lsoda to solve the differential equations numerically. The syntax should be
    ## lsoda(initial values, time points, function, parameters)
  
    simulation <- as.data.frame(lsoda(inits, dt, SIR, parms=parms)) # this way our set 'parms' will be used as default
  
    ### PLOT THE OUTPUT
  
    # If you remove the # before pdf(...) and dev.off(), the output will be written in a pdf file,
    # in the working directory. If you don't, a window containing your graph will just pop up.
  
    
    # par(cex=1.7)
    # Plot S according to time, in blue, and add the graph I and R according to time,
    # in red and dark green respectively. Call help(plot) for further details.
  
    attach(simulation) # this command allows you to refer to the columns of the data frame directly.
    
    #itt ?rom be, hogy ki?rja nekem a beta-kat meg r-eket is a plotra
    #ITT ADOM MEG A TITLE-T!!! meg a subtitle-t
    title = paste("Normálisak")
    ffn = parms[6]
    arany = paste("arány: ",(1-q), sep = "")
    kev = paste("keveredés: ",p,sep = "")
    oltotts = paste("eredeti oltottság: ", fn, sep = "")
    rata = paste("oltási ráta: ", ffn, sep = "")
    subtitle = paste(arany,kev,oltotts, rata, sep = ", ")
    #^^Teh?t ez a r?sz itt csak annyi volt, hogy megadtam neki milyen sz?veg mell? mely ?rt?keket kell irogatnia
    #a for ciklusban. Pl hogy a beta az a parms elso parm?tere, amit ?gy ?rjon ki nekem hogy "beta = az ?rt?ke" ?s
    #hogy ezt hogy v?lassza el. Vagyishogy n?lam sehogy: ""
    #Azt?n ugyanezt erre is megcsin?ltam
    #Majd csin?ltam nekik subtitle v?ltoz?t, amit meg vesszovel v?lasszon el
    
    #plot script
    p_norm = ggplot(simulation, aes(x = time)) + #simulation a data ugye
      geom_line(aes(y = Sn, colour = "Susceptible"), size=3, alpha=0.45) + #egy vonal beh?z?sa S-re, a sz?nn?l meg
      #azt adom meg, amit alul defini?lok a colorban, ami a legendet csin?lja, vagyis Susceptible. 
      #Amit ott be?rsz a sz?n mell?, azt a sz?veget jelen?ti meg a sz?n mell?
      #a size az a vonal vastags?ga, az alpha meg, hogy mennyire legyen ?tl?tsz?
      geom_line(aes(y = In, colour = "Infected"), size=3.1, alpha=0.45) + 
      geom_line(aes(y = Rn, colour = "Recovered"), size=3.2, alpha=0.45) +
      # ylim(0, NA) +
      ylab(label="Number of individuals") + #x tengely felirata
      xlab("Time") + #y tengely felirata
      ggtitle(title, subtitle = subtitle) + #ebben mondod meg, hogy mi legyen a title meg a subtitle
      #de ha ugyanazokat az elnevez?seket haszn?lod, mint ?n ?s csak ott ?tirogatod a dolgokat, akkor ehhez nem kell
      #hozz?ny?lnod
      scale_color_manual(name = " ",
                         values = c("Susceptible" = "turquoise4", "Infected" = "tomato3", "Recovered" = "goldenrod2")) +
      #itt add meg, hogy melyik vonalad milyen sz?nu legyen, ?s hogy h?vj?k azt a v?ltoz?t. De szerintem ezt is meg-
      #tarthatod, hogy k?vetheto legyen ?s mindig ugyanolyan legyen mindennek a sz?ne
      #INNENTOL LESZAROD A K?DOT!!! NEKED NEM KELL!! EZT ?GY HAGYOD!!!! itt mindenf?le h?tt?rsz?nez?st meg
      #form?z?st adok meg, amitol csill?mp?ni sz?p lesz az eg?sz, de ez mindenhol UGYANAZ!!!
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
    print(p_norm)
    
    #itt ?rom be, hogy ki?rja nekem a beta-kat meg r-eket is a plotra
    #ITT ADOM MEG A TITLE-T!!! meg a subtitle-t
    title = paste("Antivaxxerek")
    ffw = parms[5]
    arany = paste("arány: ",q, sep = "")
    kev = paste("keveredés: ",p,sep = "")
    oltotts = paste("eredeti oltottság: ", fw, sep = "")
    rata = paste("oltási ráta: ", ffw, sep = "")
    subtitle = paste(arany,kev,oltotts, rata, sep = ", ")
    #^^Teh?t ez a r?sz itt csak annyi volt, hogy megadtam neki milyen sz?veg mell? mely ?rt?keket kell irogatnia
    #a for ciklusban. Pl hogy a beta az a parms elso parm?tere, amit ?gy ?rjon ki nekem hogy "beta = az ?rt?ke" ?s
    #hogy ezt hogy v?lassza el. Vagyishogy n?lam sehogy: ""
    #Azt?n ugyanezt erre is megcsin?ltam
    #Majd csin?ltam nekik subtitle v?ltoz?t, amit meg vesszovel v?lasszon el
    
    #plot script
    p_wald = ggplot(simulation, aes(x = time)) + #simulation a data ugye
      geom_line(aes(y = Sw, colour = "Susceptible"), size=3, alpha=0.45) + #egy vonal beh?z?sa S-re, a sz?nn?l meg
      #azt adom meg, amit alul defini?lok a colorban, ami a legendet csin?lja, vagyis Susceptible. 
      #Amit ott be?rsz a sz?n mell?, azt a sz?veget jelen?ti meg a sz?n mell?
      #a size az a vonal vastags?ga, az alpha meg, hogy mennyire legyen ?tl?tsz?
      geom_line(aes(y = Iw, colour = "Infected"), size=3.1, alpha=0.45) + 
      geom_line(aes(y = Rw, colour = "Recovered"), size=3.2, alpha=0.45) +
      # ylim(0, NA) +
      ylab(label="Number of individuals") + #x tengely felirata
      xlab("Time") + #y tengely felirata
      ggtitle(title, subtitle = subtitle) + #ebben mondod meg, hogy mi legyen a title meg a subtitle
      #de ha ugyanazokat az elnevez?seket haszn?lod, mint ?n ?s csak ott ?tirogatod a dolgokat, akkor ehhez nem kell
      #hozz?ny?lnod
      scale_color_manual(name = " ",
                         values = c("Susceptible" = "turquoise4", "Infected" = "tomato3", "Recovered" = "goldenrod2")) +
      #itt add meg, hogy melyik vonalad milyen sz?nu legyen, ?s hogy h?vj?k azt a v?ltoz?t. De szerintem ezt is meg-
      #tarthatod, hogy k?vetheto legyen ?s mindig ugyanolyan legyen mindennek a sz?ne
      #INNENTOL LESZAROD A K?DOT!!! NEKED NEM KELL!! EZT ?GY HAGYOD!!!! itt mindenf?le h?tt?rsz?nez?st meg
      #form?z?st adok meg, amitol csill?mp?ni sz?p lesz az eg?sz, de ez mindenhol UGYANAZ!!!
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
    print(p_wald)
    
    #itt ?rom be, hogy ki?rja nekem a beta-kat meg r-eket is a plotra
    #ITT ADOM MEG A TITLE-T!!! meg a subtitle-t
    title = paste("Teljes populáció")
    arany = paste("Antiwaxxerek aránya arány: ",q, sep = "")
    kev = paste("keveredés: ",p,sep = "")
    oltotts = paste("eredeti oltottság: ", fw, sep = "")
    subtitle = paste(arany,kev,oltotts, sep = ", ")
    #^^Teh?t ez a r?sz itt csak annyi volt, hogy megadtam neki milyen sz?veg mell? mely ?rt?keket kell irogatnia
    #a for ciklusban. Pl hogy a beta az a parms elso parm?tere, amit ?gy ?rjon ki nekem hogy "beta = az ?rt?ke" ?s
    #hogy ezt hogy v?lassza el. Vagyishogy n?lam sehogy: ""
    #Azt?n ugyanezt erre is megcsin?ltam
    #Majd csin?ltam nekik subtitle v?ltoz?t, amit meg vesszovel v?lasszon el
    
    #plot script
    p_popu = ggplot(simulation, aes(x = time)) + #simulation a data ugye
      geom_line(aes(y = (Sn+Sw), colour = "Susceptible"), size=3, alpha=0.45) + #egy vonal beh?z?sa S-re, a sz?nn?l meg
      #azt adom meg, amit alul defini?lok a colorban, ami a legendet csin?lja, vagyis Susceptible. 
      #Amit ott be?rsz a sz?n mell?, azt a sz?veget jelen?ti meg a sz?n mell?
      #a size az a vonal vastags?ga, az alpha meg, hogy mennyire legyen ?tl?tsz?
      geom_line(aes(y = (In+Iw), colour = "Infected"), size=3.1, alpha=0.45) + 
      geom_line(aes(y = (Rn+Rw), colour = "Recovered"), size=3.2, alpha=0.45) +
      # ylim(0, NA) +
      ylab(label="Number of individuals") + #x tengely felirata
      xlab("Time") + #y tengely felirata
      ggtitle(title, subtitle = subtitle) + #ebben mondod meg, hogy mi legyen a title meg a subtitle
      #de ha ugyanazokat az elnevez?seket haszn?lod, mint ?n ?s csak ott ?tirogatod a dolgokat, akkor ehhez nem kell
      #hozz?ny?lnod
      scale_color_manual(name = " ",
                         values = c("Susceptible" = "turquoise4", "Infected" = "tomato3", "Recovered" = "goldenrod2")) +
      #itt add meg, hogy melyik vonalad milyen sz?nu legyen, ?s hogy h?vj?k azt a v?ltoz?t. De szerintem ezt is meg-
      #tarthatod, hogy k?vetheto legyen ?s mindig ugyanolyan legyen mindennek a sz?ne
      #INNENTOL LESZAROD A K?DOT!!! NEKED NEM KELL!! EZT ?GY HAGYOD!!!! itt mindenf?le h?tt?rsz?nez?st meg
      #form?z?st adok meg, amitol csill?mp?ni sz?p lesz az eg?sz, de ez mindenhol UGYANAZ!!!
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
    print(p_popu)
    
    
   #  par(mfrow=c(3,1))
   #  plot(dt, Sn, type="l", col="blue", 
   #       ylim=c(0,sum(inits)), 
   #       main = (1-q), xlab="time", ylab="number of individuals",lwd=3)
   #  lines(dt, In, type="l", col="red",lwd=3)
   #  lines(dt, Rn, type="l", col="darkgreen",lwd=3)
   #  
   #  # Add a legend to the graph
   #  legend(70,400, legend=c("Sn","In","Rn"), col=c("blue", "red", "darkgreen"), lty=1,lwd=2)
   #  #dev.off()
   #  
   #  plot(dt, Sw, type="l", col="blue", 
   #       ylim=c(0, sum(inits)), 
   #       main = q, xlab="time", ylab="number of individuals",lwd=3)
   #  lines(dt, Iw, type="l", col="red",lwd=3)
   #  lines(dt, Rw, type="l", col="darkgreen",lwd=3)
   #  
   #  # Add a legend to the graph
   #  legend(70,400, legend=c("Sw","Iw","Rw"), col=c("blue", "red", "darkgreen"), lty=1,lwd=2)
   #  #dev.off()
   #  
   #  plot(dt, (Sn+Sw), type="l", col="blue", 
   #       ylim=c(0, sum(inits)), 
   #       main = (p), sub=p, xlab="time", ylab="number of individuals",lwd=3)
   #  lines(dt, (In+Iw), type="l", col="red",lwd=3)
   #  lines(dt, (Rn+Rw), type="l", col="darkgreen",lwd=3)
   #  
   #  # Add a legend to the graph
   #  legend(70,400, legend=c("S","I","R"), col=c("blue", "red", "darkgreen"), lty=1,lwd=2)
   #  
  # 
   detach(simulation) # clean up the search path
  }  
}
# dev.off()
