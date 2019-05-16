##ZSÓFI's ggplot graphs
library(ggplot2)

#itt írom be, hogy kiírja nekem a beta-kat meg r-eket is a plotra
#ITT ADOM MEG A TITLE-T!!! meg a subtitle-t
title = paste("R_0 = ",R_0, sep = "")
beta = parms[1]
r = parms[2]
beta2 = paste("beta = ",beta, sep = "")
r2 = paste("r = ",r,sep = "")
subtitle = paste(beta2,r2, sep = ", ")
#^^Tehát ez a rész itt csak annyi volt, hogy megadtam neki milyen szöveg mellé mely értékeket kell irogatnia
#a for ciklusban. Pl hogy a beta az a parms elso parmétere, amit úgy írjon ki nekem hogy "beta = az értéke" és
#hogy ezt hogy válassza el. Vagyishogy nálam sehogy: ""
#Aztán ugyanezt erre is megcsináltam
#Majd csináltam nekik subtitle változót, amit meg vesszovel válasszon el

#plot script
plot = ggplot(simulation, aes(x = time)) + #simulation a data ugye
  geom_line(aes(y = S, colour = "Susceptible"), size=3, alpha=0.45) + #egy vonal behúzása S-re, a színnél meg
  #azt adom meg, amit alul definiálok a colorban, ami a legendet csinálja, vagyis Susceptible. 
  #Amit ott beírsz a szín mellé, azt a szöveget jeleníti meg a szín mellé
  #a size az a vonal vastagsága, az alpha meg, hogy mennyire legyen átlátszó
  geom_line(aes(y = I, colour = "Infected"), size=3.1, alpha=0.45) + 
  geom_line(aes(y = R, colour = "Recovered"), size=3.2, alpha=0.45) +
  ylab(label="Number of individuals") + #x tengely felirata
  xlab("Time") + #y tengely felirata
  ggtitle(title, subtitle = subtitle) + #ebben mondod meg, hogy mi legyen a title meg a subtitle
  #de ha ugyanazokat az elnevezéseket használod, mint én és csak ott átirogatod a dolgokat, akkor ehhez nem kell
  #hozzányúlnod
  scale_color_manual(name = " ",
                     values = c("Susceptible" = "turquoise4", "Infected" = "tomato3", "Recovered" = "goldenrod2")) +
  #itt add meg, hogy melyik vonalad milyen színu legyen, és hogy hívják azt a változót. De szerintem ezt is meg-
  #tarthatod, hogy követheto legyen és mindig ugyanolyan legyen mindennek a színe
  #INNENTOL LESZAROD A KÓDOT!!! NEKED NEM KELL!! EZT ÍGY HAGYOD!!!! itt mindenféle háttérszínezést meg
  #formázást adok meg, amitol csillámpóni szép lesz az egész, de ez mindenhol UGYANAZ!!!
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