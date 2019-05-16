##ZS�FI's ggplot graphs
library(ggplot2)

#itt �rom be, hogy ki�rja nekem a beta-kat meg r-eket is a plotra
#ITT ADOM MEG A TITLE-T!!! meg a subtitle-t
title = paste("R_0 = ",R_0, sep = "")
beta = parms[1]
r = parms[2]
beta2 = paste("beta = ",beta, sep = "")
r2 = paste("r = ",r,sep = "")
subtitle = paste(beta2,r2, sep = ", ")
#^^Teh�t ez a r�sz itt csak annyi volt, hogy megadtam neki milyen sz�veg mell� mely �rt�keket kell irogatnia
#a for ciklusban. Pl hogy a beta az a parms elso parm�tere, amit �gy �rjon ki nekem hogy "beta = az �rt�ke" �s
#hogy ezt hogy v�lassza el. Vagyishogy n�lam sehogy: ""
#Azt�n ugyanezt erre is megcsin�ltam
#Majd csin�ltam nekik subtitle v�ltoz�t, amit meg vesszovel v�lasszon el

#plot script
plot = ggplot(simulation, aes(x = time)) + #simulation a data ugye
  geom_line(aes(y = S, colour = "Susceptible"), size=3, alpha=0.45) + #egy vonal beh�z�sa S-re, a sz�nn�l meg
  #azt adom meg, amit alul defini�lok a colorban, ami a legendet csin�lja, vagyis Susceptible. 
  #Amit ott be�rsz a sz�n mell�, azt a sz�veget jelen�ti meg a sz�n mell�
  #a size az a vonal vastags�ga, az alpha meg, hogy mennyire legyen �tl�tsz�
  geom_line(aes(y = I, colour = "Infected"), size=3.1, alpha=0.45) + 
  geom_line(aes(y = R, colour = "Recovered"), size=3.2, alpha=0.45) +
  ylab(label="Number of individuals") + #x tengely felirata
  xlab("Time") + #y tengely felirata
  ggtitle(title, subtitle = subtitle) + #ebben mondod meg, hogy mi legyen a title meg a subtitle
  #de ha ugyanazokat az elnevez�seket haszn�lod, mint �n �s csak ott �tirogatod a dolgokat, akkor ehhez nem kell
  #hozz�ny�lnod
  scale_color_manual(name = " ",
                     values = c("Susceptible" = "turquoise4", "Infected" = "tomato3", "Recovered" = "goldenrod2")) +
  #itt add meg, hogy melyik vonalad milyen sz�nu legyen, �s hogy h�vj�k azt a v�ltoz�t. De szerintem ezt is meg-
  #tarthatod, hogy k�vetheto legyen �s mindig ugyanolyan legyen mindennek a sz�ne
  #INNENTOL LESZAROD A K�DOT!!! NEKED NEM KELL!! EZT �GY HAGYOD!!!! itt mindenf�le h�tt�rsz�nez�st meg
  #form�z�st adok meg, amitol csill�mp�ni sz�p lesz az eg�sz, de ez mindenhol UGYANAZ!!!
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