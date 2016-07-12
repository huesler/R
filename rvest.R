
library(rvest)
library(stringi)
library(stringr)





##########################################
#Iteration ueber die Seiten
for (i in 1:10) {
##########################################

  
##########################################
#Beachtung des Umstandes, dass erste Seite andere URL hat
if (i == 1) {
  adresse <- "https://uefaeurotippspiel.srf.ch/app/rangliste.jsp?typ=gesamt&gruppe=SPIELER"
  } else {
  adresse <-   paste("https://uefaeurotippspiel.srf.ch/app/rangliste.jsp?umgebung=false&gruppe=SPIELER&tippgemeinschaft=&typ=gesamt&wettrunde=&seite=",i, sep="")
  }
##########################################


##########################################
html   <- read_html(adresse)
rang   <- html_nodes(html, "td.first.position")
name   <- html_nodes(html, ".spieler-name span")
punkte <- html_nodes(html, ".typ-gesamt td.punkte")

rang     <- as.numeric(str_sub(str_trim(html_text(rang)), 1, -2))
name     <- html_text(name)
punkte   <- as.numeric(html_text(punkte))
##########################################



##########################################
#Encoding-Reparatur kann auch am Schluss gemacht werden
#name     <- repair_encoding(name)
##########################################

df <- data.frame(rang, name, punkte)

##########################################
#Beim ersten Loop wird Tabelle "gesamt" erstellt. Danach "df" angehaengt
if (i == 1) {
  gesamt <- df
  } else {
    gesamt <- rbind(gesamt, df)
    }
##########################################



##########################################
#Iteration ueber die Seiten
}
##########################################






