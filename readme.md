## STRUCTURAL-GLASS

# Cosas que faltan:
- IMPLEMENTARE TRANSITION STATES USANDO L'ALGORITMO DI HEUER.
  È un upgrade al programma chunksect.
  Per ogni chunk ho la traiettoria.
  La modifica avviene in BisectChunk.py
  Il meccanismo di bisezione, ma ad ogni IS aggiungo un TS.

  1) Comincio da una IS, xi0=xi(t0), che viene da una configurazione
  termica x0=x(t0)

  2) Se la IS successiva è la stessa (definire "stessa"), non facciamo
  niente e passiamo a quella dopo. Se è diversa, cominciano le danze.

  3) Chiamo xi0 una IS t.c. la successiva è diversa.

  A questo punto lancio la funzione

  def FindTransitionState(xi0, xi1, x0, x1)
  

  4) Faccio un'interpolazione lineare tra x0 e x1. Attraverso una
  bisezione posso trovare x0' e x1', su quella linea, ma piu vicine
  tra loro, t.c. xi0'=xi0, mentre xi1'=xi1 (questo sarebbe da vedere).
  Questo sarebbe un'altra funzione
  def InterpolateConfs(x0,x1,alpha) trova l'interpolazione
  def BisectInterpolation(x0,x1) usa l'interpolazione per trovare una
  coppia migliore di configurazioni.

  5) Minimizzo sia x0' che x1' per pochi passi (guardare la risposta
  che mi hanno dato nel forum).

  6) Se la distanza tra le due e` maggiore di epsilon, rifaccio la bisezione.

  7) Itero varie volte 5 e 6.

  8) Minimizzo il gradiente fino a convergenza. Questo è il TS.

  9) Come debug mi calcolo la hessiana (deve avere solo un autovalore negativo).

- Una volta che so individuare la IS successiva, serializzare e
  individuare quella successiva a ogni IS.
- Chiedersi se la minimizzazione con FIRE potrebbe dare risultati
  differenti dalla steepest descent.
  FIRE, se impongo che alpha=1 tutto il tempo, diventa uno steepest descent.
  Come alternativa potrei minimizzare facendo NVT a temperatura nulla.
- Deteccion de cristalizacion con metodo de Garrahan y Chandler.
- Medias sobre samples de las Fkt, msd y tau para ver si el grupo de
  samples se puede considerar termalizado a pesar de que la muestra
  individual parezca que no.
- Una vez que esté todo listo y pasado al cluster, habrá que agregar
  los sets completos de parametros: todas las temperaturas, con
  tiempos de termalización, ecc...
- 

# Bugs que corregir:
Nada que señalar.

## Per sincronizzare questo repositorio
git remote add origin https://github.com/mbaityje/STRUCTURAL-GLASS.git
git push -u origin master

## Notas útiles
- k*T=0.00831445986144858*T en las unidades de la simulación


# Descrizione delle directory
- Mettere una descrizione di ognuno dei files presenti nelle cartelle di output
- Mettere una descrizione di ognuno dei programmi
