## STRUCTURAL-GLASS

# Cosas que faltan:
- Tutorial che calcola: 
-- gradiente quadro G2 senza passare per il gradiente.
-- derivata di G2
-- Minimizza l'energia con steepest descent
-- Minimizza G2 con steepest descent
-- Calcola l'Hessiana dell'energia e i suoi autovalori

- Completare algoritmo di ricerca Transition States in BisectChunk.py
-- Verificare l'ordine della sella trovata, che dovrebbe essere 1.

- Confrontare Bisect FIRE con Minimizzazione steepest descent

- Deteccion de cristalizacion con metodo de Garrahan y Chandler.

- Medias sobre samples de las Fkt, msd y tau para ver si el grupo de
  samples se puede considerar termalizado a pesar de que la muestra
  individual parezca que no.


# Bugs que corregir:
Il programma SelfIntermediateScatteringFunctions viene lanciato da dentro la cartella di output. Visto che questa viene raggiunta attraverso un symlink, la ricostruzione successiva delle cartelle viene fuori sbagliata.

In quello stesso programma trajFreq=1 il che risulta in traiettorie enormi con troppi time steps.

I jobs su talapas non stanno correndo sulle GPU.


## Per sincronizzare questo repositorio
git remote add origin https://github.com/mbaityje/STRUCTURAL-GLASS.git
git push -u origin master

## Notas útiles
- k*T=0.00831445986144858*T en las unidades de la simulación


# Descrizione delle directory
- Mettere una descrizione di ognuno dei files presenti nelle cartelle di output
- Mettere una descrizione di ognuno dei programmi
