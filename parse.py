import os
import pandas as pd
import numpy as np

euristiche = []

values_names = ["Media passi","Media fissati","Media max nodi fissati",\
        "Media nodi no fissati","Max passi","Max fissati","Max max nodi fissati",\
        "Min nodi no fissati"]

heur = ["aux","auxMaxPeso","auxMinPeso","maxArchiDaFissare","maxPeso",\
        "minArchiDaFissare","minPeso","noAuxNoThorn","noAuxNoThornMaxPeso",\
        "noAuxNoThornMinPeso","thorn","thornMaxPeso","thornMinPeso"]

############# FUNCTIONS ##################

def writeToCSV(values):

    df = pd.DataFrame(data=values, columns=values_names)
    df.insert(0,'Euristica',heur)
    df.to_csv("tabs.csv",sep=",", index=False)



a = []
for filename in sorted(os.listdir('summaries/')):
    text = open('summaries/' + filename, 'r').read()
    text = text.splitlines()
    toKeep = [0, 1, 2, 4, 12, 13, 14, -2]
    toExport = []
    for i in toKeep:
        toExport.append(text[i].split(': ')[1])
    a.append(toExport)

writeToCSV(a)
