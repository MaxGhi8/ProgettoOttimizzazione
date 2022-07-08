
from math import sqrt, ceil
import numpy as np
# Import the NetworkX library per fare i grafi
import networkx as nx

from pyomo.environ import ConcreteModel, Var, Objective, SolverFactory
from pyomo.environ import minimize, Binary, RangeSet, ConstraintList

import matplotlib.pyplot as plt

# funzione per leggere il file
def ParseFile(filename):
    # OUTPUT: Ls lista con i nomi delle squadre e le due coordinate 
    doc = open(filename, 'r', encoding = 'utf-8')
    # for _ in range(50):
    doc.readline() # Salto la prima linea con le intestazioni
    # Leggo i circoli e faccio una lista di liste, ogni sottolista contiene
    # ordinatamente nome squadra come stringa e le due coordinate come float
    Ls = []
    for row in doc:
        row = row.split(',')
        # print(row)
        Ls.append( [row[1], float(row[2]), float(row[3])] )
    return Ls 


def CostMatrix(Ls):
    # INPUT: Ls è la lista in output da ParseFile
    # OUTPUT: Cm matrice con i costi Cm(i,j) = costo da squadra i a squadra j
    n = len(Ls)
    Cm = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            costo = (Ls[i][1] - Ls[j][1])**2 + (Ls[i][2] - Ls[j][2])**2 # costi al quadrato
            # Delta_phi = (Ls[i][1] - Ls[j][1])/2
            # costo = 2*6360*np.arcsin(sqrt( np.sin(Delta_phi)**2 +
            #                          ( 1-np.sin(Delta_phi)**2 - np.sin((Ls[i][1] + Ls[j][1])/2 )**2 )*
            #                          np.sin((Ls[i][2] - Ls[j][2])/2)**2))
            Cm[i,j] = costo
    return Cm
   
    
def Tau(n, M):
    # Costruzione di tau
    tau = list(M for _ in range(ceil(n/M)))
    resto = sum(tau) - n
    for i in range(resto):
        tau[i] = tau[i] - 1
    return tau
    
    
def ILP(M, Ls):
    # INPUT: - M cardinalità gironi (?ACTUNG: M divide o no?)
    #        - Ls lista delle squadre
    
    n = len(Ls) #numero di dati (=numero di squadre)
    Cm = CostMatrix(Ls)
    tau = Tau(n, M)
    kappa = len(tau) # numero di cluster (= gironi)

    # Build ILP Model
    model = ConcreteModel()
    
    # Indici
    model.N = RangeSet(n)
    model.K = RangeSet(kappa)

    # Variabili
    model.y = Var(model.N, model.K, domain = Binary)
    model.z = Var(model.N, model.N, model.K, domain = Binary)

    # Objective function of arc variables:
    # Distanze al quadrato
    # model.obj = Objective(expr = sum(Cm[i-1,j-1]*model.z[i,j,k] for i in model.N for j in model.N for k in model.K), 
    #                       sense = minimize)
    
    # Distanze
    model.obj = Objective(expr = sum(sqrt(Cm[i-1,j-1])*model.z[i,j,k] for i in model.N for j in model.N for k in model.K), 
                          sense = minimize)

    # Vincoli (y_{i,k} = 1, y_{j,k} = 1) <--> (z_{i,j,k} = 1)
    model.iff = ConstraintList()            
    for i in model.N:
        for j in model.N:
            for k in model.K:
                model.iff.add(expr = model.z[i,j,k] >= model.y[i,k] + model.y[j,k] - 1)
                model.iff.add(expr = model.z[i,j,k] <= model.y[i,k])
                model.iff.add(expr = model.z[i,j,k] <= model.y[j,k])
       
    # Vincolo sulla cardinalità dei gironi
    model.gironi = ConstraintList() 
    for k in model.K:
        model.gironi.add(expr = sum(model.y[i,k] for i in model.N) == tau[k-1])
        
    # Ogni squadra deve comparire in un solo girone
    model.unico = ConstraintList() 
    for i in model.N:
        model.unico.add(expr = sum(model.y[i,k] for k in model.K) == 1)
    
    # Risolviamo il ILP con gurobi
    solver = SolverFactory('gurobi')    
    sol = solver.solve(model, tee=False)

    # check feasibility
    sol_json = sol.json_repn()
    # Check solution status
    if sol_json['Solver'][0]['Status'] != 'ok':
        return None    
    
    # print per debugging
    # for i in model.N:
    #     for j in model.N:
    #         for k in model.K:
    #             if model.z[i,j,k]() > 0.5:
    #                 print('da {} a {} in girone {} --> {}\n'.format(i, j, k, model.z[i,j,k]()))
    
    costo_totale = model.obj()/2
    
    dict_gironi = {}
    for k in model.K:
        girone_k = []
        for i in model.N:
            if model.y[i,k]() > 0.5:
                girone_k.append( (Ls[i-1][0], i) )
        dict_gironi['Girone' + str(k)] = girone_k
        
    return dict_gironi, costo_totale
    
def PlotSolution(dict_gironi, lista_coord):
    for girone in dict_gironi:
        squadre_girone = dict_gironi[girone]
        plt.scatter( list(map(lambda x: lista_coord[x[1]-1][0], squadre_girone)),
                    list(map(lambda x: lista_coord[x[1]-1][1], squadre_girone)) )
    plt.show()
    
def CalcolaCosto(dict_gironi, lista_coord):
    costo = 0
    for girone in dict_gironi:
        squadre_girone = dict_gironi[girone]
        for squadra1 in squadre_girone:
            for squadra2 in squadre_girone:
                distanza = sqrt((lista_coord[squadra1[1]-1][0] - lista_coord[squadra2[1]-1][0])**2 +
                                        (lista_coord[squadra1[1]-1][1] - lista_coord[squadra2[1]-1][1])**2)
                costo = costo + distanza
    return costo/2

def CalcolaCosto2(dict_gironi, lista_coord):
    costo = 0
    for girone in dict_gironi:
        squadre_girone = dict_gironi[girone]
        for squadra1 in squadre_girone:
            for squadra2 in squadre_girone:
                distanza_al_quadrato = ((lista_coord[squadra1[1]-1][0] - lista_coord[squadra2[1]-1][0])**2 +
                                        (lista_coord[squadra1[1]-1][1] - lista_coord[squadra2[1]-1][1])**2)
                costo = costo + distanza_al_quadrato
                # print('Distanza al quadrato da {} a {} è {}'.format(squadra1[1], squadra2[1], distanza_al_quadrato ))
    return costo/2

def CalcolaCosto_minmax(dict_gironi, lista_coord):
    costo_tot = 0
    for girone in dict_gironi:
        squadre_girone = dict_gironi[girone]
        costo_max_girone = 0
        for squadra1 in squadre_girone:
            for squadra2 in squadre_girone:
                distanza_al_quadrato = ((lista_coord[squadra1[1]-1][0] - lista_coord[squadra2[1]-1][0])**2 +
                                        (lista_coord[squadra1[1]-1][1] - lista_coord[squadra2[1]-1][1])**2)
                if distanza_al_quadrato > costo_max_girone:
                    costo_max_girone = distanza_al_quadrato         
        costo_tot = costo_tot + costo_max_girone
    return costo_tot

# -----------------------------------------------
#   MAIN function
# -----------------------------------------------
if __name__ == "__main__":
    # filename = 'Squadre_D1_Maschile.csv'
    filename = 'Squadre_D1_Femminile.csv'
    
    lista_dati = ParseFile(filename)
    # print(lista_dati[0:2])
    # for i, squadra in enumerate(lista_dati):
    #     print('Squadra {}: {}'.format(i+1, lista_dati[i][0]))
    # print('numero di squadre = {}\n'.format(len(lista_dati))) # stampo a video il numero di squadre
    
    lista_coord = [(x,y) for _,x,y in lista_dati] # lista di tuple con solo le coordinate
    # print(lista_coord[0:2])
    
    matrice_costi = CostMatrix(lista_dati)
    # print(matrice_costi)
    
    M = 6
    [dict_gironi, costo_totale] = ILP(M, lista_dati)
    print(dict_gironi)
    
    print("Costo totale somme distanze gironi ottimi = {}".format(CalcolaCosto(dict_gironi, lista_coord)))
    print("Costo totale somme distanze al quadrato gironi ottimi = {}".format(CalcolaCosto2(dict_gironi, lista_coord)))
    print('Somma delle distanze massime al quadrato per i gironi ottimali = {}'.format(CalcolaCosto_minmax(dict_gironi, lista_coord)))
    
    PlotSolution(dict_gironi, lista_coord)
    
    # Stampo file in output
    output_name = 'Gironi_ILP_D1_m.txt'
    output = open(output_name, 'w')
    for key in dict_gironi:
        output.write("{}:\n".format(key))
        girone_k = dict_gironi[key]
        for squadra in girone_k:
            output.write("\t{}\n".format(squadra[0]))
    output.close()
    
    dict_D1_f = {}
    dict_D1_f['Girone1'] = [("NTC LURAGO D'ERBA",18), 
                             ('CIS CHIAVENNA ', 17),
                             ('T.C. COLICO ', 26),
                             ('T.C. LECCO ', 28),
                             ('ALTE GROANE', 3),
                             ('GFG SPORT – David Lloyd', 6)]
    dict_D1_f['Girone2'] = [('SAN GIUSEPPE', 27), 
                             ("OLTREPO' TENNIS ACADEMY", 5),
                             ('SPORTIVA AMP', 2),
                             ('T.C. PAVIA', 9),
                             ('CANOTTIERI TICINO PAVIA ', 13),
                             ('DEMA SPORT 2001', 7)]
    dict_D1_f['Girone3'] = [('JUNIOR TENNIS TRAINING', 23), 
                             ('ECO SPORT', 4),
                             ('T.C. GEMONIO', 20),
                             ("CIRCOLO TENNIS CANTU' ", 25),
                             ('BIEFFESPORT', 21),
                             ('JUNIOR TENNIS MILANO', 8)]
    dict_D1_f['Girone4'] = [('TENNIS CARPENEDOLO', 19), 
                             ('TORBOLE CASAGLIA ', 15),
                             ('OLIMPICA TENNIS REZZATO', 24),
                             ('CANOTTIERI L. BISSOLATI', 16),
                             ('CANOTTIERI FLORA', 12)]
    dict_D1_f['Girone5'] = [('BITETENNIS ', 11), 
                             ('T.C. BERGAMO', 22),
                             ('CANOTTIERI ADDA 1891', 1),
                             ('TENNIS CLUB LODI', 10),
                             ('GTENNIS ', 14)]
    
    print('Costo totale somme distanze gironi reali = {}'.format(CalcolaCosto(dict_D1_f, lista_coord)))
    print('Costo totale somme distanze al quadrato gironi reali = {}'.format(CalcolaCosto2(dict_D1_f, lista_coord)))
    print('Somma delle distanze massime al quadrato per i gironi reali = {}'.format(CalcolaCosto_minmax(dict_D1_f, lista_coord)))
    
    PlotSolution(dict_D1_f, lista_coord)

    
    """
    lista_gironi_D1_f_ufficiale = [[(18,17), (17,28), (28,26), (26,3), (3,6), (6,18)],
                                  [(27,5), (5,2), (2,9), (9,13), (13,7), (7,27)],
                                  [(23,4), (4,20), (20,25), (25,21), (21,8), (8,23)],
                                  [(19,15), (15,24), (24,16), (16,12), (12,19)],
                                  [(11,22), (22,1), (1,10), (10,14), (14,11)]]
    
    costo_tot_ufficiale = CalcolaCosto2(lista_gironi_D1_f_ufficiale, lista_coord)
    print('Costo totale gironi reali = {}'.format(costo_tot_ufficiale))
    
    nostra_sol = [[(28,26), (26,22), (22,18), (18,17), (17,28)],
                  [(27,25), (25,23), (23,20), (20,4), (4, 3), (3,27)],
                  [(21,13), (13,9), (9,5), (5,2), (2,21)],
                  [(12,15), (15,16), (16,19), (19,24), (24,12)],
                  [(1,6), (6,8), (8,10), (10,11), (11,14), (14,1)]]
    
    c = CalcolaCosto2(nostra_sol, lista_coord)
    print('Costo totale gironi reali = {}'.format(c))
    
    # costo_tot = CalcolaCosto2(lista_gironi_D1_f_ottima, lista_coord)
    # print('Costo totale gironi soluzione ottima = {}'.format(costo_tot))
    """