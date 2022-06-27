# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 10:05:28 2022

@author: alice e max
"""

from math import sqrt, ceil
import numpy as np

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
 
def Tau(m,M):
    # Costruzione di tau
    tau = list(M for _ in range(ceil(m/M)))
    resto = sum(tau) - m
    for i in range(resto):
        tau[i] = tau[i] - 1
    return tau

def cKM(data, centre, max_it, M):
# Constrained K-Means's Algorithm for the k-means clustering
# INPUT: data points, initial centres, max number of iterations
# OUTPUT: cluster indices for data points, clusters' centres, n. iterations
    
    # data: output di ParseFile
    # centre: lista di [(latitudine, longitudine), ...]
    m = len(data) # number of data points
    k = len(centre) # number of clusters
    tau = Tau(m, M)
    
    it = 0 # number of iteration performed so far
    convergence = False # 'true' when no improvement has been obtained
    while (not convergence) and (it < max_it):
        it = it + 1
        print('iterata: {}\n'.format(it))
        # PASSO 1        
        T = ClusterAssignment(centre, data, tau)
        
        # T = [ [], T_h, ..., [] ], dove T_h = [[], T_i, ..., []]
        # PASSO 2
        centre_new = [] # lista con i centri
        for h in range(k):
            if sum(T[h]) > 0:
                centre_new.append( ( (sum( list(map(lambda x,y: x[1]*y, data, T[h])) ))/(sum(T[h])) ,
                                  (sum( list(map(lambda x,y: x[2]*y, data, T[h])) ))/(sum(T[h])) ))
            else: 
                centre_new.append( centre[h] )
        
        if centre_new == centre: # convergence has been reached:
            convergence = True
        else:
            centre = centre_new
            
    # Assemblaggio cluster
    cluster = []
    for h in range(k):
        cluster_h = []
        for i in range(m):
            if T[h][i] > 0.5:
                cluster_h.append(data[i])
        cluster.append(cluster_h)
        
    return cluster, centre, it

def Distanza(dato, centro):
    distanza = (dato[1] - centro[0])**2 + (dato[2] - centro[1])**2
    return distanza

def ClusterAssignment(centre, data, tau):
    
    m = len(data) # number of data points
    k = len(centre) # number of clusters

    # Build ILP Model
    model = ConcreteModel()
    
    # Indici
    model.M = RangeSet(m)
    model.K = RangeSet(k)

    # Variabili
    model.x = Var(model.K, model.M, domain=Binary)

    # Objective function of arc variables
    model.obj = Objective( expr=sum(model.x[h,i]*0.5*Distanza(data[i-1], centre[h-1]) for i in model.M for h in model.K), 
                          sense = minimize)
    
    # Vincoli archi uscenti
    model.tau = ConstraintList()
    for h in model.K:
        model.tau.add(expr=sum(model.x[h,i] for i in model.M) >= tau[h-1])
        
    model.elemento = ConstraintList()
    for i in model.M:
        model.elemento.add(expr=sum(model.x[h,i] for h in model.K) == 1)
        
    solver = SolverFactory('glpk') 
    
    sol = solver.solve(model, tee=False)

    # Get a JSON representation of the solution
    sol_json = sol.json_repn()
    # Check solution status
    if sol_json['Solver'][0]['Status'] != 'ok':
        return None
    
    T = []
    for h in model.K:
        T_h = []
        for i in model.M:
            if model.x[h,i]() > 0.5:
                T_h.append(1)
            else:
                T_h.append(0)
        T.append(T_h)
        
    return T

def PlotSolution(cluster):
    k = len(cluster)
    for girone in cluster:
        plt.scatter(list(map(lambda x: x[1], girone)), list(map(lambda x: x[2], girone)))
    plt.show()
    

# -----------------------------------------------
#   MAIN function
# -----------------------------------------------
if __name__ == "__main__": 
    filename = 'Squadre_D1_Maschile.csv'
    # filename = 'Squadre_D1_Femminile.csv'
    
    lista_dati = ParseFile(filename)
    m = len(lista_dati)
    print('Dati: {}\n'.format(lista_dati[0:2]))
    # print('numero di squadre = {}\n'.format(len(lista_dati))) # stampo a video il numero di squadre
    
    lista_coord = [(x,y) for _,x,y in lista_dati] # lista di tuple con solo le coordinate
    print('Coordinate:{}\n'.format(lista_coord))
    
    M = 6
    print(Tau(m,M))
    index_centre = np.random.choice(list(range(m)), len(Tau(m,M)), False)
    centre = []
    for i in index_centre:
        centre.append(lista_coord[i])
    print(centre)
    max_it = 50
    cluster, centre, it = cKM(lista_dati, centre, max_it, M)
    print('Gironi: {}\n'.format(cluster))
    # cluster è lista di liste di liste
    
    PlotSolution(cluster)

    output = open('Gironi_D1_Maschile.txt', 'w')
    # output = open('Gironi_D1_Femminile.txt', 'w')
    for i, girone in enumerate(cluster):
        output.write("Girone {}:\n".format(i+1))
        for squadra in girone:
            output.write("\t {}\n".format( squadra[0] ))
    output.close()