
from math import sqrt, ceil
import numpy as np
# Import the NetworkX library per fare i grafi
import networkx as nx

from pyomo.environ import ConcreteModel, Var, Objective, SolverFactory
from pyomo.environ import Binary, RangeSet, ConstraintList

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


def CostList(Ls):
    # INPUT: Ls è la lista in output da ParseFile
    # OUTPUT: Cs lista dei costi fatta come [ ('partenza', 'arrivo', costo), (), ...]
    n = len(Ls)
    Cs = []
    for i in range(n):
        for j in range(n):
            if j != i:
                costo = sqrt((Ls[i][1] - Ls[j][1])**2 + (Ls[i][2] - Ls[j][2])**2)
                # Delta_phi = (Ls[i][1] - Ls[j][1])/2
                # costo = 2*6360*np.arcsin(sqrt( np.sin(Delta_phi)**2 +
                #                          ( 1-np.sin(Delta_phi)**2 - np.sin((Ls[i][1] + Ls[j][1])/2 )**2 )*
                #                          np.sin((Ls[i][2] - Ls[j][2])/2)**2))
            else:
                costo = 10000 # evitiamo che le squadre si scontrino con sè stesse
            Cs.append( (i+1, j+1, costo) )
            ### ACTUNG: nome i-esimo associato alla posizione (i+1)-esima in Ls ###
    return Cs


def BuildGraph(Cs):
    # INPUT: Cs lista dei costi, output di CostList
    # OUTPUT: G grafo indiretto con nodi=squadre, lati=costi 
    
    G = nx.DiGraph() # costruisce grafo diretto    
    # aggiungo i nodi 
    G.add_nodes_from(list(map(lambda x: x[0], Cs)))
    # aggiungo i lati con i costi associati
    G.add_weighted_edges_from(Cs) # se esiste già non lo riaggiunge
    return G
   
    
def Tau(n, M):
    # Costruzione di tau
    tau = list(M for _ in range(ceil(n/M)))
    resto = sum(tau) - n
    for i in range(resto):
        tau[i] = tau[i] - 1
    return tau
    
    
def VRPCut(M, Ls, itermax = 50):
    # INPUT: - M cardinalità gironi (?ACTUNG: M divide o no?)
    #        - Ls lista delle squadre
    
    # Build a directed graph out of the data
    n = len(Ls) #numero di dati (numero di squadre)
    Cs = CostList(Ls)
    G = BuildGraph(Cs)
    tau = Tau(n, M)
    kappa = len(tau) # numero di cluster (= gironi)

    # Build ILP Model
    model = ConcreteModel()
    
    model.N = RangeSet(n)
    model.K = RangeSet(kappa)

    # TODO: introduce only usefull variables (no arc, no variable)
    model.x = Var(model.N, model.N, model.K, domain=Binary)

    # Objective function of arc variables
    model.obj = Objective(expr=sum(G[i][j]['weight']*model.x[i, j, k] for i, j in G.edges() for k in model.K))

    # Vincoli archi uscenti
    model.outdegree = ConstraintList()
    for i in model.N:
        model.outdegree.add(expr=sum(model.x[i, j, k] for j in model.N for k in model.K) == 1)

    # Vincoli archi entranti    
    model.indegree = ConstraintList()
    for j in model.N:
        model.indegree.add(expr=sum(model.x[i, j, k] for i in model.N for k in model.K) == 1)
        
    # Vincoli cicli
    model.cicli = ConstraintList()
    for i in model.N:
        for j in model.N:
            for k in model.K:
                model.cicli.add(expr = sum(model.x[j, s, k] for s in model.N) - model.x[i, j, k] >= 0)
        
    # Vincoli unicità
    model.unic = ConstraintList()
    for i in model.N:
        for j in model.N:
            model.unic.add(expr = sum(model.x[i, j, k] for k in model.K) <= 1)
            
    # Vincolo tau
    model.tau = ConstraintList()
    for k in model.K:
        model.tau.add(expr = sum(model.x[i, j, k] for i in model.N for j in model.N) == tau[k-1])
    
    # Vincoli no coppie
    model.onedir = ConstraintList()
    if M > 3:
        for k in model.K:
            for i, j in G.edges():
                if i < j:
                    model.onedir.add(model.x[i, j, k] + model.x[j, i, k] <= 1)
                
    # Vincoli subtours
    model.subtours = ConstraintList()            
    
    solver = SolverFactory('gurobi')    

    it = 0
    while it <= itermax:
        it += 1

        sol = solver.solve(model, tee=False)

        # Get a JSON representation of the solution
        sol_json = sol.json_repn()
        # Check solution status
        if sol_json['Solver'][0]['Status'] != 'ok':
            return None
        
        counter = 0
        selected = []
        for k in model.K:
            print('iterazione: {}\t k = {}\n'.format(it, k))
            
            selected_k = []
            for i in model.N:
                for j in model.N:
                    if model.x[i, j, k]() > 0.5:
                        selected_k.append((i, j))
            
            # lista_coord = [(x,y) for _,x,y in Ls]
            # PlotSolution(lista_coord, selected)

            Cuts_k = SubtourElimin(selected_k, M, M-1)
            print("Cuts_k: {}\n".format(Cuts_k))

            if Cuts_k == []:
                counter = counter + 1
                
            for faccia in model.K:
                for S in Cuts_k:
                    model.subtours.add(sum(model.x[i, faccia] for i in S) <= len(S)-1)
            
            selected.append(selected_k)
            
        if counter == k:
            print("Optimal solution found !!")
            break
        
    return selected

def SubtourElimin(selected, M, N=-1):
    G = nx.Graph()

    for i, j in selected:
        G.add_edge(i, j)

    # Cys = nx.simple_cycles(G)
    Cys = nx.cycle_basis(G)
    # print('Selected={}\n'.format(selected))
    # print('Cys={}\n'.format(Cys))

    SubCycles = []
    for cycle in Cys:
        lun = len(cycle)
        if lun != M and lun != N:
            SubCycles.append(cycle)
    
    SUBTOURS = []
    for Ciclo in SubCycles:
        Subtours = []
        Subtours_reversed = []
        for nodo in Ciclo:
            for lato in selected:
                if lato[0] == nodo:
                    Subtours.append(lato)
                    Subtours_reversed.append( (lato[1], lato[0]))
        SUBTOURS.append(Subtours)
        SUBTOURS.append(Subtours_reversed)

    return SUBTOURS

def DisegnaPunto(A, ax):
    """
    Disegna un punto nel piano
    """
    ax.plot([A[0]], [A[1]], 'bo', alpha=0.5)

def DisegnaSegmento(A, B, ax):
    """ 
    Disegna un segmento nel piano dal punto A a al punto B
    Vedi manuale a: http://matplotlib.org/api/pyplot_api.html
    """
    # Disegna il segmento
    ax.plot([A[0], B[0]], [A[1], B[1]], 'b', lw=0.75)
    # Disegna gli estremi del segmento
    DisegnaPunto(A, ax)
    DisegnaPunto(B, ax)

def PlotSolution(Xs, Es):
    fig, ax = plt.subplots()
    for i, j in Es:
        DisegnaSegmento(Xs[i-1], Xs[j-1], ax)

    ax.scatter([i for i, j in Xs[1:]], [j for i, j in Xs[1:]], alpha=0.3, cmap='viridis')

    for i in range(len(Xs[:])):
        ax.annotate(str(i+1), Xs[i])

    plt.axis('square')
    plt.axis('off')

# -----------------------------------------------
#   MAIN function
# -----------------------------------------------
if __name__ == "__main__":
    # filename = 'Squadre_D1_Maschile-giusto.csv'
    filename = 'Squadre_D1_Femminile.csv'
    
    lista_dati = ParseFile(filename)
    # print(lista_dati[0:2])
    print('numero di squadre = {}\n'.format(len(lista_dati))) # stampo a video il numero di squadre
    
    lista_coord = [(x,y) for _,x,y in lista_dati] # lista di tuple con solo le coordinate
    # print(lista_coord[0:2])
    
    lista_costi = CostList(lista_dati)
    # print(lista_costi)
    # print('numero di coppie = {}\n'.format(len(lista_costi))) # per check
    
    G = BuildGraph(lista_costi)
    # nx.draw(G) 
    # print( 'numero di nodi del grafo = {}\n'.format(nx.number_of_nodes(G)) )
    # print( 'numero di lati del grafo = {}\n'.format(nx.number_of_edges(G)) )
    # print(G[1][1]['weight'])
    
    M = 6
    lista_gironi = VRPCut(M, lista_dati, 10000)
    # print(lista_gironi)
    # print(SubtourElimin(lista_gironi, M+1))
    PlotSolution(lista_coord, lista_gironi)

    
    output = open('Gironi.txt', 'w')
    for i, girone in enumerate(SubtourElimin(lista_gironi, M+1)):
        output.write("Girone {}:\n".format(i+1))
        for lato in girone:
            output.write("\t {}\n".format( lista_dati[lato[0]-1][0] ))
    output.close()
    
    