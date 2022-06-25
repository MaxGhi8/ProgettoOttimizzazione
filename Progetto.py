
from math import sqrt

# Import the NetworkX library per fare i grafi
import networkx as nx

import pyomo
from pyomo.environ import ConcreteModel, Var, Objective, Constraint, SolverFactory
from pyomo.environ import maximize, Binary, RangeSet, PositiveReals, ConstraintList

import pylab as pl
from matplotlib import collections as mc

# funzione per leggere il file
def ParseFile(filename):
    # OUTPUT: Ls lista con i nomi delle squadre e le due coordinate 
    doc = open(filename, 'r', encoding = 'utf-8')
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


def PlotTour(Ps, Ls, values):
    lines = [[Ps[i], Ps[j]] for i,j in Ls]
    fig, ax = pl.subplots()

    lc = mc.LineCollection(lines, linewidths=[1.5 if x > 0.501 else 1 for x in values],
                           colors=['blue' if x > 0.501 else 'orange' for x in values])
    
    ax.add_collection(lc)
    ax.scatter([i for i,j in Ps], [j for i,j in Ps], 
                s=20, alpha=0.8, color='red')
    
    ax.autoscale()
    ax.margins(0.1)
    ax.axis('equal')
    pl.show()
    
def VRPCut(M, Ls):
    # INPUT: - M cardinalità gironi (?ACTUNG: M divide o no?)
    #        - Ls lista delle squadre
    
    # Build a directed graph out of the data
    Cs = CostList(Ls)
    G = BuildGraph(Cs)

    # Build ILP Model
    model = ConcreteModel()
    
    n = len(Ls) # numero di nodi (squadre)
    model.N = RangeSet(n)

    # TODO: introduce only usefull variables (no arc, no variable)
    model.x = Var(model.N, model.N, domain=Binary)

    # Objective function of arc variables
    model.obj = Objective(expr=sum(G[i][j]['weight']*model.x[i, j] for i, j in G.edges()))

    # Vincoli archi uscenti
    model.outdegree = ConstraintList()
    for i in model.N:
        model.outdegree.add(expr=sum(model.x[v, w] for v, w in G.out_edges(i)) == 1)

    # Vincoli archi entranti
    model.indegree = ConstraintList()
    for j in model.N:
        model.indegree.add(expr=sum(model.x[v, w] for v, w in G.in_edges(j)) == 1)
    
    # Solve the model
    sol = SolverFactory('glpk').solve(model, tee=True)

    # # Get a JSON representation of the solution
    # sol_json = sol.json_repn()
    # # Check solution status
    # if sol_json['Solver'][0]['Status'] != 'ok':
    #     return None

    selected = []
    for i, j in model.x:
        if model.x[i, j]() > 0:
            selected.append((i, j))

    return selected

# -----------------------------------------------
#   MAIN function
# -----------------------------------------------
if __name__ == "__main__":
    filename = 'Squadre_D1_Maschile-giusto.csv'
    
    lista_dati = ParseFile(filename)
    # print(lista_dati[0:2])
    print('numero di squadre = {}\n'.format(len(lista_dati))) #stampo a video il numero di squadre
    
    lista_coord = [(x,y) for _,x,y in lista_dati] #lista di tuple con solo le coordinate
    # print(lista_coord[0:2])
    
    lista_costi = CostList(lista_dati)
    # print(lista_costi)
    print('numero di coppie = {}\n'.format(len(lista_costi))) # per check
    
    G = BuildGraph(lista_costi)
    print( 'numero di nodi del grafo = {}\n'.format(nx.number_of_nodes(G)) )
    print( 'numero di lati del grafo = {}\n'.format(nx.number_of_edges(G)) )
    print(G[0][1]['weight'])
    
    Es = VRPCut(2, lista_dati)
    
    n = len(lista_dati)
    values = [1 for _ in range(n)]
    PlotTour(lista_coord, Es, values)