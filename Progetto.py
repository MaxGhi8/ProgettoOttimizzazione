# funzione per leggere il file
def ParseFile(filename):
    doc = open(filename, 'r', encoding = 'utf-8')
    doc.readline() # Salto la prima linea con le intestazioni
    # Leggo i circoli e faccio una lista di liste, ogni sottolista contiene
    # ordinatamente nome squadra come stringa e le due coordinate come float
    D = []
    for row in doc:
        row = row.split(',')
        # print(row)
        D.append( [row[1], float(row[2]), float(row[3])] )
    return D

from math import sqrt

def CostList(Ls):
    # INPUT: Ls è la lista in output da ParseFile
    # OUTPUT: Cs lista dei costi fatta come [ ('partenza', 'arrivo', costo), (), ...]
    n = len(Ls)
    Cs = []
    for i in range(n):
        for j in range(n):
            if j > i:
                costo = sqrt((Ls[i][1] - Ls[j][1])**2 + (Ls[i][2] - Ls[j][2])**2)
                Cs.append( (Ls[i][0], Ls[j][0], costo) )
    return Cs


# Import the NetworkX library per fare i grafi
import networkx as nx

def BuildGraph(Cs):
    # INPUT: Cs lista dei costi, output di CostList
    # OUTPUT: G grafo indiretto con nodi=squadre, lati=costi 
    
    G = nx.Graph() # costruisce grafo indiretto    
    # aggiungo i nodi 
    G.add_nodes_from(list(map(lambda x: x[0], Cs)))
    # aggiungo i lati con i costi associati
    G.add_weighted_edges_from(Cs) #se esiste già non lo riaggiunge
    return G

import pylab as pl
from matplotlib import collections as mc

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
    
from pyomo.environ import ConcreteModel, Var, Objective, Constraint, SolverFactory
from pyomo.environ import maximize, Binary, RangeSet, PositiveReals, ConstraintList

def TSP(C):
    # Dimension of the problem
    n, n = C.shape
    
    print('city:', n)
    
    model = ConcreteModel()

    # The set of indces for our variables
    model.I = RangeSet(1, n)   # NOTE: it is different from "range(n)"
    # the RangeSet is in [1,..n], the second is [0,n(
    
    model.J = RangeSet(1, n)
    
    # Variable definition
    model.X = Var(model.I, model.J, within=Binary)
    
    # Variables for the MTZ subtour constraints
    model.U = Var(model.I, within=PositiveReals)
    
    # Objective function
    model.obj = Objective(
        expr = sum(C[i-1,j-1] * model.X[i,j] for i,j in model.X),
        sense = maximize
    )

    
    # The constraints
#    model.Indegree = Constraint(model.I, 
#                                rule =  lambda m, i: sum(m.X[i,j] for j in m.J) == 1)

    def OutDegRule(m, i):
        return sum(m.X[i,j] for j in m.J) == 1
    model.Outdegree = Constraint(model.I, rule = OutDegRule)
    
    def InDegRule(m, j):
        return sum(m.X[i,j] for i in m.I) == 1
    model.Indegree = Constraint(model.J, rule = InDegRule)
    
    # easily for forbding "pairs" tour
    model.pairs = ConstraintList()
    for i in model.I:
        for j in model.J:
            model.pairs.add(expr=model.X[i,j] + model.X[j,i] <= 1)
    
    # MTZ constraints for forbidding subtours
    model.subtour = ConstraintList()
    for i in model.I:
        for j in model.J:
            if i > 1 and j > 1 and i != j:
                model.subtour.add(model.U[i] - model.U[j] + (n-1)*model.X[i,j] <= (n-1) - 1)
    
    # Solve the model
    sol = SolverFactory('gurobi').solve(model, tee=True)
    
    sol_json = sol.json_repn()
    #print(sol_json)
    if sol_json['Solver'][0]['Status'] != 'ok':
        print("qualcosa è andato storto")
        return None

    # Retrieve the solution: as a list of edges in the optimal tour
    return [(i-1, j-1) for i, j in model.X if model.X[i,j]() > 0.0]    
    
    

# -----------------------------------------------
#   MAIN function
# -----------------------------------------------
if __name__ == "__main__":
    filename = 'Squadre_D1_Maschile-giusto.csv'
    
    lista_dati = ParseFile(filename)
    # print(lista_dati[0:2])
    print('numero di squadre = {}\n'.format(len(lista_dati))) #stampo a video il numero di squadre
    
    # lista_coord = [(x,y) for _,x,y in lista_dati] #lista di tuple con solo le coordinate
    # print(lista_coord[0:2])
    
    lista_costi = CostList(lista_dati)
    # print(lista_costi)
    print('numero di coppie diverse = {}\n'.format(len(lista_costi))) # per check
    
    G = BuildGraph(lista_costi)
    print( 'numero di nodi del grafo = {}\n'.format(nx.number_of_nodes(G)) )
    print( 'numero di lati del grafo = {}\n'.format(nx.number_of_edges(G)) )
    
    #PlotTour(lista_coord, RandomTour, values)