from math import ceil
import numpy as np
# Libreria pyomo per il modello di ottimizzazione
from pyomo.environ import ConcreteModel, Var, Objective, SolverFactory
from pyomo.environ import minimize, Binary, RangeSet, ConstraintList
# Per plottare i risultati
import matplotlib.pyplot as plt

# funzione per leggere il file
def ParseFile(filename):
    """
    Questa funzione prende in input il file dei dati creato e restituisce
    Ls = lista con i nomi delle squadre e le due coordinate 
    """
    doc = open(filename, 'r', encoding = 'utf-8')
    doc.readline() # Salto la prima linea con le intestazioni
    # Leggo i circoli e faccio una lista di liste, ogni sottolista contiene
    # ordinatamente nome squadra come stringa e le due coordinate come float
    Ls = []
    for row in doc:
        row = row.split(',') # poiché il file è un csv
        Ls.append( [row[1], float(row[2]), float(row[3])] )
    return Ls 
 
def Tau(m,M):
    """
    Parameters:
    ----------
    n : numero intero
        Numero totale di squadre nel torneo
    M : numero intero
        Numero di squadre massimo in ogni girone

    Returns
    -------
    tau : lista
          La lista alla posizione i-esima contiene il numero di squadre del 
          girone i-esimo
    """
    tau = list(M for _ in range(ceil(m/M))) # ceil = approssimazione superiore
    resto = sum(tau) - m
    for i in range(resto):
        tau[i] = tau[i] - 1
    return tau

def cKM(data, centre, max_it, M):
    """
    Constrained K-Means's Algorithm for the k-means clustering

    Parameters
    ----------
    data : list
        output di ParseFile
    centre : list
        lista di [(latitudine, longitudine), ...] = initial centres
    max_it : int
        numero massimo di iterazioni
    M : int
        numero massimo di squadre per girone

    Returns
    -------
    cluster : list
        lista contenente liste, ogni sottolista contiene un girone
    centre : list
        lista dei centri
    it : int
        numero di iterazioni impiegate all'algoritmo per convergere
    """

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
    """
    Risolutore lineare del sottoproblema per il constrained k-means cluster
    
    Parameters
    ----------
    centre : list
        lista di [(latitudine, longitudine), ...] = initial centres
    data : list
        output di ParseFile
    tau : list
        output di Tau(n,M)

    Returns
    -------
    T : list
        lista contenente k sottoliste, ogni sottolista contiene dei valori in 
        {0,1}: 1 se la variabile i è nel cluster k e 0 altrimenti

    """
    
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
        
    solver = SolverFactory('gurobi') 
    
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

def PlotSolution(cluster, centre):
    # k = len(cluster)
    for girone in cluster:
        plt.scatter(list(map(lambda x: x[1], girone)), list(map(lambda x: x[2], girone)))
    for centro in centre:
        plt.scatter(centro[0],centro[1], c='black', marker='*')
    plt.show()
    
def PlotCentre(centre):
    for centro in centre:
        plt.scatter(centro[0],centro[1])
    plt.show()
    
def CalcolaCosto(cluster):
    """ Somma delle distanze al quadrato """
    costo = 0
    for girone in cluster:
        for i, squadra1 in enumerate(girone):
            for j, squadra2 in enumerate(girone):
                if i > j:
                    costo = costo + (squadra1[1] - squadra2[1])**2 + (squadra1[2] - squadra2[2])**2
    return costo

def MaxDistanza(coord, k):
    """ 
    INPUT: - coord: lista di coordinate
           - k: numero totale di cluster
    OUTPUT: - centre inizializzati massimizzando la distanza massima
    """
    centre = []
    m = len(coord) # numero dei punti
    index_centre = np.random.choice(list(range(m)), 1)
    centre.append(coord[int(index_centre)]) # scegliamo il primo punto casuale
    coord.remove(centre[0])
    for i in range(k-1):
        dist_max = 0
        centre_new = coord[0]
        for squadra in coord:
            dist = 0
            for centro in centre:
                dist = dist + (squadra[0] - centro[0])**2 + (squadra[1] - centro[1])**2
            if dist > dist_max:
                dist_max = dist
                centre_new = squadra
        centre.append(centre_new)
        coord.remove(centre_new)
    return centre

# -----------------------------------------------
#   MAIN function
# -----------------------------------------------
if __name__ == "__main__": 
    genere = 'femminile'
    # genere = 'maschile'
    if genere == 'maschile':
        filename = 'Squadre_D1_Maschile.csv' # D1 maschile
        output_name = 'Gironi_ILP_D1_m_kmeans.txt'
    if genere == 'femminile':
        filename = 'Squadre_D1_Femminile.csv' # D1 femminile
        output_name = 'Gironi_ILP_D1_f_kmeans.txt'
    
    lista_dati = ParseFile(filename)
    m = len(lista_dati) # numero di squadre totali
    # Per debugging
    # print('Dati: {}\n'.format(lista_dati[0:2]))
    # print('numero di squadre = {}\n'.format(len(lista_dati))) # stampo a video il numero di squadre
    
    lista_coord = [(x,y) for _,x,y in lista_dati] # lista di tuple con solo le coordinate
    # Per debugging
    # print('Coordinate:{}\n'.format(lista_coord))
    
    M = 6
    print('Tau = {}\n'.format(Tau(m,M)))
    k = len(Tau(m,M)) # numero di gironi 
    
    ### Scelta di centri casuali
    # index_centre = np.random.choice(list(range(m)), k, False)
    # centre = []
    # for i in index_centre:
    #     centre.append(lista_coord[i])
    # print(centre)
    
    # OSSERVAZIONE: il metodo farthest_distance comporta un costo più stabile 
    # rispetto a quello che segue dalla scelta dei centri casuali.
    
    ### Scelta di centri secondo la massima distanza
    centre = MaxDistanza(lista_coord, k)
    print('Centri: {}\n'.format(centre))
    PlotCentre(centre)
    
    # Assegnazione dei gironi
    max_it = 50
    cluster, centre, it = cKM(lista_dati, centre, max_it, M)
    print('Gironi: {}\n'.format(cluster))
    # cluster è lista di liste di liste
    
    PlotSolution(cluster, centre)
    
    print('Costo: {}\n'.format(CalcolaCosto(cluster)))

    output = open(output_name, 'w')
    for i, girone in enumerate(cluster):
        output.write("Girone {}:\n".format(i+1))
        for squadra in girone:
            output.write("\t {}\n".format( squadra[0] ))
    output.close()
