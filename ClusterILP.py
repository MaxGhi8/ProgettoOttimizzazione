from math import sqrt, ceil
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
    Ls = lista cdi tuple di elementi (nome squadra, le due coordinate) 
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

def CostMatrix(Ls):
    """
    INPUT: Ls è la lista in output da ParseFile
    OUTPUT: Cm matrice con i costi Cm(i,j) = costo da squadra i a squadra j
    Modificare questa funzione se si vuole considerare dei pesi diversi, ad 
    esempio le reali distanze tra comuni
    """
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
    tau = list(M for _ in range(ceil(n/M))) # ceil = approssimazione superiore
    resto = sum(tau) - n
    for i in range(resto):
        tau[i] = tau[i] - 1
    return tau
    
    
def ILP(M, Ls):
    """
    Risolve il problema lineare intero per trovare la soluzione ottima che 
    minimizza le distanze delle squadre all'interno dei gironi (o le distanze
    al quadrato).

    Parameters
    ----------
    M : intero
        Numero di squadre massimo in ogni girone
    Ls : lista
        lista di tuple (output di parsefile)

    Returns
    -------
    dict_gironi : dizionario
        dizionario con chiave 'Gironek' e come valore associato una lista con 
        i nomi el'indice delle squadre associate a quel girone
    costo_totale : float
        (valore della funzione obbiettivo)/2 = costo totale del cluastering trovato

    """    
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
    """
    Fa il plot del cluster trovato

    Parameters
    ----------
    dict_gironi : dizionario
        output della funzione ILP
    lista_coord : lista
        lista di tuple, con al posto i-esimo la tupla con le coordinate
        associate alla squadra i-esima

    Returns
    -------
    None.

    """
    for girone in dict_gironi:
        squadre_girone = dict_gironi[girone]
        plt.scatter( list(map(lambda x: lista_coord[x[1]-1][0], squadre_girone)),
                    list(map(lambda x: lista_coord[x[1]-1][1], squadre_girone)) )
    plt.show()
    
def CalcolaCosto(dict_gironi, lista_coord):
    """
    Parameters
    ----------
    dict_gironi : dizionario
        output della funzione ILP
    lista_coord : lista
        lista di tuple, con al posto i-esimo la tupla con le coordinate
        associate alla squadra i-esima

    Returns
    -------
    type = float
      La somma delle distanze tra le squadre all'interno dei cluster e poi 
    somma sui cluster (una sorta di costo del cluster)
    """
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
    """
    Parameters
    ----------
    dict_gironi : dizionario
        output della funzione ILP
    lista_coord : lista
        lista di tuple, con al posto i-esimo la tupla con le coordinate
        associate alla squadra i-esima

    Returns
    -------
    type = float
      La somma delle distanze al quadrato tra le squadre all'interno dei 
    cluster e poi somma sui cluster (una sorta di costo del cluster)
    """
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
    """
    Parameters
    ----------
    dict_gironi : dizionario
        output della funzione ILP
    lista_coord : lista
        lista di tuple, con al posto i-esimo la tupla con le coordinate
        associate alla squadra i-esima

    Returns
    -------
    type = float
      Calcola le massime distanze al quadrato tra le squadre all'interno dei 
    cluster e poi somma sui cluster (una sorta di costo del cluster)
    """
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
    # genere = 'femminile'
    genere = 'maschile'
    if genere == 'maschile':
        filename = 'Squadre_D1_Maschile.csv' # D1 maschile
        output_name = 'Gironi_ILP_D1_m.txt'
        # SOLUZIONE PROPOSTA DALLA FIT
        dict_D1 = {}
        dict_D1['Girone1'] = [('TENNIS E SPORTS OPEN', 17),
                              ('C.S. RANICA SEZIONE TENNIS', 62),
                              ('POLISPORTIVA CURNO', 35),
                              ('TENNIS CLUB BERGAMO', 2),
                              ('TENNIS CLUB STEZZANO', 9),
                              ('U.S. SCANZOROSCIATE', 18)]
        dict_D1['Girone2'] = [('TENNIS TIME CHAMPION', 28),
                              ('FAGNANO TENNIS', 31),
                              ('TENNIS CLUB MILANO – A', 42),
                              ('TENNIS CLUB LUINO', 6),
                              ('T.C. GALLARATE - A ', 30),
                              ('ECO SPORT', 7)]
        dict_D1['Girone3'] = [('CANOTTIERI BALDESIO', 24),
                              ('ASDR TENNIS SUZZARA', 29),
                              ('SANZENOTENNIS 2011', 15),
                              ('T.C. CREMA', 49),
                              ('T.T. 2021 ERIDANEA – VIADANA', 56),
                              ('ELNIK SSD', 44)]
        dict_D1['Girone4'] = [('CIRCOLO TENNIS PARABIAGO', 10),
                              ('ECOSPORT – PONTENUOVO', 43),
                              ('GGIOIASPORT', 3),
                              ('TENNIS CLUB MILANO -B', 1),
                              ('MOLINETTO COUNTRY CLUB', 37),
                              ('T.C. GALLARATE - B', 41)]
        dict_D1['Girone5'] = [('TENNIS COMO', 23),
                              ('TENNIS CLUB ROVELLASCA', 13),
                              ('JUNIOR TENNIS TRAINING', 45),
                              ('TENNIS CLUB MORBEGNO', 58),
                              ('CIS CHIAVENNA ', 14),
                              ('ALTE GROANE', 11)]
        dict_D1['Girone6'] = [('TENNIS CLUB SESTO', 47),
                              ('ASD TENNISTIME BG', 40),
                              ('FORTENNIS', 26),
                              ('TENNIS CLUB BAGNATICA', 50),
                              ('TENNIS CLUB MILANO – C', 51),
                              ('ALBANI TENNIS PROJECT', 60)]
        dict_D1['Girone7'] = [('LA DOMINANTE', 34),
                              ('CLUB TENNIS CERIANO', 52),
                              ('CT CARATE', 57),
                              ('TENNIS CLUB ARCORE', 61),
                              ('TENNIS SENAGO', 25),
                              ('TENNIS CLUB MILANO – D', 36)]
        dict_D1['Girone8'] = [('TENNIS CLUB OGGIONO', 19),
                              ('TENNIS CLUB LECCO', 39),
                              ('CIRCOLO TENNIS ROSEDA', 21),
                              ('T.C. 3 STELLE NAUTICA BOLIS', 38),
                              ('ASD TC TIRANO', 59)]
        dict_D1['Girone9'] = [('GFG SPORT – David Lloyd', 32),
                              ('PRO TENNIS SCHOOL', 16),
                              ('TORBOLE CASAGLIA', 54),
                              ('TIMING', 27),
                              ('TENNIS FUTURA (S.FILIPPO)', 22)]
        dict_D1['Girone10'] = [('OFFICINE TENNIS BRESSO', 48),
                              ('T.C. 30 PARI', 8),
                              ('CANOTTIERI MILANO ', 5),
                              ('JUNIOR TENNIS MILANO', 20),
                              ('TENNIS CLUB TREZZO', 33)]
        dict_D1['Girone11'] = [('SELVA ALTA ', 53),
                              ('CANOTTIERI TICINO PAVIA', 4),
                              ('HARBOUR CLUB MILANO', 46),
                              ('MILAGO TENNIS CENTER', 55),
                              ('SPORTREND TENNIS', 12)]
    if genere == 'femminile':
        filename = 'Squadre_D1_Femminile.csv' # D1 femminile
        output_name = 'Gironi_ILP_D1_f.txt'
        # SOLUZIONE PROPOSTA DALLA FIT
        dict_D1 = {}
        dict_D1['Girone1'] = [("NTC LURAGO D'ERBA",18), 
                                 ('CIS CHIAVENNA ', 17),
                                 ('T.C. COLICO ', 26),
                                 ('T.C. LECCO ', 28),
                                 ('ALTE GROANE', 3),
                                 ('GFG SPORT – David Lloyd', 6)]
        dict_D1['Girone2'] = [('SAN GIUSEPPE', 27), 
                                 ("OLTREPO' TENNIS ACADEMY", 5),
                                 ('SPORTIVA AMP', 2),
                                 ('T.C. PAVIA', 9),
                                 ('CANOTTIERI TICINO PAVIA ', 13),
                                 ('DEMA SPORT 2001', 7)]
        dict_D1['Girone3'] = [('JUNIOR TENNIS TRAINING', 23), 
                                 ('ECO SPORT', 4),
                                 ('T.C. GEMONIO', 20),
                                 ("CIRCOLO TENNIS CANTU' ", 25),
                                 ('BIEFFESPORT', 21),
                                 ('JUNIOR TENNIS MILANO', 8)]
        dict_D1['Girone4'] = [('TENNIS CARPENEDOLO', 19), 
                                 ('TORBOLE CASAGLIA ', 15),
                                 ('OLIMPICA TENNIS REZZATO', 24),
                                 ('CANOTTIERI L. BISSOLATI', 16),
                                 ('CANOTTIERI FLORA', 12)]
        dict_D1['Girone5'] = [('BITETENNIS ', 11), 
                                 ('T.C. BERGAMO', 22),
                                 ('CANOTTIERI ADDA 1891', 1),
                                 ('TENNIS CLUB LODI', 10),
                                 ('GTENNIS ', 14)]
    
    lista_dati = ParseFile(filename)
    # Per debugging
    # print(lista_dati[0:2])
    # for i, squadra in enumerate(lista_dati):
    #     print('Squadra {}: {}'.format(i+1, lista_dati[i][0]))
    # print('numero di squadre = {}\n'.format(len(lista_dati))) # stampo a video il numero di squadre
    
    lista_coord = [(x,y) for _,x,y in lista_dati] # lista di tuple con solo le coordinate
    # Per debugging
    # print(lista_coord[0:2])
    
    matrice_costi = CostMatrix(lista_dati)
    # print(matrice_costi)
    
    M = 6
    [dict_gironi, costo_totale] = ILP(M, lista_dati) # cluster ottimale
    
    # stampo a video i gironi
    print(dict_gironi)
    # stampo a video i costi
    print("Costo totale somme distanze gironi ottimi = {}".format(CalcolaCosto(dict_gironi, lista_coord)))
    print("Costo totale somme distanze al quadrato gironi ottimi = {}".format(CalcolaCosto2(dict_gironi, lista_coord)))
    print('Somma delle distanze massime al quadrato per i gironi ottimali = {}'.format(CalcolaCosto_minmax(dict_gironi, lista_coord)))
    # plotto la soluzione ottimale
    PlotSolution(dict_gironi, lista_coord)
    # Stampo file con i gironi ottimi in output 
    output = open(output_name, 'w')
    for key in dict_gironi:
        output.write("{}:\n".format(key))
        girone_k = dict_gironi[key]
        for squadra in girone_k:
            output.write("\t{}\n".format(squadra[0]))
    output.close()
    
    # SOLUZIONE PROPOSTA DALLA FIT    
    # stampo a video i costi
    print('Costo totale somme distanze gironi reali = {}'.format(CalcolaCosto(dict_D1, lista_coord)))
    print('Costo totale somme distanze al quadrato gironi reali = {}'.format(CalcolaCosto2(dict_D1, lista_coord)))
    print('Somma delle distanze massime al quadrato per i gironi reali = {}'.format(CalcolaCosto_minmax(dict_D1, lista_coord)))
    # Plot soluzione FIT
    PlotSolution(dict_D1, lista_coord)