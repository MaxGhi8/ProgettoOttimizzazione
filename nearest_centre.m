function cluster = nearest_centre(data,centre)
% INPUT: - data: matrice dei punti da classificare (punto in riga)
%        - centre: matrice dei centri (centro in riga)
% OUTPUT: - cluster: vettore degli indici dei cluster a cui i punti
%           appartengono (cluster(i) = j tc data(i,:) appartiene al cluster
%           j-esimo)

n = size(data,1);       % Numero di punti 
k = size(centre,1);     % Numero di centri
cluster = ones(n,1);    % Inizializzazione del vettore che indica il 
                        % cluster di appartenenza di ogni punto
for i = 1:n % Ciclo su ogni punto
    best_dist = norm(data(i,:)-centre(1,:));  % Inizializzazione della distanza
    for j = 2:k % Ciclo su ogni centro diverso da quello di riferimento
        dist = norm(data(i,:)-centre(j,:));   % Calcolo la distanza tra il 
                                              % punto fissato i e il centro 
                                              % fissato k
        if (dist < best_dist)   % Se il centro j è più vicino al centro precedente:
            best_dist = dist;   % 1.aggiorno la distanza migliore,
            cluster(i) = j;     % 2.aggiorno l'indice del cluster a j.
        end
    end
end
end

