function [centroid,empty] = cluster_centroid(data,cluster)
% INPUT: - data matrice dei punti (punto in riga)
%        - cluster vettore degli indici dei cluster corripondenti ai punti
% OUTPUT: - centroid matrice dei centroidi dei clusters 
%         - empty vettore-indicatore the cluster vuoti (ha in entrata 1 in 
%           corrispondenza dei clusters vuoti)

% Assunzione: i clusters sono numerati da 1 a k
% Inizializzazione
k = max(cluster);                   % Numero dei clusters
n = size(data,1);                   % Numero dei punti
sum_data = zeros(k,size(data,2));   % Ogni riga rappresenterà la somma dei  
                                    % punti nello stesso cluster
n_data = zeros(k,size(data,2));     % Ogni riga rappresenterà il numero di 
                                    % punti nel cluster
% Completamento
for i = 1:n % Ciclo su ogni punto
    % Sommo il punto i° nella riga del suo cluster 
    sum_data(cluster(i),:) = sum_data(cluster(i),:) + data(i,:);
    % Conto il punto i° nella cardinalità del suo cluster
    n_data(cluster(i)) = n_data(cluster(i)) + 1;
end

% Inizializzazione
centroid = zeros(k,size(data,2));   % Matrice dei centroidi
empty = zeros(k,1);                 % Vettore-indicatore dei clusters vuoti

% Completamento
for j = 1:k % Ciclo sui cluster
    if (n_data(j) > 0)                           % Se il cluster j° non è vuoto
        centroid(j,:) = sum_data(j,:)/n_data(j); % Calcolo il centroide dei 
                                                 % punti del cluster j
    else
        empty(j) = 1;                            % Segno la riga j-esima
    end
end
end

