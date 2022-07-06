function [cluster,centre] = farthest_traversal(data,k)
% k-means++ algorithm (minmax)
% INPUT: - data: matrice dei punti (punto in riga)
%        - k: numero dei cluster
% OUTPUT: - cluster: vettore degli indici dei cluster associati a ogni
%           punto
%         - centre: matrice dei k centri

n = size(data,1);           % Numero totale di punti
d = size(data,2);           % Dimensione dei punti
centre_index = zeros(k,1);  % Inizializzazione del vettore degli indici dei centri
centre = zeros(k,d);        % Inizializzazione della matrice dei centri (punto in riga)

% Inserimento del primo centro scelto a caso tra i punti del dataset
centre_index(1) = randi(n); 
centre(1,:) = data(centre_index(1),:);

for j = 2:k % Ciclo su tutti i centri diversi dal primo
    % Inizializzo l'indice e la distanza dal punto più lontano dal primo
    % (j-i)-esimo centro
    farther_index = 0;
    max_dist = -1;

    for i = 1:n % Ciclo su ogni punto
        if (ismember(data(i,:),centre) == 0) % Se il punto i° non è un centro
            % nc indice del cluster più vicino al punto i°
            nc = nearest_centre(data(i,:),centre(1:j-1,:)); 
            % dist distanza tra il punto i° e l centro del cluster nc-esimo
            dist = norm(data(i,:)-data(centre_index(nc),:));
            if (dist > max_dist) % Se il punto i° è più lontano del precedente
                % Aggiornamento della distanza e dell'indice corrispondente
                max_dist = dist;
                farther_index = i;
            end
        end
    end

    % Alla fine del ciclo su tutti i punti si individua quello più lontano
    % dal centro j°, quindi lo si registra
    centre_index(j) = farther_index;
    centre(j,:) = data(farther_index,:);
end 

% Assegno a ogni punto il cluster associato al centro più vicino
cluster = nearest_centre(data,centre); 
end