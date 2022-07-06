function [cluster,centre] = kmeanspp(data,k)
% INPUT: - data: matrice di punti (punto in riga)
%        - k: numero di clusters
% OUTPUT: - cluster: vettore degli indici del cluster a cui è associato
%           ogni punto
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
    weight = zeros(n,1);    % Vettore che conterrà le distanze di ogni punto  
                            % dal centro più vicino
    for i = 1:n % Ciclo su ogni punto
        if (ismember(data(i,:),centre) == 0) % Se l'i° punto non è un centro
            % Calcolo la distanza tra l'i° pto e il centro di riferimento
            weight(i) = norm(data(i,:)-data(centre_index(1),:))^2;
            for m = 2:j-1 % Ciclo su tutti i centri precenti al j°
                % Prendo la distanza minima tra i centri finora considerati
                weight(i) = min(weight(i), norm(data(i,:)-data(centre_index(m),:))^2);
            end
        else % Se l'i° pto è un centro non lo considero
            weight(i) = 0;
        end
    end

    % Alla fine del ciclo su tutti i punti si individua quello più lontano
    % dal centro j°, quindi lo si registra
    sum_weight = sum(weight);
    prob = weight / sum_weight;
    r = rand();
    sum_prob = prob(1);
    data_index = 1;
    while (r > sum_prob)
        data_index = data_index + 1;
        sum_prob = sum_prob + prob(data_index);
    end
    centre_index(j) = data_index;
    centre(j,:) = data(data_index,:);
end 

% Assegno a ogni punto il cluster associato al centro più vicino
cluster = nearest_centre(data,centre); 
end