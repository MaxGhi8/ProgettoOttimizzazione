function [cluster,centre,iter] = lloyd(data,centre,max_iter,plotting)
% Lloyd's Algorithm for the k-means clustering
% INPUT: - data matrice dei punti (punti in riga), 
%        - centre matrice dei centri iniziali,
%        - max_iter massimo numero di iterazioni
%        - plotting può essere 'true' (quindi si plotta il clustering dopo 
%          6 iterazioni) oppure 'false (non crea alcun disegno)
% OUTPUT: - cluster vettore di indici per i punti 
%         - centre matrice dei centri dei clusters
%         - iter numero di iterazioni

n = size(data,1);       % Numero dei punti 
k = size(centre,1);     % Numero dei clusters
iter = 0;               % Contatore delle iterazioni
convergence = false;    % 'true' quando non ci sono stati aggiornamenti

% Disegna a seconda delle istruzioni date in input
if (plotting)
    figure
end

% Eseguire finché ci sono aggiornamenti e non superiamo max_iter
while (~convergence && iter < max_iter)
    % Vettore degli indici dei clusters
    cluster = nearest_centre(data,centre); 
    [centroid,empty] = cluster_centroid(data,cluster);
    for i = 1:k % Ciclo sui clusters
        if (empty(i) == 1)          % Se l'i° cluster è vuoto
            % Assegno il nuovo centro a caso
            centroid(i,:) = data(randi(n),:);
            % How can we improve this?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
    if (norm(centre-centroid) == 0) % Se c'è convergenza
        convergence = true;         % Segno convergence
    else
        centre = centroid;          % Altrimenti aggiorno il centro
    end
    iter = iter + 1;    % Conto l'iterazione
    
     % Plot delle prime 6 iterazioni
    if (iter <= 6 && plotting)
        subplot(2,3,iter)
        plot_clusters(data,cluster,centre)
        title(['k=', num2str(k), ' clusters'],[num2str(iter), ' iterazioni'])
    end
end
end

