clear all
close all
clc

%% Confronto tra k-centre e k-means

filename = 'Squadre_D1_Maschile.csv';
data = readmatrix(filename,'Range','C2:D63');
M = 6;                              % Numero di squadre in ogni girone
m = size(data,1);                   % Numero di squadre totali
k = size(Tau(m,M),1);               % Numero di gironi
fprintf('Numero dei gironi: %d \n', k)

% Punti senza clustering
figure()
scatter(data(:,1), data(:,2), '.')
title('data')

%% Centri scelti a caso
centre_index = randperm(m,k);   
centre_rand = data(centre_index,:);
max_iter = 50;
plotting = false;

% Lloyd con scelta dei centri random
[cluster_rand,centre_rand,iter] = lloyd(data,centre_rand,max_iter,plotting);
[costo_min] = costi(cluster_rand, data, 'max');
[costo_quad] = costi(cluster_rand, data, 'quad');
[costo_taxi] = costi(cluster_rand, data, 'taxi');
fprintf('Costo con starter farthest_traversal: %f \n',costo_min)
fprintf('Costo con starter farthest_traversal: %f \n',costo_quad)
fprintf('Costo con starter farthest_traversal: %f \n',costo_taxi)
figure()
plot_clusters(data,cluster_rand,centre_rand)
title('Lloyd con scelta dei centri random')

fprintf('Numero di iterazioni con Lloyd e starter random: %d \n', iter)
% OSS: Questo metodo non prevede vincoli sulla cardinalità dei clusters,
% quindi gli elementi più lontani dal gruppo sono associati a cluster
% individuali


%% Centri scelti con farthest_traversal
figure()
% 1. k-centre
% Assegno a ogni punto un cluster con l'algoritmo del k-menas++ (minmax)
[clusters,centre] = farthest_traversal(data,k); 
% Plot dei clusters in colori diversi
subplot(2,2,1)
plot_clusters(data,clusters,centre);
title('k-centre')
% OSS: Similmente a prima gli elementi più lontani sono associati a cluster
% individuali

% 2. Algoritmo di Lloyd e farthest_traversal
[cluster,centre,iter] = lloyd(data,centre,max_iter,plotting);
[costo_min] = costi(cluster, data, 'max');
[costo_quad] = costi(cluster, data, 'quad');
[costo_taxi] = costi(cluster, data, 'taxi');
fprintf('Costo con starter farthest_traversal: %f \n',costo_min)
fprintf('Costo con starter farthest_traversal: %f \n',costo_quad)
fprintf('Costo con starter farthest_traversal: %f \n',costo_taxi)
fprintf('Numero di iterazioni con starter farthest_traversal: %d \n',iter)
subplot(2,2,2)
plot_clusters(data,cluster,centre);
title("Lloyd + farthest-traversal")

% 3. k-means++
[~,centre] = kmeanspp(data,k);
subplot(2,2,3)
plot_clusters(data,clusters,centre);
title('k-means++')

% 4. Algoritmo di Lloyd e k-means++
[cluster,centre,iter] = lloyd(data,centre,max_iter,plotting);
[costo_min] = costi(cluster, data, 'max');
[costo_quad] = costi(cluster, data, 'quad');
[costo_taxi] = costi(cluster, data, 'taxi');
fprintf('Costo con starter farthest_traversal: %f \n',costo_min)
fprintf('Costo con starter farthest_traversal: %f \n',costo_quad)
fprintf('Costo con starter farthest_traversal: %f \n',costo_taxi)
fprintf('Numero di iterazioni con starter k-means++: %d \n',iter)
% Cluster finale
subplot(2,2,4)
plot_clusters(data,cluster,centre);
title("Lloyd + k-means++")


