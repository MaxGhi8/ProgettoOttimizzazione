clear all 
close all
clc

%%
filename = 'Squadre_D1_Maschile.csv';
data = readmatrix(filename,'Range','C2:D63');

% format long
% data{1,3}

M = 6;
m = size(data,1);
k = size(Tau(m,M),1);
centre_index = randperm(m,k);
centre = data(centre_index,:);
max_iter = 50;
plotting = true;
[cluster,centre,iter] = lloyd(data,centre,max_iter,plotting);