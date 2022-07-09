function [costo] = costi(cluster, data, flag)
% calcola il costo di ogni girone
% INPUT: - cluster: matrice degli indici dei cluster dei punti
%        - data: matrice dei punti 2xn
%        - flag: indica quale metodo di calcolo del costo (max, quad, taxi)
% OUTPUT: costo: numero reale del costo della clusterizzazione

n = size(data, 1);
k = max(cluster);
costo = 0;

if strcmp(flag,'max')
    for i = 1:k
        dist_max = 0;
        for j = 1:n
            for m = 1:j-1
                if cluster(j) == i && cluster(m) == i
                    dist = norm(data(j,:)-data(m,:))^2;
                    dist_max = max(dist_max, dist);
                end
            end
        end
        costo = costo + dist_max;
    end
end

if strcmp(flag,'quad')
    for i = 1:n
        for j = 1:i-1
            if cluster(i) == cluster(j)
                dist = norm(data(i,:)-data(j,:))^2;
                costo = costo + dist;
            end
        end
    end
end

if strcmp(flag,'taxi')
    for i = 1:n
        for j = 1:i-1
            if cluster(i) == cluster(j)
                dist = norm(data(i,:)-data(j,:));
                costo = costo + dist;
            end
        end
    end
end

end