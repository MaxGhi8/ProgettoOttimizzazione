function [centroid,empty] = cluster_centroid(data,cluster)
% Given a set of data (in rows) and their cluster indices, this function
% returns the centroids of the clusters and the vector empty that is 1 in
% the positions of the empty cluster (when they exist)

k = max(cluster); % we assume clusters are enumerated from 1 to k
n = size(data,1); % number of data points
sum_data = zeros(k,size(data,2)); % each rows is the sum of the points of the
                                  % cluster
n_data = zeros(k,size(data,2)); % each rows is the number of points of the
                                % cluster
for i = 1:n
    sum_data(cluster(i),:) = sum_data(cluster(i),:) + data(i,:);
    n_data(cluster(i)) = n_data(cluster(i)) + 1;
end
centroid = zeros(k,size(data,2));
empty = zeros(k,1); % vector that indicates the empty clusters
for j=1:k
    if n_data(j)>0
       centroid(j,:) = sum_data(j,:)/n_data(j); % compute the average
    else
        empty(j) = 1; % j-th cluster is empty
    end
end
end

