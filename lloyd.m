function [cluster,centre,iter] = lloyd(data,centre,max_iter,plotting)
% Lloyd's Algorithm for the k-means clustering
% INPUT: data points, initial centres, max number of iterations, plotting
% can be set to 'true' to plot the clustering after the first 6 iterations
% OUTPUT: cluster indices for data points, clusters' centres, n. iterations
n = size(data,1); % number of data points
k = size(centre,1); % number of clusters
iter = 0; % number of iteration performed so far
convergence = false; % 'true' when no improvement has been obtained
if(plotting)
    figure
end
while(~convergence && iter < max_iter)
    cluster = nearest_centre(data,centre);
    [centroid,empty] = cluster_centroid(data,cluster);
    for i=1:k
        if(empty(i)==1)
            % random selection of a new centre
            centroid(i,:) = data(randi(n),:);
            % How can we improve this?
        end
    end
    if(norm(centre-centroid)==0) % convergence has been reached
        convergence = true;
    else
        centre = centroid;
    end
    iter = iter + 1;
    % plot the first 6 iterations
    if(iter<=6 && plotting)
        subplot(2,3,iter)
        plot_clusters(data,cluster,centre)
    end
end
end

