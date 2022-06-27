function cluster = nearest_centre(data,centre)
% For each point in the rows of data return the index of the nearest centre
% (both data and centre points are in the rows of their matrix)
n = size(data,1);
k = size(centre,1);
cluster = ones(n,1); % initialisation
for i=1:n
    best_dist = norm(data(i,:)-centre(1,:)); % initialisation
    for j=2:k
        dist = norm(data(i,:)-centre(j,:));
        if(dist < best_dist) % update best centre
            best_dist = dist;
            cluster(i) = j;
        end
    end
end
end

