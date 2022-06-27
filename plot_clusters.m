function plot_clusters(data,cluster,centre)
hold on
%first 10 colours are pre-selected
colours = [1 0 0; 0 1 1; 0 1 0; 1 0 1; 0 0 1; 1 1 0; .7 .7 .7;...
    .2 .5 .1; .5 .1 .8; .8 .5 .4];
for i=1:size(cluster,1)
    if i<=10
        c = colours(i,:);
    else %for k > 10 colours are randomly generated
        r = rand();
        g = rand();
        b = rand();
        c = [r g b];
    end
    plot(data(cluster==i,1),data(cluster==i,2),'*','Color',c,'MarkerSize',3)
end
%plot the centres
plot(centre(:,1),centre(:,2),'k+','LineWidth',3,'MarkerSize',12)
end