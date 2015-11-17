
function [count]= binding_simulation_constNoise(nC, nR, fig);

%% Select Parameters

a = 0;      % lower XY limit
b = 1;      % upper XY limit

% Number of clusters, cluster density (cluster/?m2)

% nC = 10;        % Number of clusters
% nR = 100;       % Number of receptors per cluster
nTot = 2000;      % Total number of receptors

%% Generate Clusters without noise --> store in variable clusters


% Generate Cluster centers

centersx = (b-a).*rand(nC,1) + a; % generate centers
centersy = (b-a).*rand(nC,1) + a; % generate centers

% Generate Clusters

clusterx=[];
clustery=[];

for count=1:length(centersx);
   
    clustercx = (randn(nR,1)/50)+centersx(count);
    clustercy = (randn(nR,1)/50)+centersy(count);
    clusterx=[clusterx; clustercx];
    clustery=[clustery; clustercy];
    clear clustercx cluster cy
    count=count+1;
   
end

% Generate List of remaining unclustered receptors

Rec(:,1) = (b-a).*rand(1000,1) + a;
Rec(:,2) = (b-a).*rand(1000,1) + a;

clusters(:,1)=clusterx;
clusters(:,2)=clustery;

% Combine all receptors 

allRec=[Rec; clusters];

% Plot cluster centers and clusters

if fig==1;

figure('Position',[100 500 900 250])

subplot(1,3,1)
scatter(centersx, centersy); hold on;
title('Cluster centers')
axis([a b a b]);

subplot(1,3,2)
scatter(clusterx(:), clustery(:),1);hold on;
scatter(Rec(:,1), Rec(:,2),1);
axis([a b a b]);
title('Clusters')

else end

%% Select a random point and find Localization around a cricle 

rV = 0.05; % Radius Virus

% Generate location of the virus X,Y

locV(:,1) = ((b-rV)-a).*rand(1,1) + a;
locV(:,2) = ((b-rV)-a).*rand(1,1) + a;

% Draw circle around the virus

if fig==1;

subplot(1,3,3)
scatter(clusterx(:), clustery(:),1);hold on;
scatter(Rec(:,1), Rec(:,2),1);
scatter(locV(:,1),locV(:,2),'filled');
viscircles(locV,rV)
axis([a b a b]);
title('Clusters')

else end

count=0;

for i=1:length(allRec);
    
    if  sqrt(((allRec(i,1)-locV(1,1))^2)+((allRec(i,2)-locV(1,2))^2)) <= rV;
        
        count=count+1;
        
    else 
        
    end
    
end

end

