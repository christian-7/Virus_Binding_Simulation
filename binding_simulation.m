%% 2D Particle binding simulation
% Author: Christian Sieben, EPFL 
% sieben.christian@gmail.com
% April 2020

% Simulation of a 2D surface (1x1µm) with
% 
% nTot    = total number of receptors
% nC      = number of clusters
% fract   = fraction of receptors in clusters
% fig     = make figure yes/no

% fract defines how many receptors are free !!
% receptors are gradually shifted into clusters

function [count]= binding_simulation(nTot, nC, fract, fig);

%% Select Parameters

% clear, clc, close all

a = 0;          % lower XY limit, nm
b = 1000;       % upper XY limit, nm

% Number of clusters, cluster density (cluster/um2)

% fig     = 1;
% nC      = 20;        % Number of clusters
% fract   = 0.5;       % Number of receptors per cluster
% nTot    = 2000;      % Total number of receptors

% Generate random cluster centers
%% Generate Clusters without noise --> store in variable clusters

clusters = [];

clusters(:,1) = a+(b-a).*rand(nC,1); % generate cluster centers X
clusters(:,2) = a+(b-a).*rand(nC,1); % generate cluster centers Y

% Generate random cluster radius

for i = 1:nC
    
    clusters(i,3) = 40+(150-40).*rand(1,1); % cluster radius in nm from STORM
    
end

% Find number of AFs per cluster

for i = 1:nC
    
clusters = [];

clusters(:,1) = a+(b-a).*rand(nC,1); % generate centers
clusters(:,2) = a+(b-a).*rand(nC,1); % generate centers

% Find cluster radius

for i = 1:nC
    
    clusters(i,3) = 40+(150-40).*rand(1,1); % cluster radius in nm from STORM
    
end

% Find number of rec per cluster

for i = 1:nC
    
   clusters(i,4) = (clusters(i,3)/(sum(clusters(:,3))));
   clusters(i,5) = round(clusters(i,4)*(fract*nTot)); 
   
end

% Generate Clusters

clusterx=[]; clustercx = [];
clustery=[]; clustercx = [];

for i=1:nC;
   

%     clustercx = normrnd(clusters(i,1),clusters(i,3)/3,[clusters(i,5),1]); % locs as normal distribution
%     clustercy = normrnd(clusters(i,2),clusters(i,3)/3,[clusters(i,5),1]); % locs as normal distribution
        
    tempX       = clusters(i,3)*sqrt(rand(clusters(i,5),1));
    tempY       = 2*pi*rand(clusters(i,5),1);
    
    clustercx   = tempX.*cos(tempY)+clusters(i,1);
    clustercy   = tempX.*sin(tempY)+clusters(i,2);
    
    clusterx  = [clusterx; clustercx];
    clustery  = [clustery; clustercy];
    
    clear clustercx clustercy tempX tempY  
    clustercx = normrnd(clusters(i,1),clusters(i,3)/3,[clusters(i,5),1]);
    clustercy = normrnd(clusters(i,2),clusters(i,3)/3,[clusters(i,5),1]);
    clusterx  = [clusterx; clustercx];
    clustery  = [clustery; clustercy];
    clear clustercx cluster 

    
end

% Generate List of remaining unclustered receptors

Rec      = []; rec_clusters = [];
Rec(:,1) = (b-a).*rand(nTot-(fract*nTot),1) + a;
Rec(:,2) = (b-a).*rand(nTot-(fract*nTot),1) + a;

rec_clusters(:,1) = clusterx;
rec_clusters(:,2) = clustery;

rec_clusters(:,1)=clusterx;
rec_clusters(:,2)=clustery;


% Combine all receptors 

allRec = [Rec; rec_clusters];

% Plot cluster centers and clusters

if fig==1;

figure('Position',[100 500 800 250])

subplot(1,3,1)
scatter(clusters(:,1), clusters(:,2)); hold on;
title('Cluster centers')
box on
axis([a b a b]);
axis square
xlabel('x [nm]');ylabel('y [nm]');

subplot(1,3,2)
scatter(clusterx(:), clustery(:),1);hold on;
scatter(Rec(:,1), Rec(:,2),1);
axis([a b a b]);
box on
title('Clusters ')

axis square
xlabel('x [nm]');ylabel('y [nm]');

else end
 
%% Select a random point and find receptors around a cricle 

rV = 50; % Radius Virus, nm

% Generate location of the virus X,Y

locV(:,1) = ((b-rV)-a).*rand(1,1) + a;
locV(:,2) = ((b-rV)-a).*rand(1,1) + a;

% Draw circle around the virus

if fig==1;

subplot(1,3,3)
scatter(clusterx(:), clustery(:),1);hold on;
scatter(Rec(:,1), Rec(:,2),1);
scatter(locV(:,1),locV(:,2),'filled');
viscircles(locV,rV);
axis([a b a b]);
title('Clusters with virus binding site')
axis square
xlabel('x [nm]');ylabel('y [nm]');box on;

else end

count=0;

for i=1:length(allRec);
    
    if  sqrt(((allRec(i,1)-locV(1,1))^2)+((allRec(i,2)-locV(1,2))^2)) <= rV;
        
        count=count+1;
        
    else 
        
    end
    
end

end
