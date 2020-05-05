%% 2D Particle binding simulation

% Author: Christian Sieben, EPFL 
% sieben.christian@gmail.com
% April 2020

% Simulate particle binding attempts using 
% binding_simulation and binding_simulation_constNoise
% Plot CCDF curves and binding probability  

%% CCDF for binding a given number of receptors at the pathogen-cell interface

clear, clc, close all

nR      = (0:1:50);  % number of receptors per cluster  --> (locs per cluster)
nTot    = 100;       % Total number of receptors        --> (molecules per µm2)
nC      = 2;         % number of clusters               --> (cluster per µm2)
iter    = 1000;      % number of iterations

%% CCDF for having a given number of receptors at the pathogen-cell interface

clear, clc, close all

nR      = (0:1:100); % number of receptors per cluster
nTot    = 3000;      % Total number of receptors
nC      = 5;         % number of clusters
iter    = 200;       % number of iterations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

successRate = {}; p = 0; rng shuffle

for n = [1, 5, 10, 20, 50]; % minimum number of receptors  

for n = [10, 20, 50, 80, 100]; % minimum number of receptors  

    
p = p + 1;

l = 1; aver = [0,0];
    

for m = [0,0.1,0.2,0.4,0.8,1] %                     --> (% in cluster)
    
count=[]; success = 0; 

for k=1:iter; % iterate of each degree of clustering
    
    count(k,1) = binding_simulation(nTot, nC,m,0); % # of receptors, # of clusters, degree of clustering, figure 
    
    if count(k,1) > n;
        success = success+1; % count as successful attempt
    else end
                                             
end

pd          = fitdist(count,'Normal'); % distribution of bound receptors per attempt
ccdf(:,l)   = 1-cdf(pd,nR); % complementary CDF 

aver(l,1) = m;                  % in cluster
aver(l,2) = mean(count);        % mean bound receptors
aver(l,3) = std(count);         % std bound receptors
aver(l,4) = (success/iter)*100; % success rate in %
l         = l+1;

clc, n

end

successRate{p, 1} = n;          % min number of receptors
successRate{p, 2} = aver;   
successRate{p, 3} = ccdf;
=======
for m = [0,0.1,0.2,0.4,0.8,1] % --> % clustered
    
count=[]; success = 0; 

for k=1:iter
    
    count(k,1)=binding_simulation(nTot, nC,m,0); % Number of clusters 
                                           % Number of receptors per cluster
                                           % Figure  
    if count(k,1) > n;
        success = success+1;
    else end
                                             
end

pd          = fitdist(count,'Normal'); % distribution of receptors bound per attempt
ccdf(:,l)   = 1-cdf(pd,nR); % complementary CDF 

aver(l,1) = m;
aver(l,2) = mean(count);
aver(l,3) = std(count);
aver(l,4) = (success/iter)*100; % success rate in %
l         = l+1;

clc, n

end

successRate{p, 1} = n;
successRate{p, 2} = aver;
successRate{p, 3} = ccdf;

end

display('Done');

figure('Position',[400 50 600 250])
subplot(1,2,1)
>>>>>>> da30444d5e56e8cfc95dd0811700f23f6904b9af

stairs(successRate{4, 3},'LineWidth',1);hold on;
line([50 50], [1e-3 1e0], 'Color','black','LineWidth',2);
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
axis([1e0 1e2 1e-3 1e0]);
% title('CCDF for having n bound receptors');
xlabel('number of receptors [n]');
ylabel('binding probability [%]'); % probablity to bind n receptors
box on
leg = legend('show');
title(leg,'Degree of clustering %')
legend('random','10','20','40','80','100','Location','southwest');

subplot(1,2,2)
% scatter(aver(:,1),aver(:,2)); hold on; 
% errorbar(aver(:,1),aver(:,2),aver(:,3));
% title('Mean number of bound receptors');
% ylabel('mean bound receptors');
% xlabel('degree of clustering');
% box on
% 
% plot(aver(:,1),aver(:,4),'-o'); hold on; 
% ylabel('binding success rate (%)');
% xlabel('degree of clustering (%)');
% box on  
% axis square
% leg = legend('show');
% title(leg,'Minimum # of receptors')
% legend(num2str(n));
% box on

% figure
for i = 1:size(successRate,1);
   
plot(successRate{i, 2}(:,1),successRate{i, 2}(:,4),'-o','LineWidth',1); hold on; 
ylabel('binding success rate [%]');
xlabel('degree of clustering [%]');
box on  
axis square
leg = legend('show');
title(leg,'Minimum # of receptors')
legend('10','20', '50', '80', '100');
    
end

<<<<<<< HEAD
display('Done');

figure('Position',[400 50 600 250])
subplot(1,2,1)

stairs(successRate{4, 3},'LineWidth',1);hold on;
line([50 50], [1e-3 1e0], 'Color','black','LineWidth',2);
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
axis([1e0 1e2 1e-3 1e0]);
% title('CCDF for having n bound receptors');
xlabel('number of receptors [n]');
ylabel('binding probability [%]'); % probablity to bind n receptors
box on
leg = legend('show');
title(leg,'Degree of clustering %')
legend('random','10','20','40','80','100','Location','southwest');

subplot(1,2,2)
scatter(aver(:,1),aver(:,2)); hold on; 
errorbar(aver(:,1),aver(:,2),aver(:,3));
title('Mean number of bound receptors');
ylabel('mean bound receptors');
xlabel('degree of clustering');
box on

plot(aver(:,1),aver(:,4),'-o'); hold on; 
ylabel('binding success rate (%)');
xlabel('degree of clustering (%)');
box on  
axis square
leg = legend('show');
title(leg,'Minimum # of receptors')
legend(num2str(n));
box on

figure('Position',[400 50 300 300])

for i = 1:size(successRate,1);
   
plot(successRate{i, 2}(:,1),successRate{i, 2}(:,4),'-o','LineWidth',1); hold on; 
ylabel('binding success rate [%]');
xlabel('degree of clustering [%]');
box on  
axis square
leg = legend('show');
title(leg,'Minimum # of receptors')
legend('1','5', '10', '20', '50');
    
end

%% Iterate over spatial parameters (cluster size, density)

Steps   = 10;
iter    = 100;

nR = linspace(1,100,Steps)';
cR = linspace(1,100,Steps)';

[p,q] = meshgrid(nR, cR);
a3 = [p(:) q(:)];

for i = 1:length(a3);   % choose the parameter combination

    success = 0; 
    
    for k = 1:iter;     % run for n iterations
    
    count(k,1) = binding_simulation_constNoise(10, a3(i,1), a3(i,2), 0);
    
    if count(k,1) > 10;
        success = success+1; % count as successful attempt
    else end
    
    end
    
    a3(i,3) = mean(count);
    a3(i,4) = (success/iter)*100; % binding probability 
    clc;
    display(['Done with: ' num2str(a3(i,1:2))])
    
end

display('Done');

%% Plot the result


% x = a3(:,1)./(pi*(a3(:,2).^2)); % receptors/cluster area (nm2)
x = a3(:,1);
y = a3(:,2);
z = a3(:,4);

  xi = linspace(min(x),max(x),100);
  yi = linspace(min(y),max(y),100);
  [XI YI]=meshgrid(xi,yi);
  ZI = griddata(x,y,z,XI,YI);
 
% cmap = jet(32); % or 256, 64, 32 or whatever.

figure('Position',[400 50 300 300])
contourf(XI,YI,ZI,3) 
axis square
h = colorbar();
ylabel(h,'Binding probability [%]','FontSize', 12)

xlabel('receptors per cluster','FontSize', 12)
ylabel('cluster radius [nm]','FontSize', 12)



%% 

figure('Position',[400 50 600 250])

subplot(1,2,1)

for i = 1:size(successRate,1);

aver    = [0,0];
l       = 1;
nR      = [0:1:100];
p       = 0;
BindingWNoise = {};

for m=[1,10,20,50,100]; % number of receptors per cluster
    
count=[]; p = p + 1;

for k=1:200;
    
    count(k,1) = binding_simulation_constNoise(20,m,0);
    
end

pd      = fitdist(count,'Normal');
ccdf    = 1-cdf(pd,nR);

aver(l,1) = m;
aver(l,2) = mean(count);
aver(l,3) = std(count);
l         = l+1;

BindingWNoise{p, 1} = m;
BindingWNoise{p, 2} = aver;
BindingWNoise{p, 3} = ccdf;

end

display('Done');
%% 

figure('Position',[400 50 600 250])

subplot(1,2,1)

for i = 1:size(successRate,1);

stairs(BindingWNoise{i, 3},'LineWidth',1);hold on;
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
axis([1e0 1e2 1e-2 1e0]);
% title('CCDF for having n bound receptors');
xlabel('number of receptors [n]');
ylabel('Probability [%]');
legend('random','10 rec/cluster','20 rec/cluster','50 rec/cluster','100 rec/cluster','Location','southwest')
box on
end
line([50 50], [1e-3 1e0], 'Color','black','LineWidth',2);

subplot(1,2,2)
plot(BindingWNoise{5, 2}(:,1),BindingWNoise{5, 2}(:,2),'-o','LineWidth',1, 'Color','black'); hold on; 
% title('Mean number of bound receptors');
ylabel('mean bound receptors');
xlabel('receptors per cluster');
box on

