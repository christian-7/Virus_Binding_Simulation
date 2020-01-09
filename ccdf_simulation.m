%% CCDF for having a given number of receptors at the pathogen-cell interface

clear, clc, close all

nR      = (0:1:100); % number of receptors per cluster
nTot    = 3000;      % Total number of receptors
nC      = 5;         % number of clusters
iter    = 200;       % number of iterations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

successRate = {}; p = 0; rng shuffle

for n = [10, 20, 50, 80, 100]; % minimum number of receptors  
    
p = p + 1;

l = 1; aver = [0,0];
    
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


%% ccdf for having a given number of receptors at the pathogen-cell interface for different cluster densities

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
