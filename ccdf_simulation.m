%% CCDF for having a given number of receptors at the pathogen-cell interface
% nR = number of receptors per cluster

clear, clc, close all, figure('Position',[400 50 800 600])
aver=[0,0];
l=1;
nR=[0:1:50];


for m=[1,10,20,40];
    
count=[];

for k=1:100;
    
    count(k,1)=binding_simulation(50,m,0);
    
end

% ecdf(count);hold on;  
%  ccdf=1-ecdf(count);

pd = fitdist(count,'Normal');
ccdf =1-cdf(pd,nR);

subplot(2,2,1)
% loglog(ccdf); hold on;
stairs(ccdf);hold on;
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
axis([1e0 1e2 1e-2 1e0]);
title('CCDF for having n bound receptors');
xlabel('number of bound receptors [n]');
ylabel('Probability [%]');

legend('random','25 %','50 %','100 %','Location','southwest','FontSize',14)


aver(l,1)=m;
aver(l,2)=mean(count);
l=l+1;

end

subplot(2,2,2)
scatter(aver(:,1),aver(:,2))
title('Mean number of bound receptors');
ylabel('mean bound receptors');
xlabel('receptors per cluster');

%% ccdf for having a given number of receptors at the pathogen-cell interface for different cluster densities

aver=[0,0];
l=1;
nR=[0:1:50];

for m=[1,5,10,20,50];
    
count=[];

for k=1:100;
    
    count(k,1)=binding_simulation_constNoise(20,m,0);
    
end

% ecdf(count);hold on;  
% ccdf=1-ecdf(count);

pd = fitdist(count,'Normal');
ccdf = 1-cdf(pd,nR);

subplot(2,2,3)
% loglog(ccdf); hold on;
stairs(ccdf);hold on;
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
axis([1e0 1e2 1e-2 1e0]);
title('CCDF for having n bound receptors');
xlabel('number of receptors [n]');
ylabel('Probability [%]');
legend('random','5 rec/cluster','10 rec/cluster','20 rec/cluster','50 rec/cluster','Location','southwest','FontSize',14)


aver(l,1)=m;
aver(l,2)=mean(count);
l=l+1;

end

subplot(2,2,4)
scatter(aver(:,1),aver(:,2))
title('Mean number of bound receptors');
ylabel('mean bound receptors');
xlabel('receptors per cluster');




