%illustrates the GPseudoRank algorithm on simulated data sets
%settings for plots
set(0,'DefaultAxesFontName', 'Arial')
set(0,'defaultfigurecolor',[1 1 1]);

%simulate parameters and a data set and save them
%50 genes
logSw2 = randn(1,50)*0.1 + log(sqrt(1));
sw2 = exp(2*logSw2);
logL = randn(1,50).*0.1+log(0.4);
l = exp(logL);
logSe2 = randn(1,50)*0.1 + log(sqrt(0.5));
se2 = exp(2*logSe2);
ppar = [logSw2;logL;logSe2];
csvwrite('Params_Example.csv',ppar);
%the 90 pseudotimes are simulated from a uniform distribution
simGPIndiv(90,sw2,l,se2,'sim_Example.csv','tau_Example.csv',false,'uniform');
tau = csvread('tau_Example.csv');
simData = csvread('sim_Example.csv');
%plot one gene from the simulated data set
figure();
subplot(1,2,1);
c1 = plot(tau(1:30),simData(1,1:30),'mo','MarkerFaceColor','m');
hold on;
c2 = plot(tau(31:60),simData(1,31:60),'co','MarkerFaceColor','c');
hold on;
c3 = plot(tau(61:90),simData(1,61:90),'bo','MarkerFaceColor','b');
xlabel('time','FontSize',14);
ylabel('expression level','FontSize',14);
ylim([-4,4]);
set(gca,'box','off');
ax=gca;
ax.FontSize = 14;
leg = legend([c1,c2,c3],'capture time 1','capture time 2','capture time 3');
leg.FontSize = 14;

subplot(1,2,2);
c1 = plot(ones(1,30),simData(1,1:30),'mo','MarkerFaceColor','m');
hold on;
c2 = plot(repmat(2,1,30),simData(1,31:60),'co','MarkerFaceColor','c');
hold on;
c3 = plot(repmat(3,1,30),simData(1,61:90),'bo','MarkerFaceColor','b');
xlabel('capture time','FontSize',14);
ylabel('expression level','FontSize',14);
ylim([-4,4]);
set(gca,'box','off');
ax=gca;
ax.FontSize = 14;set(gcf, 'PaperUnits', 'centimeters');
x_width=17 ;y_width=10;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
set(gcf, 'PaperSize', [x_width y_width]);
print('simData1','-dpdf','-r300');
%now run the sampler, with GP parameters assumed known, the order and
%pseudotime locations assumed unknown, the starting order is random
ppar = csvread('Params_Example.csv');
fHandle          = @pseudoRankIndivParams;  
fileName         = 'sim_Example.csv';
nSamples         = 50000;  
logPriorGPMeans = ppar;
logPriorGPSds = zeros(3,50);    
paramSamplingFreq = 10^6; %fixed parameters
verbose          = true; 
initialise       = true;  
thinningFreq     = 10;    
inputSeed        = NaN;    
delta            = [1/1000,1/1000];
permuteData      = true;
captureTimes     = [zeros(1,30),ones(1,30),repmat(2,1,30)];
regInt           = true;
stepSize = ones(3,50);
burnIn = 0;%this is required only for using the adaptive version of the proposal 
%distribution, so we set this to 0 here
pp = [0.2495 0.2495 0.2495 0.2495 0.002];
permutationFileName = NaN;
feval(fHandle, fileName, 1, nSamples, logPriorGPMeans,logPriorGPSds, ...
    verbose, initialise, thinningFreq, paramSamplingFreq,stepSize, inputSeed, delta, pp,permuteData,captureTimes,...
    regInt,permutationFileName,false,NaN,burnIn,false,false)

%%%%now repeat this with fewer capture times
figure();
subplot(1,2,1);

c1 = plot(tau(1:30),simData(1,1:30),'mo','MarkerFaceColor','m');
hold on;
c2 = plot(tau(31:90),simData(1,31:90),'co','MarkerFaceColor','c');
xlabel('time','FontSize',14);
ylabel('expression level','FontSize',14);
ylim([-4,4]);
set(gca,'box','off');
ax=gca;
ax.FontSize = 14;
leg = legend([c1,c2],'capture time 1','capture time 2');
leg.FontSize = 14;

subplot(1,2,2);
c1 = plot(ones(1,30),simData(1,1:30),'mo','MarkerFaceColor','m');
hold on;
c2 = plot(repmat(2,1,60),simData(1,31:90),'co','MarkerFaceColor','c');
xlabel('capture time','FontSize',14);
ylabel('expression level','FontSize',14);
ylim([-4,4]);
set(gca,'box','off');
ax=gca;
ax.FontSize = 14;
set(gcf, 'PaperUnits', 'centimeters');
x_width=17 ;y_width=10;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
set(gcf, 'PaperSize', [x_width y_width]);
print('simData2','-dpdf','-r300');

captureTimes     = [zeros(1,30),ones(1,60)];
%use 2 as identifier
feval(fHandle, fileName, 2, nSamples, logPriorGPMeans,logPriorGPSds, ...
    verbose, initialise, thinningFreq, paramSamplingFreq,stepSize, inputSeed, delta, pp,permuteData,captureTimes,...
    regInt,permutationFileName,false,NaN,burnIn,false,false)


%%Third simulation: more noise
%simulate parameters and a data set and save them
%50 genes
logSw2 = randn(1,50)*0.1 + log(sqrt(1));
sw2 = exp(2*logSw2);
logL = randn(1,50).*0.1+log(0.4);
l = exp(logL);
logSe2 = randn(1,50)*0.1 + log(sqrt(1));
se2 = exp(2*logSe2);
ppar = [logSw2;logL;logSe2];
csvwrite('Params_Example.csv',ppar);
%the 90 pseudotimes are simulated from a uniform distribution
simGPIndiv(90,sw2,l,se2,'sim_Example3.csv','tau_Example3.csv',false,'uniform');
tau = csvread('tau_Example3.csv');
simData = csvread('sim_Example3.csv');
%plot one gene from the simulated data set
figure();
subplot(1,2,1);
c1 = plot(tau(1:30),simData(1,1:30),'mo','MarkerFaceColor','m');
hold on;
c2 = plot(tau(31:60),simData(1,31:60),'co','MarkerFaceColor','c');
hold on;
c3 = plot(tau(61:90),simData(1,61:90),'bo','MarkerFaceColor','b');
xlabel('time','FontSize',14);
ylabel('expression level','FontSize',14);
ylim([-4,4]);
set(gca,'box','off');
ax=gca;
ax.FontSize = 14;
leg = legend([c1,c2,c3],'capture time 1','capture time 2','capture time 3');
leg.FontSize = 14;

subplot(1,2,2);

c1 = plot(ones(1,30),simData(1,1:30),'mo','MarkerFaceColor','m');
hold on;
c2 = plot(repmat(2,1,30),simData(1,31:60),'co','MarkerFaceColor','c');
hold on;
c3 = plot(repmat(3,1,30),simData(1,61:90),'bo','MarkerFaceColor','b');
xlabel('capture time','FontSize',14);
ylabel('expression level','FontSize',14);
ylim([-4,4]);
set(gca,'box','off');
ax=gca;
ax.FontSize = 14;
set(gcf, 'PaperUnits', 'centimeters');
x_width=17 ;y_width=10;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
set(gcf, 'PaperSize', [x_width y_width]);
print('simData3','-dpdf','-r300');
captureTimes     = [zeros(1,30),ones(1,30),repmat(2,1,30)];
feval(fHandle, 'sim_Example3.csv', 1, nSamples, logPriorGPMeans,logPriorGPSds, ...
    verbose, initialise, thinningFreq, paramSamplingFreq,stepSize, inputSeed, delta, pp,permuteData,captureTimes,...
    regInt,permutationFileName,false,NaN,burnIn,false,false)

delete  *_Example*%uncomment, if you want to keep the samples
