%This script repeates the inference on the Shalek data, pretending not to
%be able to distiguish the first and second capture times, to test the
%GPseudoRank method further.
maxNumCompThreads = 12;
fHandle          = @pseudoRank; 
fileName         = 'Shalek.csv';
nSamples         = 500000;  
priorLogSDLength = 0.01;%standard deviation of prior of log-length scale, recommended value
paramSamplingFreq =10; %sampling parameters for every kth order
verbose          = false; % Whether or not to print output to screen
initialise       = true;  % If you have to stop the sampler, and then wish to rerun at a later date, set this to false 
thinningFreq     = 10;     % If set to X, records every X-th sample
inputSeed        = NaN;  %use clock seed with chain-dependent offset
delta            = [1/4000,1/4000];
pp               = [0.2495 0.2495 0.2495 0.2495 0.002];
permuteData      = true;
%we assume that we cannot distinguish the first two capture times as a test
captureTimes     = [repmat(0,[1,124]),repmat(2,[1,65]),repmat(4,[1 60]),repmat(6,[1 58])];
regInt           = false;%fixed recommended setting
permutationFileName = NaN;
stepSize = [0.1,0.1];%fixed recommended setting
n0 = max(floor(307/4),1);
        n3 = max(1,floor(307/20));
        n3a = max(3,floor(307/12));
        kk = 1;
        jj = 1;

%uncomment to run the MCMC sampler
%feval(fHandle, 'Shalek.csv', 1000, nSamples, priorLogSDLength, verbose, initialise, thinningFreq, paramSamplingFreq,...
 %   stepSize,inputSeed, delta, pp,permuteData,captureTimes,regInt,permutationFileName,5000,n0,n3,n3a,jj,kk);

tic
nSamples = 50000;
feval(fHandle, 'Shalek.csv', 111, nSamples, priorLogSDLength, verbose, initialise, thinningFreq, paramSamplingFreq,...
    stepSize,inputSeed, delta, pp,permuteData,captureTimes,regInt,permutationFileName,5000,n0,n3,n3a,jj,kk);
toc
%Elapsed time is 254.386482 seconds.
%compute the cell position frequencies
burnIn = 5000;%number of thinned samples to be discarded as burn-in
%Run the following lines only after running the MCMC sampler.
%computePositionsPseudoRank('Shalek.csv',1000,1000,burnIn,captureTimes,'ShalekUncertB.mat');

%%%%%%%
%plot mean and sd of position of each cell
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 8)
set(0,'defaultfigurecolor',[1 1 1]);
set(0, 'DefaultLineLineWidth', 0.7);
%for the plot, we colour the cells by their actually known capture times
captureTimes = [repmat(1,[1,49]),repmat(2,[1,75]),repmat(3,[1,65]),repmat(4,[1 60]),repmat(5,[1 58])];
%we use the true known capture times here, to identify them in the plot,
%not for computation
%convert positions to pseudotimes
[x1, ~] = rankToPseudo('Shalek_Results_Orders_Chain1000.csv','Shalek.csv',...
    'Shalek_Permutation1000.csv',5000,captureTimes);
xx = mean(x1,1);
yy = sqrt(var(x1));
figure()
a1=plot(xx(250:307),yy(250:307),'ko','MarkerSize',4);
hold on;
a2=plot(xx(190:249),yy(190:249),'ro','MarkerSize',4);
hold on;
a3=plot(xx(125:189),yy(125:189),'bo','MarkerSize',4);
hold on;
a4=plot(xx(50:124),yy(50:124),'go','MarkerSize',4);
hold on;
a5=plot(xx(1:49),yy(1:49),'mo','MarkerSize',4);
xlim([0,1]);
set(gca,'box','off');
ax=gca;
xlabel('mean pseudotime','FontSize',10);
ylabel('sd of pseudotime','FontSize',10);
leg=legend([a5,a4,a3,a2,a1],'0h', '1h', '2h', '4h','6h','Orientation','horizontal',...
    'Location','NorthOutSide');
leg.FontSize = 8;
legend boxoff 
set(gcf, 'PaperUnits', 'centimeters');
x_width = 8.8;
y_width = 5;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
 set(gcf,'PaperSize',[x_width y_width]);
print('ShalekPosMeanSdB','-dpdf','-r600');
%
