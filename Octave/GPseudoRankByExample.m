%% GPseudoRank by example
% This vignette provides examples on
%
% 1)  how to sample distributions of orders with GPseudoRank.
%
% 2) how to preprocess a larger data set using the mini-cluster approximation
% and how to use the setDefaultParams function (example 2).
%
% 3) how to select genes for pseudotime ordering (example 3).
%
%
%% Example 1: Application to medium-sized scRNA-seq data set without mini-cluster approximation
clear all;
close all;
graphics_toolkit gnuplot
pkg load statistics
nChains = 1;% 12 are ideal for detailed convergence analysis; 1 sufficient for general 
% use
%%
fHandle          = @pseudoRank; 
fileName         = 'Shalek.csv';
%%
nSamples         = 70000;
%% 
% Informative priors of the model are necessary to ensure it concentrates 
% probability mass around the correct order and to avoid that a sampling or estimation 
% algorithm gets trapped in local modes.
%%
priorLogSDLength = 0.01;%standard deviation of prior of log-length scale, 
% recommended value
paramSamplingFreq =10; %sampling parameters for every kth order
verbose          = true; % Whether or not to print output to screen
initialise       = true;  % If you have to stop the sampler, and then 
% %wish to rerun at a later date, set this to false 
thinningFreq     = 10;     % If set to X, records every X-th sample
inputSeed        = NaN;  %use clock seed with chain-dependent offset
%% 
% delta is a concentration parameter for the distribution to choose pairs 
% of cells for moves 2 (swaps of cells close in terms of L1-distance) and move 
% 3 (reversal of the path between such two cells); the lower delta, the less the 
% concentration. Decrease delta in case of very high acceptance rates.
%%
delta            = [1/4000,1/4000];
%% 
% We apply each move of the proposal distribution with the same probability, 
% except the reversal of the entire path, which is applied with a probability 
% of 0.002.
%%
pp               = [0.2495 0.2495 0.2495 0.2495 0.002];
%% 
% Start with a random permutation:
%%
permuteData      = true;
%% 
% The times at which the cells were captured:
%%
captureTimes = [repmat(0,[1,49]),repmat(1,[1,75]),repmat(2,[1,65]),...
    repmat(4,[1 60]),repmat(6,[1 58])];
%% 
% With regInt = false, the pseudotime distances between cells are adapted 
% to different speeds of biological development over the course of the reaction.
%%
regInt           = false;%fixed recommended setting
%% 
% permutationFileName = NaN: uses a random permutation, if the name of a 
% .csv file was specified here, than the permutation specified in that file would 
% be used.
%%
permutationFileName = NaN;
%% 
% The intial step size for the GP parameters, adapted during the first 5000 
% iterations. 
%%
stepSize = [0.1,0.1];%fixed recommended setting
%% 
% We use the default settings for the proposal distribution of the moves, 
% as described in the GPseudoRank paper. 
%%
%default
n0 = max(floor(307/4),1);
        n3 = max(1,floor(307/20));
        n3a = max(3,floor(307/12));
        kk = 1;
        jj = 1;
tic
%for j = 1:nChains
%feval(fHandle, 'Shalek.csv', j, nSamples, priorLogSDLength, verbose,...
  %  initialise,thinningFreq, paramSamplingFreq,stepSize,inputSeed, ...
   % delta, pp,permuteData,captureTimes,regInt,permutationFileName,...
  %  5000,n0,n3,n3a,jj,kk);
%end
toc

%% Uncertainty of cell positions
% The uncertainty of the cell positions changes over the course of the response 
% to the infection.
% 
% First we compute the distribution of the posterior cell positions from 
% the
% 
% output files of the GPseudoRank algorithm.
%%
burnIn = 2000;%number of thinned samples to be discarded as burn-in
% 1:nChains refers to the identifieres of the MCMC chains and the initial
% permutations, the output is saved in the file 'ShalekUncertA.mat'
computePositionsPseudoRank('Shalek.csv',1:nChains,1:nChains,burnIn,...
    captureTimes,'ShalekUncertA.mat');
%% 
% Plotting the posterior distributions of the cell positions: we see that 
% the uncertainty decreases, once the reaction of the cells to the stimulant has 
% set in, and becomes much larger again towards the end of the observed pseudotime 
% range. This is relevant for downstream analyes such as the subsequent application 
% of clustering or network methods for time series data, whose accuracy may benefit 
% from focusing on the less uncertain part of the pseudotime series. 
%%
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 10)
set(0,'defaultfigurecolor',[1 1 1]);
load('ShalekUncertA.mat')
if nChains > 1
  A = mean(posFreq,3);
  xx = mean(xx,2);
else
  A = posFreq; 
 end
[sorted_xx xxI] = sort(xx);
f = figure;
set(f, "visible", "off")
colormap(hot);
imagesc(A(xxI,:)');
xlabel('cells ordered by mean position','FontSize',12);
ylabel('position (distribution)','FontSize',12);
set(gca,'YDir','normal')
ax= gca;
c=colorbar;
print("fig1", "-dpng")


%% Example 2: mini-cluster approximation and default parameter settings for the proposal moves
% Simulating a data set with 5000 cells
%%
logSw2 = randn(1,100)*0.1 + log(sqrt(1));
sw2 = exp(2*logSw2);
logL = randn(1,100).*0.1+log(0.4);
l = exp(logL);
logSe2 = randn(1,100)*0.1 + log(sqrt(0.5));
se2 = exp(2*logSe2);
ppar = [logSw2;logL;logSe2];
simGPIndiv(5000,sw2,l,se2,'sim_Example1.csv','tau_Example1.csv',false,'uniform');
%% 
% Applying the mini-cluster approximation
%%
%assuming there are capture times
captureTimes = [zeros(1,1000),ones(1,1000),repmat(2,1,1000),repmat(3,1,1000),repmat(4,1,1000)];
miniclust('sim_Example1.csv',captureTimes);
%% 
% Finding the default parameters for this data set:
%%
[n0,n3,n3a,jj,kk,delta,pp]  = setDefaultParams('sim_Example1KM.csv',true);
%%
delete 'sim_Example1*'
%% Example 3: selecting genes
% Simulating a data set with 5000 genes
%%
logSw2 = randn(1,5000)*0.1 + log(sqrt(1));
sw2 = exp(2*logSw2);
logL = randn(1,5000).*0.1+log(0.4);
l = exp(logL);
logSe2 = randn(1,5000)*0.1 + log(sqrt(0.5));
se2 = exp(2*logSe2);
ppar = [logSw2;logL;logSe2];
simGPIndiv(60,sw2,l,se2,'sim_Example1.csv','tau_Example1.csv',false,'uniform');
%% 
% Selecting genes with capture times using anova to test difference between 
% mean expressions for different capture times
%%
captureTimes = [zeros(1,20),ones(1,20),repmat(2,1,20)];
selGenes('sim_Example1.csv',captureTimes);
%% 
% Selecting genes without capture times
%%
captureTimes = [zeros(1,20),ones(1,20),repmat(2,1,20)];
selGenes('sim_Example1.csv');
%%
delete 'sim_Example1*'
        
