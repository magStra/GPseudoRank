function [] = pseudoRank(fileName, uniqueIdentifier, nSamples, logPriorSDLength, ...
    verbose, initialise, thinningFreq, paramSamplingFreq,stepSize, inputSeed, delta, pp,permuteData,captureTimes,...
    regularInt,permutationFileName,burnInParams,n0,n3,n3a,jj,kk)
%input: 
%fileName: name of the .csv file containing the mRNA-expression matrix
%uniqueIdentifier: positive integer to identify the chain
%nSamples: number of samples to be drawn from the posterior distribution
%logPriorSDLength: standard deviation of the prior distribution of the 
%log length scale of the GP
%verbose: if true, prints acceptance rates at every 1000th iteration
%initialise: If you want to start the sampler from the beginning, 
%set this to true. If you have to stop the sampler, and then wish to rerun 
%at a later date, set this to false. 
%thinningFreq: if thinningFreq = k, we save every kth sample in the csv
%output files
%paramSamplingFreq: how often we sample the Gaussian process parameters 
%(paramSamplingFreq = k means we sample them at every kth iteration).
%stepSize: initial step size for the GP parameters, adapted during the
%first burnInParams iterations
%inputSeed: set for the random number generator, if set to NaN it is
%clock-seeded with a chain-dependent offset
%delta: vector of two hyperparameters for moves 2 and 3 of the proposal 
%distribution for the orders: 
%swaps of elements with short L1-distances, and reversals of the original order
%between those elements
%the higher delta, the less the distributions decrease with the L1-distances 
%between the cells
%pp: stochastic vector of length 5, with the probabilities of the
%individual moves being chosen
%permuteData: should generally be set to true to make the chain start from
%a random order of the cells
%captureTimes: the capture times of the cells on which the sampler is run
%regularInt: if true, the pseudotime intervals are assumed to be regular
%(faster, but this does not account for different speeds of biological
%development during the process) or irregular (approximating pseudotimes as
%sums of local Euclidean distances)
%permutationFileName: set to NaN, if you want to start from a randomly
%permutated order, to start from a particular order specified in
%'file.csv', set permutationFileName = 'file.csv'
%burnInParams: number of iterations at the beginning of the chain during
%which the proposal distribution for the sampling of the GP parameters is
%adapted
%n0, n3,n3a,kk,jj parameters for the proposal distribution (see
%documentation and )

%%%%%%%
%output files:
%Four .csv files: FileName refers to the fileName without the .csv
%extension:
    %1) FileName_PermutationuniqueIdentifier.csv: contains the permutation
    %of the original order of the cells that was used as the starting order
    %2) FileName_Results_ChainuniqueIdentifier.csv: contains the samples
    %the GP parameters and number of sampled and accepted orders for each
    %individual move of the orders, and the number of proposed and accepted
    %GP parameters
    %hyperparameters
    %3) FileName_Results_Orders_ChainuniqueIdentifier.csv: the samples of
    %the cell orders
    
if(isnan(inputSeed))
    inputSeed = sum(100*clock) + sum(100*double(uniqueIdentifier)); %%clock-seed the random number generator (with a chain-depenedent offset)
end

rng(inputSeed);

saveFileNameOrder = [strtok(fileName, '.'),'_Results_Orders_Chain', num2str(uniqueIdentifier)];
saveFileName      = [strtok(fileName, '.'),'_Results_Chain', num2str(uniqueIdentifier)];
%Lk, parameters, acceptance rates
allData = importdata(fileName, ',');
try
    data    = allData.data;
catch
    data    = allData;%if there are no gene names and cell names
end
data = (data-mean(data(:)));
[nGenes nFeatures]  = size(data);
varData             = var(data(:));
N = nGenes*nFeatures;
outFileOrder    = [saveFileNameOrder '.csv'];
outFile         = [saveFileName '.csv'];
nOrderProposals     = [0 0 0 0 0];
nOrderAcceptances   = [0 0 0 0 0];
nParamProposals     = [0 0];
nParamAcceptances  = [0 0];
if(~initialise)
    % Set up from the csv files
    B = csvread([saveFileNameOrder,'.csv'],0,0);
    currentOrder = B(end,:);
    initPerm = csvread(permutationFileName);
    data = data(:,initPerm);
    permData = data(:,currentOrder);
    if regularInt == false
        tau1 = sqrt(sum((diff(permData,1,2).^2),1));
        pseudoTimes = cumsum([0 tau1])/sum(tau1);
    else
       pseudoTimes = (0.5:1:(nFeatures-0.5))/nFeatures;
    end
    C = csvread([saveFileName,'.csv']);
    
    nOrderAcceptances   = C(end,8:12);
    nOrderProposals     = C(end,13:17);
    currentLogParams    = C(end,2:3)';
    nParamAcceptances   = C(end,4:5);
    nParamProposals     = C(end,6:7);
    logPriorGPParams = [log(sqrt(0.9*varData)),0.01; log(1/2), logPriorSDLength]; 

else
    
    % Read in data (genenames, featurenames, and the data)
    N = nGenes*nFeatures;
    if permuteData == true && isnan(permutationFileName)
        xx = 1:nFeatures;
        yy = [];
        uniqueCaptureTimes = unique(captureTimes);
        lu = length(uniqueCaptureTimes);
        for i = 1:lu
            xy = xx(captureTimes == uniqueCaptureTimes(i));
            yy = [yy randsample(xy,length(xy))];
        end
        data = data(:,yy);
    elseif permuteData == true && ~isnan(permutationFileName)
        perm = csvread(permutationFileName);
        data = data(:,perm);
        yy = perm;
    else
        yy = 1:nFeatures;
    end
    cellFileName = [strtok(fileName, '.'),'_Permutation', num2str(uniqueIdentifier)];
    CellFileName = [cellFileName '.csv'];
    csvwrite(CellFileName,yy); 
    permData = data;
    if regularInt == false
        tau1 = sqrt(sum((diff(permData,1,2).^2),1));
        pseudoTimes = cumsum([0 tau1])/sum(tau1);
    else
       pseudoTimes = (0.5:1:(nFeatures-0.5))/nFeatures;
    end
    currentOrder = 1:nFeatures;
    fid = fopen(outFileOrder, 'wt');
    fclose(fid);
    fid = fopen(outFile, 'wt');
    fclose(fid);
    logPriorGPParams = [log(sqrt(0.9*varData)),0.01; log(1/2), logPriorSDLength]; 
    currentLogParams = proposeParams(logPriorGPParams(:,1),logPriorGPParams(:,2)',varData);
    if currentLogParams(1) > log(sqrt(varData))
        currentLogParams(1) = log(sqrt(varData) - 0.05);
    end
end
currentGPPrior   = [log(normpdf(currentLogParams(1),logPriorGPParams(1,1), logPriorGPParams(1,2))), ...
                log(normpdf(currentLogParams(2),logPriorGPParams(2,1), logPriorGPParams(2,2)))];
    
[X, Y] = meshgrid(pseudoTimes);
timeDiffs               = (-(X - Y).^2);
%initial log-likelihood
currentSquaredParams = exp(2*currentLogParams);
K           = currentSquaredParams(1)*exp(timeDiffs/(2*currentSquaredParams(2))) ...
    + (varData-currentSquaredParams(1))*eye(nFeatures);

%compute the inverse and the determinant
    KChol       = chol(K);
    opts.LT     = false;
    opts.UT     = true;
    X           = linsolve(KChol,eye(nFeatures),opts);
    KInv        = X*X';
    logDetK     =  2*(sum(log(diag(KChol))));
 abc = 0;
    for i = 1:nGenes
        abc = abc + permData(i,:)*KInv*permData(i,:)';
    end

 currentLogMarginalLikelihood = -0.5*(N*log(2*pi)+logDetK+abc);

%compute the probabilities of choosing moves
distanceMatrix              = computeDistMat(data);
moveProbs = zeros(nFeatures^2,3);
moveProb1                = exp(-delta(1)* distanceMatrix.^2);
moveProb1                = moveProb1./sum(moveProb1);
moveProb1                = moveProb1.^jj;
moveProb1                = moveProb1./sum(moveProb1);
moveProbs(:,1)           = cumsum(moveProb1);
moveProb2                = exp(-delta(2)* distanceMatrix.^2);
moveProb2                = moveProb2./sum(moveProb2);
moveProb2                = moveProb2.^kk;
moveProb2                = moveProb2./sum(moveProb2);
moveProbs(:,2)           = cumsum(moveProb2 );
PP                       = cumsum(pp);



    
%     
    
    


nMcmc = nSamples;


 for sampleNumber = 1:nMcmc  
     
    %%% save down the current sample (only save every 'thinningFreq'-th 
    %%% sample)
    if(mod(sampleNumber,thinningFreq) == 0)
        dlmwrite(outFileOrder,currentOrder, '-append', 'delimiter',',');
        dlmwrite(outFile,[currentLogMarginalLikelihood currentLogParams' ...
            nParamAcceptances nParamProposals nOrderAcceptances nOrderProposals], '-append', 'delimiter',',');
    end

    if(verbose && mod(sampleNumber,1000) == 0)
            disp(['Sample number: ',num2str(sampleNumber)]);
            oAR = nOrderAcceptances./nOrderProposals;
            oAR = oAR(~isnan(oAR));
            disp(['Order acceptance rate = ', num2str(oAR)]);
            disp(['Parameter acceptance rate = ', num2str(nParamAcceptances./nParamProposals)]);
            disp(['Parameters: ', num2str([currentSquaredParams(1),sqrt(currentSquaredParams(2))])]);
    end    
    
% %sample the parameters
    if mod(sampleNumber,paramSamplingFreq)==0
        nParamProposals = nParamProposals + ones(1,2);
        proposedLogParams = proposeParams(currentLogParams,stepSize,varData);
        if regularInt == false
            [newLogLik newLogPriorOfLogHypers newLogHypers newSquaredHypers accepted] = ...
            updateHypersPseudoRank(proposedLogParams,currentSquaredParams, timeDiffs,permData,currentLogMarginalLikelihood,...
            currentGPPrior,logPriorGPParams, nFeatures, nGenes,varData);
        else
             [newLogLik newLogPriorOfLogHypers newLogHypers newSquaredHypers accepted,newKInv,newLogDet] = ...
            updateHypersPseudoRank(proposedLogParams,currentSquaredParams, timeDiffs,permData,currentLogMarginalLikelihood,...
            currentGPPrior,logPriorGPParams, nFeatures, nGenes,varData);
            if ~isnan(newLogDet)
                KInv = newKInv;
                logDetK = newLogDet;
            end
        end
        currentLogMarginalLikelihood    = newLogLik;
        currentGPPrior                  = newLogPriorOfLogHypers;
        currentLogParams        = newLogHypers;
        currentSquaredParams    = newSquaredHypers;
        nParamAcceptances       = nParamAcceptances + accepted;
        accRates                = nParamAcceptances./nParamProposals;
        if  sampleNumber <= burnInParams
            for k = 1:2
                if accRates(k) < 0.45
                    stepSize(k) = stepSize(k)*0.85;
                elseif accRates(k) > 0.5
                    stepSize(k) = stepSize(k)/0.85;
                end
            end
       end
        
    end
    %sampling the orders
    if regularInt == false
        [ppath tau X dataProposedOrder] = ProposingPseudoRanks(currentOrder,PP,data,moveProbs,nFeatures,n0,n3,n3a);   
        [XX, YY] = meshgrid(tau);
        proposedTimeDiffs   = (-(XX - YY).^2);
        [newLogLik  accepted] = updateOrdersPseudoRank([currentLogParams(1)  currentLogParams(2)....
            log(sqrt(varData-currentSquaredParams(1)))], proposedTimeDiffs,dataProposedOrder,currentLogMarginalLikelihood, nFeatures, nGenes);
    else
     [ppath X dataProposedOrder] = ProposingPseudoRanksRegular(currentOrder,PP,data,moveProbs,nFeatures,n0,n3,n3a);
     [newLogLik  accepted] = updateOrdersPseudoRankRegular(KInv, logDetK,dataProposedOrder,currentLogMarginalLikelihood, nFeatures, nGenes);    
    end
    nOrderProposals(X+1) = nOrderProposals(X+1) + 1;
    if accepted == 1
        nOrderAcceptances(X+1) = nOrderAcceptances(X+1) + 1;
        currentLogMarginalLikelihood = newLogLik;
        currentOrder = ppath;
        permData = dataProposedOrder;
        if regularInt == false
            timeDiffs = proposedTimeDiffs;
            pseudoTimes = tau;
        end
    end
    
 end   
end
