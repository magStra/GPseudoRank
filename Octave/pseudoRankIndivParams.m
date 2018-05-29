function [] = pseudoRankIndivParams(fileName, uniqueIdentifier, nSamples, logPriorGPMeans, logPriorGPSds,...
    verbose, initialise, thinningFreq,paramSamplingFreq, stepSize, inputSeed, delta, pp,permuteData,captureTimes,...
    regularInt,permutationFileName,adjustAccRate, adjustParams,burnInParams,adjustForCellSize,saveParams)
    %This is the function to compare the efficiency of the different moves
    %for the proposal of the orders
    %Unlike pseudoRank.m this function uses individual parameters for each
    %gene, not sharing them across different genes. 
initSampleNumber = 1;
if size(stepSize,2) == 3
    stepSize = stepSize';
end
if(isnan(inputSeed))
    inputSeed = sum(100*clock) + sum(100*double(uniqueIdentifier)); 
end

rng(inputSeed);

saveFileNameOrder = [strtok(fileName, '.'),'_Results_Orders_Chain', num2str(uniqueIdentifier)];
saveFileName      = [strtok(fileName, '.'),'_Results_Chain', num2str(uniqueIdentifier)];
saveFileStepSize  = [strtok(fileName, '.'),'_StepSize_Chain', num2str(uniqueIdentifier)];
saveFileNameN  = [strtok(fileName, '.'),'_N_Chain', num2str(uniqueIdentifier)];

allData = importdata(fileName, ',');
try
    data    = allData.data;
catch
    data    = allData;%if there are no gene names and cell names
end
[nGenes nFeatures] = size(data);
outFileOrder = [saveFileNameOrder '.csv'];
outFile    = [saveFileName '.csv'];
outFileStepSize = [saveFileStepSize,'.csv'];
outFileN = [saveFileNameN '.csv'];
nOrderProposals     = [0 0 0 0 0];
nOrderAcceptances   = [0 0 0 0 0];
nParamProposals     = zeros(nGenes,3);
nParamAcceptances   = zeros(nGenes,3);
if(~initialise)
    % Set up from the csv files
    B = csvread([saveFileNameOrder,'.csv'],0,0);
    currentOrder = B(end,:);
    permData = data(:,currentOrder);
    if regularInt == false
        tau1 = sqrt(sum((diff(permData,1,2).^2),1));
        pseudoTimes = cumsum([0 tau1])/sum(tau1);
    else
       pseudoTimes = (0.5:1:(nFeatures-0.5))/nFeatures;
    end
    C = csvread([saveFileName,'.csv']);
    
    nOrderAcceptances = C(end,2:6);
    nOrderProposals = C(end,7:11);
    if saveParams == true
        currentLogParams = C(end,12:end);
        currentLogParams = reshape(currentLogParams,3,nGenes);
    end
    stepSize = csvread([saveFileStepSize,'.csv']);
    NN = csvread([saveFileNameN,'.csv']);
    n0 = NN(1);
    n3 = NN(4);
    n3a = NN(5);
    kk = NN(3);
    jj = NN(2);
    initSampleNumber = size(B,1)+1;
    
else
    if (permuteData == true) && any(isnan(permutationFileName))
        xx = 1:nFeatures;
        yy = [];
        uniqueCaptureTimes = unique(captureTimes);
        lu = length(uniqueCaptureTimes);
        for i = 1:lu
            xy = xx(captureTimes == uniqueCaptureTimes(i));
            yy = [yy randsample(xy,length(xy))];
        end
        data = data(:,yy);
        
        cellFileName = [strtok(fileName, '.'),'_Permutation', num2str(uniqueIdentifier)];
        CellFileName = [cellFileName '.csv'];
        csvwrite(CellFileName,yy); 
    elseif (permuteData == true) && (~any(isnan(permutationFileName)))
        perm = csvread(permutationFileName);
        data = data(:,perm);
    end
    permData = data;
    if regularInt == false
        tau1 = sqrt(sum((diff(permData,1,2).^2),1));
        pseudoTimes = cumsum([0 tau1])/sum(tau1);
    else
       pseudoTimes = (0.5:1:(nFeatures-0.5))/nFeatures;
    end
    if isnan(adjustParams)
        n0 = max(floor(nFeatures/4),1);
        n3 = max(1,floor(nFeatures/20));
        n3a = max(3,floor(nFeatures/12));
        kk = 1;
        jj = 1;
    else
        n0 = adjustParams(1);
        n3 = adjustParams(4);
        n3a = adjustParams(5);
        kk = adjustParams(3);
        jj = adjustParams(2);
    end
    currentOrder = 1:nFeatures;
    fid = fopen(outFileOrder, 'wt');
    fclose(fid);
    fid = fopen(outFile, 'wt');
    fclose(fid);
    currentLogParams = proposeIndivParams(logPriorGPMeans,logPriorGPSds,nGenes);
end
%adjust permData for cell size
if adjustForCellSize == true
    uniqueCaptureTimes = unique(captureTimes);
    lu = length(uniqueCaptureTimes);
    permData1 = [];
    for j = 1:lu
        aa = permData(:,captureTimes==uniqueCaptureTimes(j));
        aaMod = bsxfun(@minus, aa,mean(aa,2));
        aa1 = bsxfun(@minus,aa,median(aaMod,1));
        permData1 = [permData1,aa1];
    end
    permData = permData1;
end
currentGPPrior = zeros(nGenes,3);
for j  = 1:nGenes
    currentGPPrior(j,:)  = [log(normpdf(currentLogParams(1,j),logPriorGPMeans(1,j), logPriorGPSds(1,j))), ...
                    log(normpdf(currentLogParams(2,j),logPriorGPMeans(2,j), logPriorGPSds(2,j))),...
                    log(normpdf(currentLogParams(3,j),logPriorGPMeans(3,j), logPriorGPSds(3,j)))];
end
[X, Y] = meshgrid(pseudoTimes);
timeDiffs               = abs(X - Y);

%initial log-likelihood
currentL    = exp(currentLogParams(2,:));
currentW2   = exp(2*currentLogParams(1,:));
currentE2   = exp(2*currentLogParams(3,:));
K = zeros(nFeatures,nFeatures,nGenes);
KInv = zeros(nFeatures,nFeatures,nGenes);
logDetK         = zeros(1,nGenes);
currentLogLks   = zeros(1,nGenes);
opts.UT = true;

for j = 1:nGenes
 K(:,:,j)           = currentW2(j) *exp(-timeDiffs.^2/(2*currentL(j)^2)) ...
         + currentE2(j) * eye(nFeatures);
%compute the inverse and the determinant

KChol       = chol(K(:,:,j));
X           = linsolve(KChol,eye(nFeatures),opts);
KInv(:,:,j) = X*X';
logDetK(j)  =  2*(sum(log(diag(KChol))));

abc =  permData(j,:)*KInv(:,:,j)*permData(j,:)';
currentLogLks(j) = -0.5*(nFeatures*log(2*pi)+logDetK(j)+abc);
end
currentLogLk = sum(currentLogLks);

%compute the probabilities of choosing moves
distanceMatrix1             = bsxfun(@minus,data,permute(data,[1 3 2]));
distanceMatrix              = sum(reshape(abs(distanceMatrix1),size(data,1),[]));
moveProbs = zeros(nFeatures^2,2);
moveProb1                = exp(-delta(1)* distanceMatrix.^2);
moveProb1                = moveProb1./sum(moveProb1);
moveProb1                = moveProb1.^jj;
moveProb1                = moveProb1./sum(moveProb1);
moveProbs(:,1)           = cumsum(moveProb1);
moveProb2                = exp(-delta(2)* distanceMatrix.^2);
moveProb2      = moveProb2./sum(moveProb2);
moveProb2 = moveProb2.^kk;
moveProb2      = moveProb2./sum(moveProb2);
moveProbs(:,2)   = cumsum(moveProb2 );
PP            = cumsum(pp);

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
nMcmc = nSamples;


 for sampleNumber = initSampleNumber:nMcmc  
    %%% save down the current sample (only save every 'thinningFreq'-th 
    %%% sample)
    if(mod(sampleNumber,thinningFreq) == 0)
        dlmwrite(outFileOrder,currentOrder, '-append', 'delimiter',','); 
        if saveParams == true
        dlmwrite(outFile,[currentLogLk  nOrderAcceptances ...
            nOrderProposals currentLogParams(:)'], '-append', 'delimiter',',');
        else
            dlmwrite(outFile,[currentLogLk  nOrderAcceptances ...
            nOrderProposals], '-append', 'delimiter',',');
        end
    end
    orderAccRates = nOrderAcceptances./nOrderProposals;
    if(verbose && mod(sampleNumber,10) == 0)
            disp(['Sample number: ',num2str(sampleNumber)]);
            disp(['Order acceptance rate = ', num2str(orderAccRates)]);
            if saveParams == true
                disp(['log-Parameters: ', num2str(currentLogParams(:)')]);
                disp(['Parameter acceptance rate = ', num2str(mean(nParamAcceptances./nParamProposals))]);
            end
    end    
% %sample the parameters
    if mod(sampleNumber,paramSamplingFreq)==0
        nParamProposals = nParamProposals + ones(nGenes,3);
        proposedLogParams = proposeIndivParams(currentLogParams,stepSize,nGenes);
        currentAccepted = zeros(nGenes,3);
        if regularInt == false
            for j = 1:nGenes
            [newLogLik newGPPrior newLogParams  acc] = ...
            updateHypersPseudoRankIndivParams(proposedLogParams(:,j), currentLogParams(:,j), timeDiffs,permData(j,:),currentLogLks(j),...
            currentGPPrior(j,:),logPriorGPMeans(:,j),logPriorGPSds(:,j), nFeatures);
             currentLogParams(:,j) = newLogParams;
             currentLogLks(j) = newLogLik;
             currentGPPrior(j,:) = newGPPrior;
             currentAccepted(j,:) = acc';
            end   
        else
            for j = 1:nGenes
            [newLogLik newGPPrior newLogParams  acc newKInv newLogDetK] = ...
            updateHypersPseudoRankIndivParamsRegular(proposedLogParams(:,j), currentLogParams(:,j), timeDiffs,permData(j,:),currentLogLks(j),...
            currentGPPrior(j,:),logPriorGPMeans(:,j),logPriorGPSds(:,j), nFeatures);
             currentLogParams(:,j) = newLogParams;
             currentLogLks(j) = newLogLik;
             currentGPPrior(j,:) = newGPPrior;
             currentAccepted(j,:) = acc';
             if ~isnan(newLogDetK)
                 KInv(:,:,j) = newKInv;
                 logDetK(j) = newLogDetK;
             end
            end
        end

        currentLogLk = sum(currentLogLks);
        nParamAcceptances = nParamAcceptances + currentAccepted;
        accRates = nParamAcceptances./nParamProposals;
        if  sampleNumber <= burnInParams
        for j = 1:nGenes
            for k = 1:3
                if accRates(j,k) < 0.45
                    stepSize(k,j) = stepSize(k,j)*0.95;
                elseif accRates(j,k) > 0.5
                    stepSize(k,j) = stepSize(k,j)*1.05;
                end
           end
        end
       end
        
        
    end
       
    %sampling the orders
    if regularInt == false
        [ppath tau X dataProposedOrder] = ProposingPseudoRanks(currentOrder,PP,data,moveProbs,nFeatures,n0,n3,n3a);   
        [XX, YY] = meshgrid(tau);
        proposedTimeDiffs   = abs(XX - YY);
        [newLogLks newLogLik  accepted] = updateOrdersPseudoRankIndivParams(currentLogParams,...
            proposedTimeDiffs,dataProposedOrder,currentLogLk, nFeatures, nGenes);
    else
     [ppath X dataProposedOrder] = ProposingPseudoRanksRegular(currentOrder,PP,data,moveProbs,nFeatures,n0,n3,n3a);
     
     [newLogLks newLogLik  accepted] = updateOrdersPseudoRankRegularIndivParams(KInv, logDetK,dataProposedOrder, currentLogLk,nFeatures, nGenes);    
   
    end   
    nOrderProposals(X+1) = nOrderProposals(X+1) + 1;
    if accepted == 1
        nOrderAcceptances(X+1) = nOrderAcceptances(X+1) + 1;
        currentLogLk = newLogLik;
        try
         currentLogLks = newLogLks;
        catch
        end
        currentOrder = ppath;
        permData = dataProposedOrder;
        if regularInt == false
            timeDiffs = proposedTimeDiffs;
        end
    end
        if sampleNumber <= burnInParams && adjustAccRate == true && mod(sampleNumber,50) == 0 
        if orderAccRates(1) < 0.2
            n0 = max(n0 - 1,floor(nFeatures/30));
        end
        if orderAccRates(1) > 0.3
            n0 = max(n0 - 1,floor(nFeatures/3));
        end
        if orderAccRates(4) < 0.1
            if rand < 0.5
                n3 = max(n3-1,floor(nFeatures/50));
            else
                n3a = max(n3a-1,4);
            end
        end
        if orderAccRates(4) > 0.3
            if rand < 0.5
                n3 = max(n3+1,floor(nFeatures/10));
            else
                n3a = max(n3a+1,floor(nFeatures/10));
            end
        end
        if orderAccRates(2) > 0.3
            jj = jj*0.98;
            moveProb1 = moveProb1.^jj;
            moveProb1 = moveProb1/sum(moveProb1);
            moveProbs(:,1) = cumsum(moveProb1);
        end
        if orderAccRates(2) < 0.1
            jj = jj*(1/0.98);
           moveProb1 = moveProb1.^jj;
            moveProb1 = moveProb1/sum(moveProb1);
            moveProbs(:,1) = cumsum(moveProb1);
        end
       
        if orderAccRates(3) < 0.1
            kk = kk *(1/0.98);
            moveProb2 = moveProb2.^(1/0.99);
            moveProb2 = moveProb2/sum(moveProb2);
            moveProbs(:,2) = cumsum(moveProb2);
        end
        if orderAccRates(3) > 0.3
            kk = kk*0.98;
            moveProb2 = moveProb2.^0.99;
            moveProb2 = moveProb2/sum(moveProb2);
            moveProbs(:,2) = cumsum(moveProb2);
        end
        end
        
    if sampleNumber == burnInParams
        csvwrite(outFileStepSize,stepSize);
        csvwrite(outFileN,[n0 jj kk n3 n3a]);
    end
 end 
end
