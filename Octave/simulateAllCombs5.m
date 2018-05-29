function [] = simulateAllCombs5(identifier)
ppar = csvread(sprintf('sim%d_paramsMeansBasic.csv',identifier));
fHandle          = @pseudoRankIndivParams;  
fileName         = sprintf('sim%dBasic.csv',identifier);
nSamples         = 100000; 
logPriorGPMeans = ppar;%the real parameters are the means
logPriorGPSds = zeros(3,50);    
paramSamplingFreq = 10^6; %fixed parameters
verbose          = false; 
initialise       = true;  
thinningFreq     = 10;     
inputSeed        = NaN;    
delta            = [1/1000,1/1000];
permuteData      = true;
captureTimes     = [zeros(1,30),ones(1,30),repmat(2,1,30)];
regInt           = true;
permutationFileName = NaN;
stepSize = ones(3,16);
burnInAdapt = 5000;

%move 1
pp               = [0.998 0 0 0 0.002];
tic
for k = 6:10
feval(fHandle, fileName, k+100, nSamples, logPriorGPMeans,logPriorGPSds, ...
    verbose, initialise, thinningFreq, paramSamplingFreq,stepSize, inputSeed, delta, pp,permuteData,captureTimes,...
    regInt,permutationFileName,true,NaN,burnInAdapt,false,false)
end
toc
%move 2
pp               = [0 0.998 0 0 0.002];
tic
for k = 6:10
    permutationFileName = sprintf('sim%dBasic_Permutation%d.csv',identifier,k+100);
    feval(fHandle, fileName, k+200, nSamples, logPriorGPMeans,logPriorGPSds, ...
    verbose, initialise, thinningFreq, paramSamplingFreq,stepSize, inputSeed, delta, pp,permuteData,captureTimes,...
    regInt,permutationFileName,true,NaN,burnInAdapt,false,false)
end
toc
%move3
pp               = [0 0 0.998 0 0.002];
tic
for k = 6:10
    permutationFileName = sprintf('sim%dBasic_Permutation%d.csv',identifier,k+100);
feval(fHandle, fileName, k+300, nSamples, logPriorGPMeans,logPriorGPSds, ...
    verbose, initialise, thinningFreq, paramSamplingFreq,stepSize, inputSeed, delta, pp,permuteData,captureTimes,...
    regInt,permutationFileName,true,NaN,burnInAdapt,false,false)
end
toc
%move 4
pp               = [0 0 0 0.998 0.002];
tic
for k = 6:10
    permutationFileName = sprintf('sim%dBasic_Permutation%d.csv',identifier,k+100);
feval(fHandle, fileName, k+400, nSamples, logPriorGPMeans,logPriorGPSds, ...
    verbose, initialise, thinningFreq, paramSamplingFreq,stepSize, inputSeed, delta, pp,permuteData,captureTimes,...
    regInt,permutationFileName,true,NaN,burnInAdapt,false,false)
end
toc

%moves 1 and 2
pp = [0.499 0.499 0 0 0.002];
tic
for k = 6:10
    permutationFileName =sprintf('sim%dBasic_Permutation%d.csv',identifier,k+100);
feval(fHandle, fileName, k+1200, nSamples, logPriorGPMeans,logPriorGPSds, ...
    verbose, initialise, thinningFreq, paramSamplingFreq,stepSize, inputSeed, delta, pp,permuteData,captureTimes,...
    regInt,permutationFileName,true,NaN,burnInAdapt,false,false)
end
toc
%moves 1 and 3
pp = [0.499 0 0.499 0 0.002];
tic
for k = 6:10
    permutationFileName = sprintf('sim%dBasic_Permutation%d.csv',identifier,k+100);
feval(fHandle, fileName, k+1300, nSamples, logPriorGPMeans,logPriorGPSds, ...
    verbose, initialise, thinningFreq, paramSamplingFreq,stepSize, inputSeed, delta, pp,permuteData,captureTimes,...
   regInt,permutationFileName,true,NaN,burnInAdapt,false,false)
end
toc
%moves 1 and 4
pp = [0.499 0 0 0.499 0.002];
tic
for k = 6:10
    permutationFileName = sprintf('sim%dBasic_Permutation%d.csv',identifier,k+100);
feval(fHandle, fileName, k+1400, nSamples, logPriorGPMeans,logPriorGPSds, ...
    verbose, initialise, thinningFreq, paramSamplingFreq,stepSize, inputSeed, delta, pp,permuteData,captureTimes,...
    regInt,permutationFileName,true,NaN,burnInAdapt,false,false)
end
toc
%moves 2 and 3
pp = [0 0.499 0.499 0 0.002];
tic
for k = 6:10
    permutationFileName = sprintf('sim%dBasic_Permutation%d.csv',identifier,k+100);
feval(fHandle, fileName, k+2300, nSamples, logPriorGPMeans,logPriorGPSds, ...
    verbose, initialise, thinningFreq, paramSamplingFreq,stepSize, inputSeed, delta, pp,permuteData,captureTimes,...
    regInt,permutationFileName,true,NaN,burnInAdapt,false,false)
end
toc
%moves 2 and 4
pp = [0 0.499 0 0.499 0.002];
tic
for k = 6:10
    permutationFileName = sprintf('sim%dBasic_Permutation%d.csv',identifier,k+100);
feval(fHandle, fileName, k+2400, nSamples, logPriorGPMeans,logPriorGPSds, ...
    verbose, initialise, thinningFreq, paramSamplingFreq,stepSize, inputSeed, delta, pp,permuteData,captureTimes,...
    regInt,permutationFileName,true,NaN,burnInAdapt,false,false)
end
toc
%moves 3 and 4
pp = [0 0 0.499 0.499 0.002];
tic
for k = 6:10
    permutationFileName = sprintf('sim%dBasic_Permutation%d.csv',identifier,k+100);
feval(fHandle, fileName, k+3400, nSamples, logPriorGPMeans,logPriorGPSds, ...
    verbose, initialise, thinningFreq, paramSamplingFreq,stepSize, inputSeed, delta, pp,permuteData,captureTimes,...
    regInt,permutationFileName,true,NaN,burnInAdapt,false,false)
end
toc
%moves 1,2,3
pp = [0.357 0.357 0.357 0 0.0019];
tic
for k = 6:10
    permutationFileName = sprintf('sim%dBasic_Permutation%d.csv',identifier,k+100);
feval(fHandle, fileName, k+12300, nSamples, logPriorGPMeans,logPriorGPSds, ...
    verbose, initialise, thinningFreq, paramSamplingFreq,stepSize, inputSeed, delta, pp,permuteData,captureTimes,...
   regInt,permutationFileName,true,NaN,burnInAdapt,false,false)
end
toc
%moves 1,2,4
pp = [0.357 0.357 0 0.357 0.0019];
tic
for k = 6:10
    permutationFileName = sprintf('sim%dBasic_Permutation%d.csv',identifier,k+100);
feval(fHandle, fileName, k+12400, nSamples, logPriorGPMeans,logPriorGPSds, ...
    verbose, initialise, thinningFreq, paramSamplingFreq,stepSize, inputSeed, delta, pp,permuteData,captureTimes,...
    regInt,permutationFileName,true,NaN,burnInAdapt,false,false)
end
toc
%moves 2,3,4
pp = [0 0.357 0.357 0.357 0.0019];
tic
for k = 6:10
    permutationFileName = sprintf('sim%dBasic_Permutation%d.csv',identifier,k+100);
feval(fHandle, fileName, k+23400, nSamples, logPriorGPMeans,logPriorGPSds, ...
    verbose, initialise, thinningFreq, paramSamplingFreq,stepSize, inputSeed, delta, pp,permuteData,captureTimes,...
    regInt,permutationFileName,true,NaN,burnInAdapt,false,false)
end
toc
%moves 1,3,4
pp = [0.357 0.357 0 0.357 0.0019];
tic
for k = 6:10
    permutationFileName = sprintf('sim%dBasic_Permutation%d.csv',identifier,k+100);
feval(fHandle, fileName, k+13400, nSamples, logPriorGPMeans,logPriorGPSds, ...
    verbose, initialise, thinningFreq, paramSamplingFreq,stepSize, inputSeed, delta, pp,permuteData,captureTimes,...
   regInt,permutationFileName,true,NaN,burnInAdapt,false,false)
end
toc
%all moves
pp = [0.2495 0.2495 0.2495 0.2495 0.002];
tic
for k = 6:10
    permutationFileName = sprintf('sim%dBasic_Permutation%d.csv',identifier,k+100);
feval(fHandle, fileName, k+123400, nSamples, logPriorGPMeans,logPriorGPSds, ...
    verbose, initialise, thinningFreq, paramSamplingFreq,stepSize, inputSeed, delta, pp,permuteData,captureTimes,...
    regInt,permutationFileName,true,NaN,burnInAdapt,false,false)
end
toc

end

