function [simData,tau] = simGPIndiv(T,sw2,l,se2,simDataFileName,tauFileName,regular,distri)
%Simulates GP data with T timepoints and given vectors of GP parameter
%values. The pseudotime points may be chosen to be regular, or randomly
%sampled from a uniform distribution or several normal distributions around
%capture times. 
%input:
%T: number of cells to simulated
%sw2: vector of values of scale parameter for the GP of each gene, vector 
%of length n, where n is the number of genes to be simulated
%l: vector of values of length parameter, of length n
%se2: vector of values of noise parameter, of length n
%simDataFileName: name of the .csv file where you save the simulated data
%tauFileName: name of output .csv file containing the simulated pseudotime
%points
%regular: if true, we do not use pseudotimes sampled from a distribution,
%but equidistant ones
%distri:  is either 'uniform' or 'normal', set to any other value if you
%use "regular=false"
%for "uniform", we draw T samples from U(0,1) and sort them
%for "normal", we use 4 capture times 0.125,0.375,0.625,0.875 
%the capture times are the means and the sds are 0.125
%output files: 
%The nxT matrix of simulated data is saved as simDataFileName.
%The simulated pseudotime points are saved as tauFileName.
if regular == false
    if distri == 'uniform'
        tau = sort(rand([1 T]));
    elseif distri == 'normal'
    tau = randn([1 T])*0.125 + [repmat(0.125,[1 floor(T/4)]),repmat(0.375,[1 floor(T/4)]),...
        repmat(0.625,[1 floor(T/4)]),repmat(0.875,[1 T-3*floor(T/4)])];
    end
else%regular pseudotime points
    tau = (0.5:(T-0.5))/T;
end
[X,Y] = meshgrid(tau,tau);
lw = length(sw2);
K = zeros(T,T,lw);
for j = 1:lw
    K(:,:,j) = sw2(j) * exp(-(X-Y).^2/l(j)) + se2(j)*eye(T);
end
simData = zeros(lw,T);
for j = 1:lw
    simData(j,:) = mvnrnd(zeros([1 T]),K(:,:,j),1);
end
csvwrite(simDataFileName,simData);
csvwrite(tauFileName,tau);
end
