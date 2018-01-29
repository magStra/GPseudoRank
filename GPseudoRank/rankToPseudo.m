function [pseudoTimes,positions] = rankToPseudo(fileNameOrders,fileName,initPermutFileName,burnIn,captureTimes)
%The function (re)computes the pseudotimes for given samples, given the sampled orders
%and the data.
%fileNameOrders: csv-file with matrix of cell orderings, each line is one sampled order
%fileName: name of csv-file containing the data matrix, genes along the
%rows and cells along the columns
%initPermutFileName: csv-file containing the starting permutation of the
%chain; this is an output file of the pseudoRank function
%burnIn: number of initial samples to be discarded
%captureTimes: capture times at which the cells were measured
%output: pseudotimes of the cells, positions of the cells in the order
permut  = csvread(initPermutFileName);
allData = importdata(fileName, ',');
T = length(permut);
orders         = csvread(fileNameOrders);
m = size(orders,1)-burnIn;
pseudoTimesA    = zeros(m,T);
pseudoTimes     = zeros(m,T);
positions       = zeros(m,T);
try
    data    = allData.data;
catch
    data    = allData;%if there are no gene names and cell names
end
data = data(:,permut);
if size(captureTimes,1) == 1
        captureTimes = captureTimes';
end
cors = corr(captureTimes,orders');
orders(cors < 0,:) = orders(cors < 0,end:-1:1);
for j = 1:m
    permData    = data(:,orders(j,:));
    tau1        = sqrt(sum((diff(permData,1,2).^2),1));
    pseudoTimesA(j,:) = cumsum([0 tau1])/sum(tau1); 
    invPerm(permut(orders(j,:))) = 1:T;
    positions(j,:) = invPerm;
    pseudoTimes(j,:)  = pseudoTimesA(j,invPerm);
end
end