function [n0,n3,n3a,jj,kk,delta,pp] = setDefaultParams(fileName,captureTimes)
%given an input .csv file containing a matrix of gene expression levels 
%this function computes default parameters for the proposal distribution
%This default will work for most, but not all data sets.
%If necessary, the parameters should be adapted as follows:
%jj,kk and delta need to be increased in case of low acceptance rate.
%For large data sets, the function miniclust.m should be run before this
%one.
%captureTimes is a logical variable indicating whether there are several
%different capture times
allData = importdata(fileName, ',');
try
    data    = allData.data;
catch
    data    = allData;%if there are no gene names and cell names
end
[nG,nC] = size(data)
if nC > 700
    'Use the mini-cluster approximation (miniclust.m) before computing the default parameters'
    n0 = NaN;
    n3 = NaN;
    n3a = NaN;
    jj = NaN;
    kk = NaN;
    delta = NaN;
else
    if nC < 350 && captureTimes == true 
        pp = [repmat(0.998/4,1,4),0.002];
        n3 = floor(nC/20);
        n3a = floor(nC/12);
    elseif nC > 30 && nC < 350
        pp = [0,repmat(0.998/3,1,3),0.002];
        n3 = floor(nC/20);
        n3a = floor(nC/12);
    elseif nC > 350 
         pp = [0,0,repmat(0.998/2,1,2),0.002];
         n3 = floor(nC/20);
         n3a = floor(nC/12);
    else

    end
    n0 = floor(nC/7);
    delta = 1/8000;
    jj = 0.1;
    kk = 0.1;
    if nC < 50
        kk = 0.01;
    end
end
end

