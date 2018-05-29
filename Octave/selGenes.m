function [] = selGenes(fileName,varargin)
%input: name of .csv file with RNA-expression matrix
%optional argument: capture times
%the function writes a .csv file with the selected rows of the
%RNA-expression matrix

allData = importdata(fileName, ',');
try
    data    = allData.data;
catch
    data    = allData;%if there are no gene names and cell names
end
if nargin == 2
    captureTimes = varargin{1};
    for j = 1:size(data,1)
    p(j) = anova1(data(j,:),captureTimes,'off');
    end
    [x inds] = sort(p);
    nGenes = size(data,1);
    dataRed = data(inds(1:min(floor(nGenes*0.8),100)),:);
    genesRed = inds(inds(1:min(floor(nGenes*0.8),100)));
else %if there are no capture times
    inds3 = [];
    for jk = 1:size(data,1)
        inds3 = [inds3;sum(data(jk,:)~=0)>101*0.7];
    end
    xx = 1:size(data,1);
    all1 = data(inds3,:);
    xx = xx(inds3);
    nGenes = size(all1,1);
    [x inds] = sort(var(all1,[],2),'descend');
    inds = inds(1:min(floor(nGenes*0.8),500));
    [y inds1] = sort(mean(all1,2),'descend');
    inds1 = inds1(1:min(floor(nGenes*0.8),500));
    inds2 = intersect(inds1,inds);
    dataRed = all1(inds2,:);
    genesRed = xx(inds2);
end
csvwrite([strtok(fileName, '.'),'Red.csv'],dataRed);
csvwrite([strtok(fileName, '.'),'GenesRed.csv'],genesRed);
end

