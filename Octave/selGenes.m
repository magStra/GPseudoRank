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
    p(j) = anova(data(j,:),captureTimes);
    end
    [x inds] = sort(p);
    nGenes = size(data,1);
    dataRed = data(inds(1:min(floor(nGenes*0.8),100)),:);
    genesRed = inds(inds(1:min(floor(nGenes*0.8),100)));
else %if there are no capture times
    inds3 = [];
    [nGenes,nCells] = size(data);
    for jk = 1:nGenes
        if sum(data(jk,:)~=0)>nCells*0.7
            inds3 = [inds3;jk];%reduce if needed
        end
    end
    xx = 1:nGenes;
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

