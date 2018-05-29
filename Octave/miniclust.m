function [] = miniclust(fileName,CT)
%mini-cluster approximation for large data sets
%file name: name of .csv file containing expression matrix, rows are genes,
%columns are cells
%CT:vector of capture times (cells all might have been capture at the
%same time)
allData = importdata(fileName, ',');
try
    A    = allData.data;
catch
    A    = allData;%if there are no gene names and cell names
end
uCT = unique(CT);
b = max(3,ceil(length(CT)/200));
c = min(floor(length(CT)/(3*length(uCT))),ceil(200/length(uCT)));
captureTimes = [];
subS = [];
assignments = {};
for j = 1:(length(uCT))
    aa = CT==uCT(j);    
    a = sum(aa);
    [B BB] = kmeans(A(:,CT==uCT(j))',max(c,floor(a/b)),'Distance','cityblock');
    subS = [subS,BB'];
    captureTimes = [captureTimes,repmat(uCT(j),1,max(c,floor(a/b)))];
    assignments{j} = B;
end

csvwrite([strtok(fileName, '.'),'KM.csv'],subS);
csvwrite([strtok(fileName, '.'),'KMCT.csv'],captureTimes);
%save also the allocations to the mini-clusters
ass = assignments{1};
for j = 2:length(assignments)
    a = max(ass);
    aa = a + assignments{j};
    ass = [ass;aa];
end
csvwrite([strtok(fileName, '.'),'KMAss.csv'],ass);
end

