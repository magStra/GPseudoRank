function posFreq = computePositionsPseudoRank(fileName,uniqueIdentifiers,...
permFileIdentifiers,burnIn,captureTimes,saveFile)
%input: 
%fileName: name of data set, e.g. 'Shalek.csv'
%uniqueIdentifiers: unique identifiers of the chains 
%permFileIdentifiers: unique identifiers of the initial permutations of the
%orders, normally uniqueIdentifiers = permFileIdentifiers
%burnIn: number of iterations to be discarded
%captureTimes: capture times of the individual single cells
%saveFile: name of mat-file saving the cell positions for each chain and their
%means and standard deviations
%output: cell positions, additionally saving their means and standard
%deviations in a mat-file
lu = length(uniqueIdentifiers);
saveFileName    = [strtok(fileName, '.'),'_Results_Orders_Chain', num2str(uniqueIdentifiers(1))];
orders          = dlmread([saveFileName '.csv'],',', burnIn-1,0);
nCells          = size(orders,2);
posFreq         = zeros([nCells nCells lu]);

for k = 1:lu
    if k > 1
        saveFileName       = [strtok(fileName, '.'),'_Results_Orders_Chain', num2str(uniqueIdentifiers(k))];
        orders             = dlmread([saveFileName '.csv'],',', burnIn-1,0);
    end
    if size(captureTimes,1) == 1
        captureTimes = captureTimes';
    end
    permFileName       = [strtok(fileName, '.'),'_Permutation', num2str(permFileIdentifiers(k))];
    permut             = csvread([permFileName '.csv']);
    orders1 = orders;
    for j = 1:size(orders,1)
        orders(j,:)  = permut(orders1(j,:));
    end
    cors = corr(captureTimes,orders');
    orders(cors < 0,:) = orders(cors < 0,end:-1:1);
    posFreq(:,:,k) = convertToMatrix(orders');   
end
xx = [];%means of cell positions
yy = [];%sds of cell positions
for k = 1 : lu      
    M = posFreq(:,:,k);
    A = repmat((1:nCells),nCells,1);
    aa = sum(M.*A,2);
    xx = [xx aa];
    bb = sqrt(sum(M.*A.^2,2) - aa.^2);
    yy = [yy bb];
end

save(saveFile,'xx','yy','posFreq');
end
