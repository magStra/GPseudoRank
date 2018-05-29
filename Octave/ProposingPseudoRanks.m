function [ppath tau X dataProposedOrder] = ProposingPseudoRanks(order,PP,data,moveProbs,nFeatures,...
    n0,n3,n3a)
%This function samples from the proposal distribution of the orders
%and returns the corresponding order and  pseudotimes

X = sum(PP<=rand);

switch X
    
    
    case 0
        % swaps of neighbours
        ppath       = order;
        aaa         = max(n0,1);
        bbb         = randi(aaa);
        ii1         = randsample(nFeatures-1,bbb);
        for j = 1:bbb
            i1 = ii1(j);
            i2 = i1+1;
            helpVar     = ppath([i1,i2]);
            ppath(i1)   = helpVar(2);
            ppath(i2)   = helpVar(1);
        end
   
    case 1
         %second case, swaps of similar elements
         Y         = sum(moveProbs(:,1) <= rand)+1;
         ppath     = order;
          m        = mod(Y,nFeatures);
         n         = ceil(Y/nFeatures);
         if m == 0
             m = nFeatures;
         end
         i1         = find(ppath == n);
         i2         = find(ppath ==  m);
         ppath(i1)  = order(i2);
         ppath(i2)  = order(i1);
      case 2% third case, reversals
         Y         = sum(moveProbs(:,2) <= rand)+1;
         ppath      = order;
         m          = mod(Y,nFeatures);
         n          = ceil(Y/nFeatures);
         if m == 0
             m = nFeatures;
             n = n-1;
         end
         i1         = find(ppath == n);
         i2         = find(ppath ==  m);
         ppath(min(i1,i2):max(i1,i2))  = order(max(i1,i2):-1:min(i1,i2)); 
     
         
    
      case 3%random permutations
        %draw the number of such permutations
        nPerm     = randi(max(1,n3));
        ppath     = order;
        %for each such permutation draw the length
        % we use a length of at least 3
        lPerms  = randi([3 max(n3a,3)],[1 nPerm]);
        %lPerms  = randi([3 8],[1 nPerm]);
        for i = 1:nPerm
        %now for each permutation pick the first index
            i1 = randi(nFeatures-lPerms(i)+1);
            i2 = i1 + lPerms(i)-1;
            xx = ppath(i1:i2);
            helpVar = xx(randperm(i2-i1+1));
            ppath(i2:-1:i1) = helpVar;
        end 
      case 4%reversal of entire path
          ppath = order(end:-1:1);
          
end 
%now convert path to pseudo-time
dataProposedOrder = data(:,ppath);
tau1 = sqrt(sum((diff(dataProposedOrder,1,2).^2),1));
tau = cumsum([0 tau1])/sum(tau1);

if any(unique(order) ~= unique(ppath))
    'error -- new order is not a permutation of old one'
end 
if any(isnan(ppath))
   'error: NaN in order'
end
end

