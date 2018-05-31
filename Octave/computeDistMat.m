function [distMat] = computeDistMat(data)
%computes the matrix of L1 distances required for 
% the computation of the probabilites of reversal moves
T = size(data,2);
distMat = zeros(T,T);
for j = 2:T
    for i = 1:(j-1)
        distMat(i,j) = sum(abs(data(:,i)-data(:,j)));
    end
end
distMat = distMat+distMat';
distMat = distMat(:)';
end

