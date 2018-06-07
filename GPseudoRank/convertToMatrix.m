function meanPermutMatrix = convertToMatrix(permuts)
%converts permutations to permutation matrices
%and returns the mean of the permutation matrices
%input:
%each column in permuts is a permutation
%output: 
%returns mean of the permutations
[lp np] = size(permuts);
I(permuts(:,1)) = 1:lp;
meanPermutMatrix = zeros(lp,lp);
permutMatrix = zeros(lp,lp);
meanPermutMatrix(:,I) = eye(lp);

for j = 2:np
    I(permuts(:,j)) = 1:lp;
    permutMatrix(:,I) = eye(lp);
    meanPermutMatrix = (permutMatrix + meanPermutMatrix*(j-1))/j;
end
end

