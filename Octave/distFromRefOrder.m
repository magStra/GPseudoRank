function [distVecL1] = distFromRefOrder(orders,refOrder)
%orders: each row is one permutation
%ref order is a row vector
[n T] = size(orders);
distVecL1 = zeros(n,1);
distVecSpearman = zeros(n,1);
for j = 1:n
    A = orders(j,:);
    pos1(A) = 1:T;
   distVecL1(j) = pdist([pos1;refOrder],'cityblock');
   %distVecSpearman(j) = pdist([pos1;refOrder],'spearman'); slow!
end
end

