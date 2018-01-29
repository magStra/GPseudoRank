function [newLogLks newLogLik  accepted] = ...
    updateOrdersPseudoRankRegularIndivParams(KInv, logDetK,dataProposedOrder,...
    logLik, nTimes, nGenes)
accepted = 0;
newLogLik = logLik;
newLogLks = zeros(1,nGenes);
for j = 1:nGenes
    newLogLks(j) = -0.5*(nTimes*log(2*pi)+logDetK(j)+dataProposedOrder(j,:)*KInv(:,:,j)*dataProposedOrder(j,:)');
end 
proposedLogMarginalLikelihood = sum(newLogLks);
if log(rand) < proposedLogMarginalLikelihood-logLik
         newLogLik = proposedLogMarginalLikelihood;
         accepted = 1;
end
end
