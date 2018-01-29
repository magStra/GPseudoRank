function [newLogLik  accepted] = ...
    updateOrdersPseudoRankRegular(KInv, logDetK,dataProposedOrder,logLik, nTimes, nGenes)
N = nGenes*nTimes;
accepted = 0;
newLogLik = logLik;

 abc = 0;
for i = 1:nGenes
    abc = abc + dataProposedOrder(i,:)*KInv*dataProposedOrder(i,:)';
end 
 proposedLogMarginalLikelihood = -0.5*(N*log(2*pi)+logDetK+abc);
 if log(rand) < proposedLogMarginalLikelihood-logLik
     newLogLik = proposedLogMarginalLikelihood;
     accepted = 1;
     
 end

end

