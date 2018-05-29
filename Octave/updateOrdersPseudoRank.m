function [newLogLik  accepted] = ...
    updateOrdersPseudoRank(logParams, proposedTimeDiffs,dataProposedOrder,logLik, nTimes, nGenes)
accepted = 0;
newLogLik = logLik;
N = nTimes*nGenes;
L           = exp(logParams(2));
sf2         = exp(2*logParams(1));
se2         = exp(2*logParams(3));

K           = sf2 * exp(proposedTimeDiffs/(2*L^2)) ...
    + se2 * eye(nTimes);
KChol       = chol(K);
opts.LT     = false;
opts.UT     = true;
X           = linsolve(KChol,eye(nTimes),opts);
KInv        = X*X';
logDetK     =  2*(sum(log(diag(KChol))));

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

