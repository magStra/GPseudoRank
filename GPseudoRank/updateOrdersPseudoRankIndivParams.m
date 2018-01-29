function [newLks newLogLik  accepted] = ...
    updateOrdersPseudoRankIndivParams(logParams, proposedTimeDiffs,dataProposedOrder,logLik, nTimes, nGenes)
accepted = 0;
newLogLik = logLik;
newLks = NaN;
opts.LT     = true;
propLogLiks     = zeros(1,nGenes);
if nCores > 1
    parfor j = 1:nGenes
    L           = exp(logParams(2,j));
    sf2         = exp(2*logParams(1,j));
    se2         = exp(2*logParams(3,j));
    K           = sf2 * exp(-proposedTimeDiffs.^2/(2*L^2)) ...
        + se2 * eye(nTimes);
  
  
    KChol       = chol(K);
    z           = linsolve(KChol',dataProposedOrder(j,:)',opts);
    logDetK     =  2*(sum(log(diag(KChol))));
    abc = z'*z;
    propLogLiks(j) =  -0.5*(nTimes*log(2*pi)+logDetK+abc);
    end 
    else
    for j = 1:nGenes
    L           = exp(logParams(2,j));
    sf2         = exp(2*logParams(1,j));
    se2         = exp(2*logParams(3,j));
    
    K           = sf2 * exp(-proposedTimeDiffs.^2/(2*L^2)) ...
        + se2 * eye(nTimes);
  
    KChol       = chol(K);
    z           = linsolve(KChol',dataProposedOrder(j,:)',opts);
    logDetK     =  2*(sum(log(diag(KChol))));
    abc = z'*z;
    propLogLiks(j) =  -0.5*(nTimes*log(2*pi)+logDetK+abc);
    end 
end

 proposedLogMarginalLikelihood = sum(propLogLiks);

 if log(rand) < proposedLogMarginalLikelihood-logLik
     newLogLik = proposedLogMarginalLikelihood;
     newLks    = propLogLiks;
     accepted = 1;
     
 end

end

