function [newLogLik newLogPriorOfLogHypers newLogHypers newSquaredHypers accepted newKInv newLogDet] = ...
    updateHypersPseudoRank(proposedLogHypers,currentSquaredHypers, timeDiffs,permData,logLik,...
    logPriorOfLogHypers,hyperPriorParams, nTimes, nGenes,varData,sparse)
N = nGenes*nTimes;
newKInv = NaN;
newLogDet = NaN;
accepted = zeros(1,2);
proposedSquaredHypers           = exp(2* proposedLogHypers);
newSquaredHypers    = currentSquaredHypers;
newLogLik = logLik;
newLogPriorOfLogHypers = logPriorOfLogHypers;
for j = 1:2
squaredHypers       = newSquaredHypers;
squaredHypers(j)    = proposedSquaredHypers(j);
sf2         = squaredHypers(1);
l2          = squaredHypers(2);
se2         = varData - squaredHypers(1);

K           = sf2*exp(timeDiffs/(2*l2)) + se2*eye(nTimes);

%compute the inverse and the determinant

KChol       = chol(K);
opts.LT     = false;
opts.UT     = true;
X           = linsolve(KChol,eye(nTimes),opts);
KInv        = X*X';
logDetK     =  2*(sum(log(diag(KChol))));

if(~isreal(logDetK))  %We should not ever enter here, but just in case
    disp('Sampling orders: Numerical error - covariance matrix may not be positive definite')
end
%computing the proposed logPriorOfLogHypers
proposedLogPriorOfLogHypers = newLogPriorOfLogHypers;
proposedLogPriorOfLogHypers(j) = log(normpdf(proposedLogHypers(j),hyperPriorParams(j,1),...
                    hyperPriorParams(j,2)));
%computing the proposed logLikelihood
 abc = 0;
    for i = 1:nGenes
        abc = abc + permData(i,:)*KInv*permData(i,:)';
    end
   
    
 proposedLogMarginalLikelihood = -0.5*(N*log(2*pi)+logDetK+abc);
 
 if log(rand) < (proposedLogMarginalLikelihood + sum(proposedLogPriorOfLogHypers) - ...
         newLogLik - sum(newLogPriorOfLogHypers))
     newLogLik = proposedLogMarginalLikelihood;
     
     newLogPriorOfLogHypers = proposedLogPriorOfLogHypers;
     accepted(j) = 1;
      newSquaredHypers(j)    = squaredHypers(j);
      newKInv = KInv;
      newLogDet = logDetK;
    
 end
end
newLogHypers = 0.5*log(newSquaredHypers);

end

