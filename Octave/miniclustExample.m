%simulate data set with 5000 cells
logSw2 = randn(1,100)*0.1 + log(sqrt(1));
sw2 = exp(2*logSw2);
logL = randn(1,100).*0.1+log(0.4);
l = exp(logL);
logSe2 = randn(1,100)*0.1 + log(sqrt(0.5));
se2 = exp(2*logSe2);
ppar = [logSw2;logL;logSe2];
simGPIndiv(2000,sw2,l,se2,'sim_Example1.csv','tau_Example1.csv',false,'uniform');
%assuming there are capture times
captureTimes = [zeros(1,500),ones(1,500),repmat(2,1,500),repmat(3,1,500)];
miniclust('sim_Example1.csv',captureTimes);

%now find default parameters
[n0,n3,n3a,jj,kk,delta,pp]  = setDefaultParams('sim_Example1KM.csv',true);

delete 'sim_Example1.csv'
delete 'sim_Example1KM.csv'
delete 'tau_Example1.csv'
