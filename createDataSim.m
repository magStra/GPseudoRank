function [] = createDataSim()
% creates the data for the simulations
%The pseudotime points are randomly sampled from U(0,1) distributions.
%This code generates a total of 32 simulated data sets, 16 with lower noise
%level, and 16 with higher noise level.
%output files: 
%'sim%d_paramsMeansBasic.csv': values of the GP parameters for the
%simulation with lower noise levels
%'sim%dBasic.csv': %d=1:16, simulated data with less noise, 50x90 matrix
%corresponding to 90 simulated cells and 50 simulated genes
%'tau%dBasic.csv': uniform simulated pseudotime points corresponding to
%simulated data set %d
%%%%%%%%%%%%%%%%%%%%%%%%%
%'sim%d_paramsMeansNoise.csv': values of the GP parameters for the
%simulation with higher noise levels
%'sim%dNoise.csv': %d=1:16, simulated data with less noise, 50x90 matrix
%corresponding to 90 simulated cells and 50 simulated genes
%'tau%dNoise.csv': uniform simulated pseudotime points corresponding to
%simulated data set %d

%basic 
for j = 1:16
    logSw2 = randn(1,50)*0.1 + log(sqrt(1));
    sw2 = exp(2*logSw2);
    logL = randn(1,50).*0.1+log(0.4);
    l = exp(logL);
    logSe2 = randn(1,50)*0.1 + log(sqrt(0.5));
    se2 = exp(2*logSe2);
    ppar = [logSw2;logL;logSe2];
    csvwrite(sprintf('sim%d_paramsMeansBasic.csv',j),ppar);
    %save the matrix of parameter values
    simGPIndiv(90,sw2,l,se2,sprintf('sim%dBasic.csv',j),...
        sprintf('tau%dBasic.csv',j),false,'uniform');
end
% more noise
for j = 1:16
    logSw2 = randn(1,50)*0.1 + log(sqrt(1));
    sw2 = exp(2*logSw2);
    logL = randn(1,50).*0.1+log(0.4);
    l = exp(logL);
    logSe2 = randn(1,50)*0.1 + log(sqrt(1));
    se2 = exp(2*logSe2);
    ppar = [logSw2;logL;logSe2];
    csvwrite(sprintf('sim%d_paramsMeansNoise.csv',j),ppar);
    simGPIndiv(90,sw2,l,se2,sprintf('sim%dNoise.csv',j),...
        sprintf('tau%dNoise.csv',j),false,'uniform');
end