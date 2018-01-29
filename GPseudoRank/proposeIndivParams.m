function proposedIndivParams = proposeIndivParams(currentLogParams,stepSize,nGenes)
%sampling parameters for pseudoRankIndivParams   
proposedIndivParams = currentLogParams + stepSize.*randn([3 nGenes]);
end