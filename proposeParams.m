function proposedParams = proposeParams(currentLogParams,stepSize,varData)
%sampling parameters for pseudoRank
    proposedParams = currentLogParams + (stepSize.*randn([1 2]))';
    if proposedParams(1) > log(sqrt(varData))
        proposedParams(1) = log(sqrt(varData) - 0.05);
    end
end