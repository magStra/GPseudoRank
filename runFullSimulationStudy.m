
createDataSim();
parpool(4);%adapt to the number of cores
parfor j = 1:16
    j
    simulateAllCombs5(j);%less noisy data, adaptive moves during burn-in phase
end
parfor j = 1:16
    j
    simulateAllCombs4(j);%using same starting paths as parfor 5, moves not adaptive
end
parfor j = 1:16
    j
    simulateAllCombs44(j);%fewer capture times, moves not adapative
end
parfor j = 1:16
    j
    simulateAllCombs55(j);%fewer capture times, adaptive, same starting paths as 44
end
parfor j = 1:16
    j
    simulateAllCombs444(j);%more noise, moves not adapative
end
parfor j = 1:16
    j
    simulateAllCombs555(j);%more noise, adaptive, same starting paths as 444
end
poolobj = gcp('nocreate');
delete(poolobj);