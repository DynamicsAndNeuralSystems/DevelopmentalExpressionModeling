numRepeats = 5;

numTimePoints = 7;
lambdaEst = zeros(numRepeats,numTimePoints);
strengthEst = zeros(numRepeats,numTimePoints);
offSetEst = zeros(numRepeats,numTimePoints);

for i = 1:numRepeats
    fprintf(1,'~~~~~~~ITERATION %u~~~~~~~\n',i);
    [maxExtent,lambdaEst(i,:),strengthEst(i,:),offSetEst(i,:)] = simulateGrid(false);
end

%-------------------------------------------------------------------------------

f = figure('color','w');
subplot(3,1,1);
errorbar(maxExtent,mean(lambdaEst),std(lambdaEst),'k')
title('Spatial scale')
subplot(3,1,2);
errorbar(maxExtent,mean(strengthEst),std(strengthEst),'k')
title('Strength')
subplot(3,1,3);
errorbar(maxExtent,mean(offSetEst),std(offSetEst),'k')
title('Offset')
xlabel('Max Distance')
