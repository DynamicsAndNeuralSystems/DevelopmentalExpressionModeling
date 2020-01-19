numRepeats = 10;

params = GiveMeDefaultParams();

numTimePoints = 7;
lambdaEst = zeros(numRepeats,numTimePoints);
lambdaErr = zeros(numRepeats,numTimePoints);
strengthEst = zeros(numRepeats,numTimePoints);
strengthErr = zeros(numRepeats,numTimePoints);
offSetEst = zeros(numRepeats,numTimePoints);
offSetErr = zeros(numRepeats,numTimePoints);

for i = 1:numRepeats
    fprintf(1,'~~~~~~~ITERATION %u~~~~~~~\n',i);
    [maxExtent,paramStruct] = simulateGrid(params,false);
    paramNames = coeffnames(paramStruct{1});

    % Mean estimates:
    lambdaEst(i,:) = cellfun(@(x)1/x.n,paramStruct);
    strengthEst(i,:) = cellfun(@(x)x.A,paramStruct);
    offSetEst(i,:) = cellfun(@(x)x.B,paramStruct);

    % CIs of estimates:
    CIs = cellfun(@(x)confint(x),paramStruct,'UniformOutput',false);
    % Invert for lambda:
    nIndex = strcmp(paramNames,'n');
    for j = 1:numTimePoints
        CIs{j}(:,nIndex) = 1./CIs{j}(:,nIndex);
    end
    % Store errors:
    strengthErr(i,:) = cellfun(@(x)diff(x(:,1))/2,CIs);
    offSetErr(i,:) = cellfun(@(x)diff(x(:,2))/2,CIs);
    lambdaErr(i,:) = cellfun(@(x)diff(x(:,3))/2,CIs);
end


%-------------------------------------------------------------------------------
f = figure('color','w');
myColors = BF_getcmap('dark2',7,1);
for i = 1:3
    subplot(1,3,i); hold('on')
    switch i
    case 1
        paramEst = lambdaEst;
        paramErr = lambdaErr;
        ylabel('Spatial scale, \lambda')
    case 2
        paramEst = strengthEst;
        paramErr = strengthErr;
        ylabel('Strength, A')
    case 3
        paramEst = offSetEst;
        paramErr = offSetErr;
        ylabel('Offset, f_0')
    end

    % Fit linear regression:
    if i==1
        ft = fittype('c*x');
        [c,Stats] = fit(maxExtent,mean(paramEst)',ft);
        f_handle = @(x) c.c*x;
    end
    xRange = linspace(0,maxExtent(end));
    plot(xRange,f_handle(xRange),':k')

    % Plot (colored) parameter variation
    for t = 1:numTimePoints
        errorbar(maxExtent(t),mean(paramEst(:,t)),mean(paramErr(:,t)),'o','color',myColors{t},'LineWidth',1.8)
        % if i==1 | i==2
        % Add empirical data:
        load('parameterFits.mat','paramMeanValues','paramErrValues');
        empirical = paramMeanValues(i,:)';
        empiricalErrs = paramErrValues(i,:)';
        smallOffset = 0.4;
        errorbar(maxExtent(t)+smallOffset,empirical(t),empiricalErrs(t),'o',...
                            'color',brighten(myColors{t},0.5),'LineWidth',1.8)
        % end
    end

    xlabel('Brain size, dmax (mm)')
end
f.Position = [1000        1116         831         222];
