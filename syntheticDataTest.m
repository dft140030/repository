function [betaOutput,syntheticProperties,testResults] = syntheticDataTest(actinPropPerTimeInt, ...
    betaInput,combineConditions,smactinFlag)
%
% SYNOPSIS: [betaOutput,syntheticProperties,testResults] = syntheticDataTest(actinPropPerTimeInt ...
%               betaInput,combineConditions,smactinFlag)
% Deryl Tschoerner, July 2018

numMovs = length(actinPropPerTimeInt);
syntheticProperties = cell(1,numMovs);
threshMet = cell(1,numMovs);
PROPORTION = 1;
for i = 1 : numMovs
    
    for j = 1 : length(actinPropPerTimeInt{i}) - 1
        
        numParticles = ceil(size(actinPropPerTimeInt{i}(j).speckleInitPos,1)*PROPORTION);
        syntheticProperties{i}(j).meanPos = actinPropPerTimeInt{i}(j).speckleInitPos(1:numParticles,:);
        
        syntheticProperties{i}(j).netSpeed = horzcat(ones(numParticles,1), ...
            actinPropPerTimeInt{i}(j).speckleSpeed(1:numParticles,:) ./ combineConditions.fsmInterval(i), ...
            actinPropPerTimeInt{i}(j).speckleDensity(1:numParticles,:)) * betaInput(:,1) + ...
            randn(size(syntheticProperties{i}(j).meanPos,1),1)*10^-4;
        
        syntheticProperties{i}(j).meanAmp = horzcat(ones(numParticles,1), ...
            actinPropPerTimeInt{i}(j).speckleSpeed(1:numParticles,:) ./ combineConditions.fsmInterval(i), ...
            actinPropPerTimeInt{i}(j).speckleDensity(1:numParticles,:)) * betaInput(:,2) + ...
            randn(size(syntheticProperties{i}(j).meanPos,1),1)*10^-4;
        
        syntheticProperties{i}(j).diffCoef = horzcat(ones(numParticles,1), ...
            actinPropPerTimeInt{i}(j).speckleSpeed(1:numParticles,:) ./ combineConditions.fsmInterval(i), ...
            actinPropPerTimeInt{i}(j).speckleDensity(1:numParticles,:)) * betaInput(:,3) + ...
            randn(size(syntheticProperties{i}(j).meanPos,1),1)*10^-4;
        
        threshMet{i}{j} = (1:length(syntheticProperties{i}(j).meanPos(1:numParticles,:)));
        
        %syntheticProperties{i}(j).netSpeed = syntheticProperties{i}(j).netSpeed(randsample(numParticles,numParticles));
        %syntheticProperties{i}(j).meanAmp = syntheticProperties{i}(j).meanAmp(randsample(numParticles,numParticles));
        %syntheticProperties{i}(j).diffCoef = syntheticProperties{i}(j).diffCoef(randsample(numParticles,numParticles));
    end
end

combineProperties = smactinPropPerTimeInt(syntheticProperties,actinPropPerTimeInt,combineConditions,threshMet,smactinFlag);
[tassm,matchindx] = smactinMat(combineProperties,combineConditions,smactinFlag);
[testResults,~,~] = smactinAggMov(tassm,matchindx,smactinFlag);
betaOutput = testResults{1}.mvrgs.coef;
end