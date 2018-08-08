function subpopPlot(spklData,smPropPerTimeInt,actinPropPerTimeInt,combinePropPerTimeInt) 
%% SUBPOPPLOT 
%  plot different subpopulations of single-molecules in comparison to speckles
%  for a particular frame
% Deryl Tschoerner, May 2018
FRAME = 3;
for i = 1;

    figure, imshow(spklData)
    hold on

    plot(smPropPerTimeInt{i,1}(FRAME).meanPos(combinePropPerTimeInt{i,1}(FRAME).smThresh ...
        (combinePropPerTimeInt{i,1}(FRAME).distSM2Speckle(:,1) <= 5,1),1), ...
        smPropPerTimeInt{i,1}(FRAME).meanPos(combinePropPerTimeInt{i,1}(FRAME).smThresh ...
        (combinePropPerTimeInt{i,1}(FRAME).distSM2Speckle(:,1) <= 5,1),2),'go')
    plot(smPropPerTimeInt{i,1}(FRAME).meanPos(combinePropPerTimeInt{i,1}(FRAME).smThresh ...
        (combinePropPerTimeInt{i,1}(FRAME).distSM2Speckle(:,1) > 5,1),1), ...
        smPropPerTimeInt{i,1}(FRAME).meanPos(combinePropPerTimeInt{i,1}(FRAME).smThresh ...
        (combinePropPerTimeInt{i,1}(FRAME).distSM2Speckle(:,1) > 5,1),2),'mo')    
    plot(actinPropPerTimeInt{i,1}(FRAME).speckleInitPos(:,1), ...
        actinPropPerTimeInt{i,1}(FRAME).speckleInitPos(:,2),'cs')
end

end