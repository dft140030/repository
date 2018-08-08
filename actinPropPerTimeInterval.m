function [actinPropPerTimeInt] = actinPropPerTimeInterval(qfsmPackPaths,actinConditions)

%% ACTINPROPPERTIMEINTERVAL
%  reads data from QFSM Package and calculates spatiotemporal properties of
%  speckle data per frame
%
%  SYNOPSIS: [actinPropPerTimeInt] = actinPropPerTimeInterval(qfsmPackPaths)
%
%     INPUT: qfsmPackPaths: vector of paths leading to QFSM Package
%
%          actinConditions: vector parameter set for this function. Each
%                           element (movie) should contain some values
%
%              .rad_density: radius used for detecting neighboring
%                            speckles
%
%                  .locmask: = 0 : there exists a folder with external
%                                  masks within path provided in
%                                  qfsmPackPaths (make sure folder name
%                                  follows a name in extlroi variable.
%                            = 1 : there exists a folder with refined
%                                  masks within the QFSMPackage folder
%                                  provided in qfsmPackPaths
%                            = 2 : there exists a folder with refined
%                                  masks within the SegmentationPackage
%                                  folder provided in qfsmPackPaths
%
%                  .maskLoc: = 0 : the mask location is in an external roi
%                                  folder
%                            = 1 : the mask location is in the QFSMPackage
%                                  folder as a refined mask
%                            = 2 : the mask location is in the
%                                  SegmentationPackage folder as a refined mask
%
%    OUTPUT: actinPropPerTimeInt: Structure array with number of entries =
%                                 fsmFrames. For each FSM Frame n, it will
%                                 have the following fields:
%
%         .speckleList      : List of speckles that exist in the frame
%
%         .kinScore         : The kinetic score from the qFSM package                            
%
%         .kinScorePos      : The positions at which we have kinetic scores
%                             (units: pixels)
%
%         .speckleInitPos   : The speckle position for a particular FSM 
%                             frame (units: pixels)
%
%         .speckleMidPos    : The speckle midpoint position for a 
%                             particular FSM frame (units: pixels)
%
%         .speckleVelocity  : The speckle velocity from a frame to the
%                             following FSM frame (units: pixels/frame)
%
%         .speckleSpeed     : The speckle speed from particular FSM frame
%                             to the following FSM frame 
%                             (units: pixels/frame)
%
%         .speckleMvmtCohere: The average comovement between a speckle and
%                             neighboring speckles
%
%         .maskDensity      : The number of speckles within a mask divided
%                             by the number of pixels comprising the mask
%
%         .indxMask         : Indices of the mask obtained from qFSM package
%
% Deryl Tschoerner, February 2018

%% definitions
numMovs = length(qfsmPackPaths);
fs = filesep;
maps = cell(numMovs,1);
actinPropPerTimeInt = cell(numMovs,1);
fsmFrames = nan(numMovs,1);

%naming possible combinations of external masking names which were used in
%the past. Fourth option may be applied most often as of June 2018
extlroi = ["external_roi","external_mask1","External_mask","External_ROI","external_ROI", ...
    "External_Segmentation"];

for i = 1 : numMovs
    
    %field names for the variables found in the kineticAnalysis folder
    maps{i,1} = dir(strcat([qfsmPackPaths{i,1} fs 'QFSMPackage' fs 'kineticAnalysis' fs '*.mat'])); 
    fsmFrames(i,1) = length(maps{i,1});
    
    %skip an element in the cell if path does not exist. 
    if isempty(qfsmPackPaths{i,1}) 
       continue
    end
    
    %given path, loading speckle tracks
    load([qfsmPackPaths{i,1} fs 'QFSMPackage' fs 'speckleTracks' fs '_tracks.mat'], 'MPM');
    
    %assign zeros to NaN and retire the position information of the first
    %frame. we only want the first fsm frame for the derivative data!!
    MPM(MPM == 0) = NaN; MPM(:,end + 1 : end + 2) = NaN; % #ok<NODEF>
    
    %store field names for readout of values from the structure fields found
    %in 'maps' variable
    temp = repmat(struct('polyMap',[],'kinScore',[],'depolyMap',[], ...
        'kinMap2C',[]),fsmFrames(i,1),1);
    
    %field names for vector of structures holding speckle properties
    actinPropPerTimeInt{i,1} = repmat(struct('kinScorePos',[],'kinScore',[], ...
        'speckleList',[],'speckleInitPos',[],'speckleMidPos',[], ...
        'speckleVelocity',[],'speckleSpeed',[],'speckleDensity',[] , ...
        'speckleMvmtCohere',[]),fsmFrames(i,1),1);
    
    %if there is a non-empty element, then readout the masks
    if ~isempty(qfsmPackPaths{i,1}) && actinConditions.maskLoc(i) ~= 0
        
        %change directory to path to look for folder
        cd(qfsmPackPaths{i,1})
        
        if actinConditions.maskLoc(i) == 1
            
            masks = dir([qfsmPackPaths{i,1} fs char(extlroi(isfolder(extlroi)))]);            
        elseif actinConditions.maskLoc(i) == 2
            
            masks = dir([qfsmPackPaths{i,1} fs 'QFSMPackage' fs 'refined_masks' fs 'refined_masks_for_channel_1']);
        elseif actinConditions.maskLoc(i) == 3
            
            masks = dir([qfsmPackPaths{i,1} fs 'SegmentationPackage' fs 'refined_masks' fs 'refined_masks_for_channel_1']);    
        else
            
            disp('error:: Folder of masks does not exist.')
        end
    end
    %only retain movieInfo features that exist within mask
    
    %% get .kinScorePos, .kinScore, .specklePosInit, .speckleVelocity, and
    %% .speckleSpeed
    for j = 1 : fsmFrames(i,1)
        
        %assigning speckle list, initial positions, velocities, and speeds
        actinPropPerTimeInt{i,1}(j,1).speckleList = find(~isnan(MPM(:,2*j)));
        
        actinPropPerTimeInt{i,1}(j,1).speckleInitPos = ...
            [MPM(actinPropPerTimeInt{i,1}(j,1).speckleList,2*j), ...
            MPM(actinPropPerTimeInt{i,1}(j,1).speckleList,2*j - 1)];
        
        if actinConditions.maskLoc(i) ~= 0
            
            if actinConditions.maskLoc(i) == 1

                %extract (y,x) coordinates for all pixels within external mask
                [actinPropPerTimeInt{i,1}(j).indxMask(:,2), ...
                    actinPropPerTimeInt{i,1}(j).indxMask(:,1)] = ...
                    find(imread([qfsmPackPaths{i,1} fs char(extlroi(isfolder(extlroi))) fs ...
                    masks(j + 2).name]));

            elseif actinConditions.maskLoc(i) == 2

               %extract (y,x) coordinates for all pixels within refined mask
               [actinPropPerTimeInt{i,1}(j).indxMask(:,2), ...
                   actinPropPerTimeInt{i,1}(j).indxMask(:,1)] = ...
                   find(imread([qfsmPackPaths{i,1} fs 'QFSMPackage' fs 'refined_masks' ...
                   fs 'refined_masks_for_channel_1' fs masks(j + 2).name]));

            elseif actinConditions.maskLoc(i) == 3

               [actinPropPerTimeInt{i,1}(j).indxMask(:,2), ...
                   actinPropPerTimeInt{i,1}(j).indxMask(:,1)] = ...
                   find(imread([qfsmPackPaths{i,1} fs 'SegmentationPackage' fs 'refined_masks' ...
                   fs 'refined_masks_for_channel_1' fs masks(j + 2).name]));
            end        
            
            %indices of positions
            pos = [actinPropPerTimeInt{i,1}(j).speckleInitPos(:,2) ...
                actinPropPerTimeInt{i,1}(j).speckleInitPos(:,1)];

            actinPropPerTimeInt{i,1}(j).speckleMaskDensity = ...
                sum(ismember(round(pos),[actinPropPerTimeInt{i,1}(j).indxMask(:,2), ...
                actinPropPerTimeInt{i,1}(j).indxMask(:,1)],'rows')) / ...
                length(actinPropPerTimeInt{i,1}(j).indxMask(:,1));
        end
        
        %polluting speckle positions with noise since positions are in
        %integers
        nextpos = [MPM(actinPropPerTimeInt{i,1}(j,1).speckleList,2*j + 2), ...
            MPM(actinPropPerTimeInt{i,1}(j,1).speckleList,2*j + 1)];
        
        actinPropPerTimeInt{i,1}(j,1).speckleInitPos = ...
            actinPropPerTimeInt{i,1}(j,1).speckleInitPos + ...
            normrnd(0,15/90,size(actinPropPerTimeInt{i,1}(j,1).speckleInitPos,1),2);
        
        %vector of directories
        temp(j) = ...
            load(strcat([qfsmPackPaths{i,1} fs 'QFSMPackage' fs 'kineticAnalysis' fs maps{i}(j).name]));
        
        %assigning kinetic scores and kinetic positions
        actinPropPerTimeInt{i,1}(j,1).kinScore = temp(j).kinScore(:,4);
        actinPropPerTimeInt{i,1}(j,1).kinScorePos = ...
            [temp(j).kinScore(:,3), temp(j).kinScore(:,2)];
        
        %speckle velocities and speeds
        actinPropPerTimeInt{i,1}(j,1).speckleVelocity = ...
            [MPM(actinPropPerTimeInt{i,1}(j,1).speckleList,2*j + 2) - ...
            actinPropPerTimeInt{i,1}(j,1).speckleInitPos(:,1), ...
            MPM(actinPropPerTimeInt{i,1}(j,1).speckleList,2*j + 1) - ...
            actinPropPerTimeInt{i,1}(j,1).speckleInitPos(:,2)];
        
        actinPropPerTimeInt{i,1}(j,1).speckleSpeed = ...
            sqrt(sum(actinPropPerTimeInt{i,1}(j,1).speckleVelocity .^ 2,2));
        
        %assigning speckle midpoint positions
        actinPropPerTimeInt{i,1}(j,1).speckleMidPos = ...
            actinPropPerTimeInt{i,1}(j,1).speckleInitPos + ...
            (actinPropPerTimeInt{i,1}(j,1).speckleVelocity / 2); 
        
        %distances of all speckles between eachother for a given frame of speckles
        distMats = createDistanceMatrix(actinPropPerTimeInt{i}(j).speckleInitPos, ...
            actinPropPerTimeInt{i}(j).speckleInitPos);
        
        for iCol = 1 : size(distMats,2)
            
            actinPropPerTimeInt{i}(j).speckleDensity(iCol,1) = ...
                sum(distMats(:,iCol) <= actinConditions.rad_density(i));
            matchindx = find(distMats(:,iCol) <= actinConditions.rad_density(i));
            matchindx = matchindx(matchindx ~= iCol);
            temp2 = nan(1,length(matchindx));
            if isempty(matchindx)
                
                actinPropPerTimeInt{i}(j).speckleMvmtCohere(iCol,1) = NaN;
                continue
            end
            
            for iAng = 1 : length(matchindx)
                
                vel{1} = actinPropPerTimeInt{i}(j).speckleVelocity(matchindx(iAng),:);
                vel{2} = actinPropPerTimeInt{i}(j).speckleVelocity(iCol,:);
                temp2(iAng) = acos(dot(vel{1},vel{2}) / (norm(vel{1}) * norm(vel{2})));               
            end
            
            actinPropPerTimeInt{i}(j).speckleMvmtCohere(iCol,1) = nanmean(real(temp2));
        end            
    end
end

end