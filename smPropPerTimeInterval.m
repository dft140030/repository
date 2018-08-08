function [smPropPerTimeInt,threshMet] = smPropPerTimeInterval(utrackPackPaths, ...
    smConditions,qfsmPackPaths)

%% SMPROPPERTIMEINTERVAL 
%  reads data from U-Track folders and calculates spatiotemporal properties
%  of single-molecule data given length of frames to perform the calculations 
%
%  SYNOPSIS: [smPropPerTimeInt] = smPropPerTimeInterval(utrackPackPaths, ...
%                windAvg,divFlag,smConditions,qfsmPackPaths)
%
% INPUT: utrackPackPaths: vector of paths leading to U-Track Package
%                         First column are paths
%
%           smConditions: vector parameter set for this function. Each
%                         element (movie) should contain some values
%           
%            .rad_density: radius used for detecting neighboring
%                          speckles
%
%           .perturbation: = 0 : speckles under control
%                          = 1 : speckles influenced by Latrunculin-A
%                          = 2 : speckles recovering from Latrunculin-A
%
%                .maskLoc: = 0 : there do not exist any masks
%                          = 1 : there exists a folder with external
%                                masks within path provided in
%                                qfsmPackPaths (make sure folder name
%                                follows a name in extlroi variable.
%                          = 2 : there exists a folder with refined
%                                masks within the QFSMPackage folder
%                                provided in qfsmPackPaths
%                          = 3 : there exists a folder with refined
%                                masks within the SegmentationPackage
%                                folder provided in qfsmPackPaths
%
%               .windAvg: vector of values to average single molecule imaging
%                         data           
%
%               .divFlag: vector of flags to perform calculation of properties
%                         = -1 : link SM frames to earlier FSM frames
%                         =  0 : link SM frames to center FSM frames
%                         =  1 : link SM frames to later FSM frames
%
%                .minLft: minimum acceptable lifetime for single-molecules
%                         for analysis
%                .minObs: minimum acceptable existing frames for
%                         single-molecules
%          qfsmPackPaths: vector of paths leading to QFSM Package
%                         only input if sm should be within qfsm mask
%
% OUTPUT:  smPropPerTimeInt: cell array of vectors of structures. Cell
%                            array is (number of movies by 1), number of
%                            structures in vector is (windAvg by 1).
%                            Has the following fields
%
%   .smFrames        : List of single molecule frames associated with FSM
%                      frame n.
%   .particleList    : List of particles that exist in the window
%                      for time-averaging
%   .meanPos         : Mean position (x/y-coordinates) per particle in
%                      the window for time-averaging.
%   .meanDisp        : Mean frame-to-frame displacement per particle in FSM
%                      frame n.
%   .lifetime        : final frame of the last observation for a SM track 
%                      minus the first frame of the first observation for
%                      a SM track within a particular interval.
%   .numObservations : number of observations for a SM track within a 
%                      particular interval.
%   .netVelocity     : SM velocity vector  per frame averaged
%                      across a particular interval.
%   .meanAmp         : mean intensity of a SM track during a particular
%                      window of time
%   .netSpeed        : the net speed for a particular SM track during a
%                      particular window of time
%   .netVelAngle     : the net velocity angle within a particular interval
%                      for a particular SM track with respect to a 
%                      horizontal axis crossing the initial position of the SM
%   .trackClass      : classification of all tracks based on moment scaling
%                      spectrum analysis
%                      = 0 : immobile
%                      = 1 : confined brownian
%                      = 2 : pure brownian (free diffusion)
%                      = 3 : directed motion
%   .diffCoef        : diffusion coefficient with row number corresponding
%                      to particle index
%   .confRad         : confinement radius of immobile, confined particle with row 
%                      number corresponding to particle index
%
% Deryl Tschoerner, February 2018

%% definitions

%number of sm movies being processed
numMovs = length(utrackPackPaths);
if isempty(utrackPackPaths)
    error('error:: single-particle tracking path was not found')
end
fs = filesep;
smPropPerTimeInt = cell(numMovs,1);
threshMet = cell(numMovs,1);
numFramesFSM = nan(numMovs,1);
numFramesSM = nan(numMovs,1);
extlroi = ["external_roi","external_mask1","External_mask","External_ROI","external_ROI", ...
    "External_Segmentation"];

for i = 1 : numMovs
    
    if isempty(utrackPackPaths{i,1})
       continue
    end
    
    %optional: qfsm map extraction to exclude potentially irrelevant sm
    if ~isempty(qfsmPackPaths{i,1}) && smConditions.maskLoc(i) ~= 0
        
        cd(qfsmPackPaths{i,1})
        if smConditions.maskLoc(i) == 1
            
            masks = dir([qfsmPackPaths{i,1} fs char(extlroi(isfolder(extlroi)))]);            
        elseif smConditions.maskLoc(i) == 2
            
            masks = dir([qfsmPackPaths{i,1} fs 'QFSMPackage' fs 'refined_masks' fs 'refined_masks_for_channel_1']);
        elseif smConditions.maskLoc(i) == 3
            
            masks = dir([qfsmPackPaths{i,1} fs 'SegmentationPackage' fs 'refined_masks' fs 'refined_masks_for_channel_1']);
        else
            
            disp('error:: Folder of masks does not exist.')
        end
    end
    
    for ii = 1 : 2
        
        if ii == 1
            trormo = 'Motion';
        else
            trormo = 'Track';
        end
        
        try
            load([utrackPackPaths{i,1} fs trormo 'Analysis' fs ...
                'channel_1.mat'], 'tracks','diffAnalysisRes');
            continue
        end        
    end
    
    trackedFeatureInfo = convStruct2MatIgnoreMS(tracks);
    classifyDiffMat = [vertcat(diffAnalysisRes.classification), ...
        catStruct(1,'diffAnalysisRes.fullDim.normDiffCoef'), ...
        catStruct(1,'diffAnalysisRes.confRadInfo.confRadius')];

    %rows = all y/x-coords from tracks, columns=frames
    xyMat(:,:,1) = trackedFeatureInfo(:,1 : 8 : end);
    xyMat(:,:,2) = trackedFeatureInfo(:,2 : 8 : end);
    ampMat = trackedFeatureInfo(:,4 : 8 : end);
    
    numFramesSM(i) = size(xyMat,2);
    numFramesFSM(i) = floor(numFramesSM(i) / smConditions.windAvg(i) + 1);
    smPropPerTimeInt{i,1} = struct;

    %% getting .smFrames
    %assigning which SM frame(s) belong to each FSM frame and which particle(s)
    %belong to which interval of frames - depending on divFlag input
    if smConditions.divFlag(i) == 0
        
        %to divide an interval for centering FSM frame
        cHalf = ceil(smConditions.windAvg(i) / 2);

        %the first interval will begin with the first FSM track
        windFirst = 1;

        %one more interval for centered
        windAdd = 1;

        %tail interval
        smPropPerTimeInt{i,1}(numFramesFSM(i),1).smFrames = ...
            (numFramesSM(i) - cHalf + 1 : numFramesSM(i))';

        %head interval
        smPropPerTimeInt{i,1}(1,1).smFrames = (1 : cHalf + 1)';

        %middle intervals
        for j = 2 : numFramesFSM(i) - 1
            smPropPerTimeInt{i,1}(j,1).smFrames = ...
            ((j - 2)*smConditions.windAvg(i) + cHalf + 1: (j - 1)*smConditions.windAvg(i) + cHalf + 1)';
        end
    else
        %no interval increment if not centered
        windAdd = 0;

        %the first interval will begin with the first or second FSM track; add
        %to this block if divFlag values changes
        if smConditions.divFlag(i) == -1, windFirst = 1; else windFirst = 2; end

        %memory allocation for frame intervals, ESPECIALLY since this final
        %interval will have an extra frame and could not be initialized below
        smPropPerTimeInt{i,1}(windFirst + (numFramesFSM(i) - 2),1).smFrames = ...
            ((numFramesFSM(i) - 2)*smConditions.windAvg(i) + 1 : ...
            (numFramesFSM(i) - 1)*smConditions.windAvg(i))';

        %initializing all other intervals
        for j = windFirst : windFirst + (numFramesFSM(i) - 3)
            smPropPerTimeInt{i,1}(j,1).smFrames = ...
            ((j - windFirst)*smConditions.windAvg(i) + 1 : ...
            (j + 1 - windFirst)*smConditions.windAvg(i) + 1)';
        end
    end
    
    %% getting .meanPos, .particleList, .meanDisp, .lifeTime, .numObservations,
    %% .netVelocity, .meanAmp, .netVelAngle, .netSpeed, .trackClass, .diffCoef,
    %% and .confRad
    windLast = windFirst + windAdd + numFramesFSM(i) - 2;
    maskCoord = cell(windLast - windFirst + 1,1);
    
    for j = windFirst : windLast
        
        if smConditions.maskLoc(i) ~= 0
            
            %extract (y,x) coordinates for all pixels within refined, external,
            %or segmented mask, respectively.
            if smConditions.maskLoc(i) == 1

                [maskCoord{j,1}(:,2),maskCoord{j,1}(:,1)] = ...
                    find(imread([qfsmPackPaths{i,1} fs char(extlroi(isfolder(extlroi))) fs ...
                    masks(j + 2).name]));
            elseif smConditions.maskLoc(i) == 2

                [maskCoord{j,1}(:,2),maskCoord{j,1}(:,1)] = ...
                    find(imread([qfsmPackPaths{i,1} fs 'QFSMPackage' fs 'refined_masks' ...
                    fs 'refined_masks_for_channel_1' fs masks(j + 2).name]));
            elseif smConditions.maskLoc(i) == 3

                [maskCoord{j,1}(:,2),maskCoord{j,1}(:,1)] = ...
                    find(imread([qfsmPackPaths{i,1} fs 'SegmentationPackage' fs 'refined_masks' ...
                    fs 'refined_masks_for_channel_1' fs masks(j + 2).name]));
            end
        end
        %particular interval of (x,y) coordinate and amplitude matrices
        xyInterval(:,1 : length(smPropPerTimeInt{i,1}(j,1).smFrames),1) = ...
            xyMat(:,smPropPerTimeInt{i,1}(j,1).smFrames,1);
        xyInterval(:,1 : length(smPropPerTimeInt{i,1}(j,1).smFrames),2) = ...
            xyMat(:,smPropPerTimeInt{i,1}(j,1).smFrames,2);
        
        %take average positions and amplitude of particles within a particular 
        %interval
        xyMean = [nanmean(xyInterval(:,:,1),2), nanmean(xyInterval(:,:,2),2)];
        ampMean = nanmean(ampMat(:,smPropPerTimeInt{i,1}(j,1).smFrames),2);
        
        %vector of existing particle indices within interval
        indxNoNan = find(~isnan(xyMean(:,1)));
        smPropPerTimeInt{i,1}(j,1).smList = indxNoNan;
        
        %store motion classification for every particle existing with interval
        %of frames
        smPropPerTimeInt{i,1}(j,1).trackClass = ...
            classifyDiffMat(indxNoNan,2);
        
        %store diffusion coefficient for sm
        smPropPerTimeInt{i,1}(j,1).diffCoef = ...
            classifyDiffMat(indxNoNan,4);
        
        %store confinement radius for confined sm
        smPropPerTimeInt{i,1}(j,1).confRad = ...
            classifyDiffMat(indxNoNan,5);
        
        %store only average position of those particles with relevance to the
        %particular interval. in (x,y) coordinates
        smPropPerTimeInt{i,1}(j,1).meanPos = [xyMean(indxNoNan,1), xyMean(indxNoNan,2)];
        smPropPerTimeInt{i,1}(j,1).meanAmp = ampMean(indxNoNan);
        
        %take the mean displacement
        smPropPerTimeInt{i,1}(j,1).meanDisp = ...
            nanmean(sqrt((xyInterval(indxNoNan,2 : end,1) - ...
            xyInterval(indxNoNan,1 : end - 1,1)) .^ 2 + ...
            (xyInterval(indxNoNan,2 : end,2) - ...
            xyInterval(indxNoNan,1 : end - 1,2)) .^ 2),2);
        
        %matrix of existing SM tracks within a particular interval
        existTracks = xyInterval(indxNoNan,:,1);
        
        distMats = createDistanceMatrix(smPropPerTimeInt{i}(j).meanPos, ...
            smPropPerTimeInt{i}(j).meanPos);
        
        for iCol = 1 : size(distMats,2)

            mtchindx = find(distMats(:,iCol) <= smConditions.rad_density(i));
            smPropPerTimeInt{i}(j).smDensity(iCol,1) = length(mtchindx);
        end
        
        %finding the lifetime, number of observations, and net velocity of a SM
        %within a particular interval
        for k = 1 : length(indxNoNan)
            
            %those SM frames that a particular SM track from matrix of existing
            %tracks within a particular interval was found to exist in
            temp = find(~isnan(existTracks(k,:)));
            
            %the lifetime for an existing SM track within a particular interval
            smPropPerTimeInt{i,1}(j,1).lifetime(k,1) = max(temp) - min(temp) + 1;
            
            %the number of observations for an existing SM track within a 
            %particular interval
            smPropPerTimeInt{i,1}(j,1).numObservations(k,1) = ...
                sum(~isnan(existTracks(k,:)));
            
            %the net velocity for an existing SM track within a particular
            %interval along y/x-directions
            smPropPerTimeInt{i,1}(j,1).netVelocity(k,:) = ...
                (xyInterval(indxNoNan(k),max(temp),:) - ...
                xyInterval(indxNoNan(k),min(temp),:)) / ...
                (smPropPerTimeInt{i,1}(j,1).lifetime(k) - 1);
            
            %the net speed for an existing SM track within a particular
            %interval along the plane
            smPropPerTimeInt{i,1}(j,1).netSpeed(k,1) = ...
                sqrt(sum(smPropPerTimeInt{i,1}(j,1).netVelocity(k,:) .^ 2));
            
            %the net velocity angle for an existing SM track within a
            %particular interval.
            %NOTE: THESE ANGLES ARE TAKEN WITH RESPECT TO IMAGE COORDINATE
            %SYSTEM; i.e. SGN(ANGLE) => -SGN(ANGLE)
            smPropPerTimeInt{i,1}(j,1).netVelAngle(k,1) = ...
                atan(smPropPerTimeInt{i,1}(j,1).netVelocity(k,2) / ...
                smPropPerTimeInt{i,1}(j,1).netVelocity(k,1));
        end
        
        %SM thresholding survivors within the mask
        if isempty(indxNoNan) && smConditions.maskLoc(i) ~= 0
            
            threshMet{i,1}{j,1} = find(smPropPerTimeInt{i,1}(j,1).lifetime ...
                >= smConditions.minLft(i) & ismember(round(smPropPerTimeInt ...
                {i,1}(j,1).meanPos), [maskCoord{j,1}(:,1),maskCoord{j,1}(:,2)],'rows'));
            
            smPropPerTimeInt{i,1}(j,1).maskDensity = ...
                length(indxNoNan(threshMet{i,1}{j,1})) / length(maskCoord{j,1}(:,1));
        else %~isempty(indxNoNan) && smConditions.maskLoc(i) == 0
            
            threshMet{i,1}{j,1} = ...
                find(smPropPerTimeInt{i,1}(j,1).lifetime >= smConditions.minLft(i));
        end
    end
    clear xyMat xyInterval
end

end