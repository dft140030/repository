function [tassm,matchindx] = smactinMat(combinePropPerTimeInt,combineConditions,...
    smactinFlag)

%% SMACTINMAT
%  Converts some combinePropPerTimeInt fields into data matrix form
%
%  SYNOPSIS: [tassm,plotIndices] = smactinMat(combinePropPerTimeInt, ...
%                combineConditions)
%  INPUT: combinePropPerTimeInt: vector of structures containing various properties
%                                of both single molecules and actin (see
%                                smactinPropPerTimeInt for more details)%
%
%         combineConditions: vector of structures with various fields imposing
%                            restrictions pertaining to both sm, speckles,
%                            and ks (see smactinPropPerTimeInt for more details)
%
%               smactinFlag: structure of flag with various fields imposing
%                            restrictions on the computation (see
%                            smactinPropPerTimeInt)
%
%                       .combFlag: vector of integers from 1 to 9 pertaining to
%                                  the particular combination of properties that
%                                  should be calculated (might want to have a matrix
%                                  of these flags if using multiple
%                                  combinations per movie)
%                                  IMPORTANT: IF WANTING TO TAKE RECEPTORS
%                                  MATCHED TO BOTH KINETIC SCORE AND
%                                  SPECKLE NEAREST NEIGHBORS, THEN INPUT A
%                                  10.
%
%  OUTPUT: tassm: (total attributes speckles and single molecules): data 
%                 matrix of single-molecule and actin properties
%
%      matchindx: indices of sm,ks,and speckles meeting requirements from
%                 combineConditions
%
% REMARKS: Properties available for concatenation are (what gets concatenated 
%                     depends on combFlag:
%
%              kinetic score                kinetic score density
%
%              speckle speed                speckle density
%              speckle movement coherence
%
%              sm mean displacement         sm mean amplitude
%              sm diffusion coefficient     sm confinement radius*
%              sm density                   sm net speed   
%              sm comovement angle          * - Applicable only for immobile/confined sms
%
% Deryl Tschoerner, May 2018

numMovs = size(combinePropPerTimeInt,1);
tassm = cell(numMovs,1);
matchindx = tassm;

if smactinFlag.match == 1    
    iMatch = 1;
elseif smactinFlag.match == 2
    iMatch = 2;
else
    error('error:: invalid match parameter')
end

%change value depending on the combineProperties function as input
for iMov = 1 : numMovs
    
    wind1st = combineConditions.firstFsmFrame(iMov);
    if isempty(combinePropPerTimeInt{iMov})
        continue
    end

    for iComb = 1 : length(smactinFlag.combFlag(iMov)) %1 : size(combinePropPerTimeInt,2)
        
        switch smactinFlag.combFlag(iMov)
        case 3
            %concatenate nn: sm, object: speckle attributes
            tassm{iMov,iComb}{1} = [ ...
                vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).smMeanDisp) ... (:,1)
                vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).smNetSpeed) ... (:,2)
                vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).smMeanAmp) ...  (:,3)
                %vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).smDensity) ... (:,4)
                vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).diffCoef) ...   (:,5)
                vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).confRad)...     (:,6)
                vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).trackClass)];  %(:,7)
            
            tassm{iMov,iComb}{2} = [ ...
                vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).speckleSpeed) ...     (:,1)
                vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).matchAngle) ...       (:,2)
                vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).speckleDensity)]; ... (:,3)
            
        case 6
            %concatenate nn: sm, object: ks attributes
            tassm{iMov,iComb}{1} = [ ...
                vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).smMeanDisp) ... (:,1)
                vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).smNetSpeed) ... (:,2)
                vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).smMeanAmp) ...  (:,3)
                %vertcat(combinePropPerTimeInt{iMov,iComb}{iNN}(wind1st : end).smDensity) ...    (:,4)
                vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).diffCoef) ...   (:,5)
                vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).confRad)...     (:,6)
                vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).trackClass)];  %(:,7)
            
            tassm{iMov,iComb}{2} = [ ...
                vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).kinScore) ...       (:,1)
                vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).kinScoreDensity)]; %(:,2)
        
        case 2
            %concatenate nn: ks, object: actin attributes
            tassm{iMov,iComb}{1} = [ ...
                vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).kinScore) ...       (:,1)
                vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).kinScoreDensity)]; %(:,2)
            
            tassm{iMov,iComb}{2} = [ ...
                vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).speckleSpeed) ...     (:,1)
                vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).matchAngle) ...       (:,2)
                vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).speckleDensity)]; ... (:,3)
            
        case 8
            %concatenate nn: sm, object: ks attributes
            tassm{iMov,iComb}{1} = [ ...
                vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).kinScore)]; ...       (:,1)
                ...vertcat(combinePropPerTimeInt{iMov,iComb}{iNN}(wind1st : end).kinScoreDensity)];  %(:,2)
            
            tassm{iMov,iComb}{2} = [ ...
                ...vertcat(combinePropPerTimeInt{iMov,iComb}{iNN}(wind1st : end).smMeanDisp) ... (:,1)
                vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).smNetSpeed) ... (:,2)
                vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).smMeanAmp) ...  (:,3)
                ...vertcat(combinePropPerTimeInt{iMov,iComb}{iNN}(wind1st : end).smDensity) ...  (:,4)
                vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).diffCoef) ...   (:,5)
                vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).confRad)...     (:,6)
                vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).trackClass)];  %(:,7)
            
        case 4
            %concatenate nn: actin, object: ks attributes
            tassm{iMov,iComb}{1} = [ ...
                vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).speckleSpeed) ...  (:,1)
                vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).matchAngle) ...    (:,2)
                vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).speckleDensity)]; %(:,3)
            
            tassm{iMov,iComb}{2} = [ ...
                vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).kinScore) ...       (:,1)
                vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).kinScoreDensity)]; %(:,2)
            
        case 7 
            %concatenate nn: actin, object: sm attributes
            tassm{iMov,iComb}{1} = [ ...
                vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).speckleSpeed), ...         (:,1)
                ...vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).speckleMvmtCohere), ... (:,2)
                vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).speckleDensity)];         %(:,3)             
            
            tassm{iMov,iComb}{2} = [ ...
                vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).smNetSpeed), ...    (:,1)
                vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).smMeanAmp), ...     (:,2)
                vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).diffCoef)];        %(:,3)
            
            if smactinFlag.synthData(iMov) ~= 1
                
                tassm{iMov,iComb}{2} = [tassm{iMov,iComb}{2} ...
                    ...vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).smMeanDisp) ...     (:,4)
                    ...vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).smDensity), ...     (:,5)
                    ...vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).smComvmtAngle), ... (:,6)
                    vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).confRad),...        (:,7)
                    vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).trackClass)];   %(:,8)
            end
            
        case 10
            %concatenate nn: actin/ks, object: sm attributes
            tassm{iMov,iComb}{1} = [tassm{iMov,8}{1} tassm{iMov,7}{1}]; ... (:,1:3)
            tassm{iMov,iComb}{2} = tassm{iMov,8}{2};                       %(:,1:5)            
     
        end
        
        if iMatch == 1
            
            %during nn analysis, bisect nns that match or unmatch under
            %combineConditions. otherwise this should not be processed
            if iComb == 10
                
                matchindx{iMov,iComb} = matchindx{iMov,7} & matchindx{iMov,8};
            else
                %%% PUT BACK THE iMATCH %%%
                matchindx{iMov,iComb} = vertcat(combinePropPerTimeInt{iMov,iComb}(wind1st : end).dist_match) <= ...
                    combineConditions.matchRadius(iMov);
            end
        end
    end
end

end