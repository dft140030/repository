function [testResults,normAggMat,aggMat,trajClass] = smactinAggMov(tassm,matchindx,smactinFlag)


%% SMACTINAGGMOV
% Data analysis of multiple single-molecule and actin movies
%
% SYNOPSIS: [testResults,normAggMat,aggMat] = smactinAggMov(tassm,matchindx,smactinFlag)
%
% INPUT: tassm: data matrix of single-molecule and actin properties
%
%    matchindx: indices of speckles,ks, or sms matching to an object
%
% OUTPUT: testResults: cells containing structures with fields of different
%                      testing results on aggMat.
%
%                .pca: field containing principal components and percentage of total 
%                      variance explained by particular component. Analysis
%                      of properties for matched object/nn, unmatched
%                      object/nn, and averaged for matched object/neighbors
%
%              .mvrgs: field containing fields of general linear model ouptuts
%                      max iterations can be changed for the estimation
%                      algorithm used
%                   .coef: general linear model coefficients
%          .covarResponse: covariance matrix of response variables
%              .residuals: matrix of residuals for the fit
%             .covarParam: covariance matrix of regression coefficients
%
%               .tsne: field containing t-distributed stochastic neighbor
%                      embedding 2D points. Perplexity = 20 and the Learning
%                      Rate = 2000 can dramatically influence coordinates. 
%                      Use scatterplot to make use of the ML clustering algorithm
%
% Deryl Tschoerner, May 2018

%set up error message if number of single-molecule attributes is unequal in
%at least one movie
numMovs = size(tassm,1);
aggMat = cell(size(tassm,2),1);
normAggMat = aggMat;    

%loop through all possible combinations in movie channel
for iComb = 1 : length(unique(smactinFlag.combFlag))
    
    nnMat = [];
    obMat = [];
    trajClass = [];
    
    %for a combination, look through all movies
    for iMov = 1 : numMovs
        
        if smactinFlag.match(iMov) == 2 && ~isempty(tassm{iMov,iComb})
            
            matchindx{iMov,iComb} = 1 : size(tassm{iMov,iComb}{2},1);
        elseif isempty(tassm{iMov,iComb})
            
            continue
        end
        
        %for a particular combination, loop through all movies to get trajectory classes matched nn 
        trajClass = vertcat(trajClass,tassm{iMov,iComb}{2}(matchindx{iMov,iComb},end)); 
        
        if smactinFlag.synthData(iMov) == 1
            
            aggMat{iComb} = vertcat(aggMat{iComb}, horzcat(vertcat(nnMat,tassm{iMov,iComb}{1}(matchindx{iMov,iComb},:)), ...
                vertcat(obMat,tassm{iMov,iComb}{2}(matchindx{iMov,iComb},1 : end))));
            %{
            if
               THIS CAN BE USED IF ANALYSIS SHOULD BE PERFORMED ON UNMATCHED SINGLE MOLECULES FROM 
               NEAREST NEIGHBOR ANALYSIS
                %for a particular combination, loop through all movies to get unmatched nn properties 
                %umat{1} = vertcat(umat{1},tassm{iMov,iComb}{iMatch}{1}(~matchindx{iMov,iComb},:));

                %for a particular combination, loop through all movies to get object properties without matchings
                %umat{2} = vertcat(umat{2},tassm{iMov,iComb}{iMatch}{2}(~matchindx{iMov,iComb},1:end-1));

            %elseif iMatch == 2 && ~isempty(tassm{iMov,iComb})

                %for a particular combination, loop through all movies to get averaged neighbor properties
                %mat1{iMatch} = vertcat(mat1{iMatch},tassm{iMov,iComb}{iMatch}{1}(matchindx{iMov,iComb},:));

                %for a particular combination, loop through all movies to get object properties (with averaged neighbors)
                %mat2{iMatch} = vertcat(mat2{iMatch},tassm{iMov,iComb}{iMatch}{2}(matchindx{iMov,iComb},1:end-1));
            %end
            %}
        elseif smactinFlag.synthData(iMov) == 0
            
            %for a particular combination, loop through all movies to get matched nn/object properties
            aggMat{iComb} = vertcat(aggMat{iComb}, horzcat(vertcat(nnMat,tassm{iMov,iComb}{1} ...
                (matchindx{iMov,iComb},:)),vertcat(obMat,tassm{iMov,iComb}{2}(matchindx{iMov,iComb},1 : end - 1))));
        end
    end
    
    for iProp = 1 : size(aggMat{iComb},2)
        
        %normalized aggregate matrix of matched properties of nn combinations
        normAggMat{iComb} = horzcat(normAggMat{iComb},(aggMat{iComb}(:,iProp) - ...
            nanmean(aggMat{iComb}(:,iProp))) / nanstd(aggMat{iComb}(:,iProp)));       
    end
    
    numIndVars = size(tassm{1,iComb}{1},2);
    
    %remove the final column, as all receptors do not have a confinement
    %radius
    [testResults{iComb,1}.pca.coef,~,~,~,testResults{iComb,1}.pca.explain] = ...
        pca(normAggMat{iComb}(:,1:end - 1));
    
    oneVect = ones(size(aggMat{iComb},1),1);
    if smactinFlag.synthData(iMov) == 1

        %as before: removing the final column, as all receptors do not
        %have a confinement radius
        [testResults{iComb,1}.mvrgs.coef,testResults{iComb,1}.mvrgs.covarResponse, ...
            testResults{iComb,1}.mvrgs.residuals,testResults{iComb,1}.mvrgs.covarParam] = ...
            mvregress([oneVect aggMat{iComb}(:,1 : numIndVars)], ...
            aggMat{iComb}(:,numIndVars + 1 : end),'maxIter',100);
    
    elseif smactinFlag.synthData(iMov) == 0
        
        [testResults{iComb,1}.mvrgs.coef,testResults{iComb,1}.mvrgs.covarResponse, ...
            testResults{iComb,1}.mvrgs.residuals,testResults{iComb,1}.mvrgs.covarParam] = ...
            mvregress([oneVect normAggMat{iComb}(:,1 : numIndVars)], ...
            normAggMat{iComb}(:,numIndVars + 1: end - 1),'maxIter',100);
    else
        
        error('Invalid input for handling synthetic data.')
    end
    
    for iDiff = 0 : 3
        
        diffInd = trajClass == iDiff;
        oneVectDC = ones(sum(diffInd),1);
        if iDiff == 0 || iDiff == 1
            
            [testResults{iComb,2}.pca{iDiff + 1}.coef,~,~,~,testResults{iComb,2}.pca{iDiff + 1}.explain] = ...
                pca(normAggMat{iComb}(diffInd,:));
            
            [testResults{iComb,2}.mvrgs{iDiff + 1}.coef,testResults{iComb,2}.mvrgs{iDiff + 1}.covarResponse, ...
                testResults{iComb,2}.mvrgs{iDiff + 1}.residuals,testResults{iComb,2}.mvrgs{iDiff + 1}.covarParam] = ...
                mvregress([oneVectDC normAggMat{iComb}(diffInd,1 : numIndVars)], ...
                normAggMat{iComb}(diffInd,numIndVars + 1 : end));

        elseif iDiff == 2 || iDiff == 3

            [testResults{iComb,2}.pca{iDiff + 1}.coef,~,~,~,testResults{iComb,2}.pca{iDiff + 1}.explain] = ...
                pca(normAggMat{iComb}(diffInd,1 : end - 1));

            [testResults{iComb,2}.mvrgs{iDiff + 1}.coef,testResults{iComb,2}.mvrgs{iDiff + 1}.covarResponse, ...
                testResults{iComb,2}.mvrgs{iDiff + 1}.residuals,testResults{iComb,2}.mvrgs{iDiff + 1}.covarParam] = ...
                mvregress([oneVectDC normAggMat{iComb}(diffInd,1 : numIndVars)], ...
                normAggMat{iComb}(diffInd,numIndVars + 1 : end - 1));
        end
    end
end
testResults = testResults';

end