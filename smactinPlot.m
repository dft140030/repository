function [f] = smactinPlot(tassm,plotIndices,combineConditions,...
    figpath,foldername)

%% PLOTSPECKLESM 
%  plots individual joint speckle and single-molecule properties
%
%  SYNOPSIS: plotSpeckleSM(tassm,plotIndices,combineConditions,figpath, ...
%                foldername)
%
%  INPUT: tassm: (total attributes speckles and single molecules): data 
%                 matrix of single-molecule and actin properties
%          
%   plotIndices: indices of single-molecules meeting actin-specific
%                requirements
%
%   combineConditions: vector of structures with various fields imposing
%                      restrictions pertaining to both actin and SM.
%                      (see smactinPropPerTimeInt for more details)
%       
%       figpath: path to import figures created
%
%    foldername: name of folder for figure importation
%
%
%  OUTPUT: f: figures
%
% Deryl Tschoerner, March 2018

%% plotting
%axis titles. units element must be final element in the cell array.
%When adding new fields to plot, be sure to add below - as the
%scatterplot for loop is dependent upon the titles listed
set(0,'DefaultFigureVisible','off')
numMovs = length(tassm);
f = repmat(struct,numMovs,1);
cd(figpath)
mkdir(foldername)
suff = '.tif';

%plotAxis holds titles and units for sm, actin, and miscellaneous
%information, respectively.
plotAxis{1} = ...
    {{'SM Mean Displacement',' (units: pixels)'}, ...
     {'SM Net Speed',' (units: pixels/smInterval)'}, ...
     {'SM Mean Amplitude',' (units: ?)'}, ...
     {'SM Diffusion Coefficient',' (units: pixels^2 / frame)'}, ...
     {'SM Confinement Radius',' (units: pixels)'}};

plotAxis{2} = ...
    {{'Kinetic Score of Kinetic Event',' (units: arbitrary)'}, ...
    {'Speed of Actin Speckle Matching',' (units: pixels/frame)'}, ...
    {'Net Angle between SM velocity',' (units: radians)'}, ...
    {'Density of Speckles about a SM',' (units: speckles/ pixels^2)'}};
%Put Speckle Matching label in 2nd cell
plotAxis{3} = ...
    {'Kinetic Score Matching','','Magnitude of Polymerization', ...
    'Magnitude of Depolymerization'};

fn = {{'polymerization','depolymerization'};{'disp','speed','amp', ...
    'diff_coef','conf_rad'};{'kin_score','speed','angle','density'}};

for iMov = numMovs
    
    cd(foldername)
    mkdir(['movie_' num2str(iMov)])
    cd(['movie_' num2str(iMov)])
    mkdir(['match_radius_' num2str(combineConditions.matchRadius(iMov))])
    cd(['match_radius_' num2str(combineConditions.matchRadius(iMov))])
%{
    %row and column-wise number of panes for the subplot and initializing
    %number of bins for actin
    cpanes = 3;
    rpanes = length(plotAxis{2});
    xBins = nan(size(tracksAttributesActinSM{i,1}{1,1},2) - 1,19,2);
    yBins = nan(length(plotAxis{2}) + 1,19);
    yBins(3,:) = linspace(-4*10^-3,4*10^-3,19);
%}
    %% histogram of distances
    f(iMov).dist(1,1) = figure;
    histogram(tassm{iMov,1}{3,1}{3,1}(plotIndices{iMov,1}{5,1}(:,1),1), ...
        'Normalization','Probability')
    legend([num2str(sum(plotIndices{iMov,1}{5,1}(:,1))) 'NNs'])
    xlabel('Distance')
    ylabel('Normalized Count');
    saveas(f(iMov).dist(1,1),['dist_speckle_to_sm' suff])
    
    f(iMov).dist(2,1) = figure;
    histogram(tassm{iMov,1}{3,1}{4,1}(plotIndices{iMov,1}{7,1}(:,1),1), ...
        'Normalization','Probability')
    legend([num2str(sum(plotIndices{iMov,1}{7,1}(:,1))) 'NNs'])
    xlabel('Distance from speckle to NN sm')
    ylabel('Normalized Count');
    saveas(f(iMov).dist(2,1),['dist_sm_to_speckle' suff])
    
    %% histograms of polymerization and depolymerization
    sgn = 1; 
    for j = 1 : 2
%{
        ysubplot(2,1,k)NormThresh = yCounts / sum(yCounts);
        stairs(yBins(k,:),[yNormThresh, 0])
%}
        f(iMov).poly(j) = figure;
        polymId = tassm{iMov,1}{2,1}(:,1) * sgn >= 0;
%%%%%%%%%%%%%%%%%%%%%%% MAKE BINS NICER %%%%%%%%%%%%%%%%%%%%%%%%%%
        histogram(tassm{iMov,1}{2,1}(polymId & ...
            plotIndices{iMov,1}{1,1}(:,1),1) * sgn,'Normalization','Probability')
        legend([num2str(sum(polymId)) ' kinetic scores'])
        xlabel(plotAxis{3}{j + 2})
        ylabel('Count');
        saveas(f(iMov).poly(j),[fn{1}{j} suff])
        sgn = -sgn;
    end

    %% figures for every single molecule attribute
    for j = 1 : length(plotAxis{1})      

        %% histogram of matchings of single molecule attribute and
        %% 2-sample Kolmogorov-Smirnov test for difference in pmfs
        xBins = nan(length(plotAxis{1}),19);
        xBins(j,:) = linspace(min(tassm{iMov,1}{1,1}(:,j)), ...
            max(tassm{iMov,1}{1,1}(:,j)),19);
        for k = 1 : 2
            
            f(iMov).sm(j,k) = figure;
            histogram(tassm{iMov,1}{1,1} ...
                (plotIndices{iMov,1}{2*k - 1,1}(:,1),j),xBins(j,:), ...
                'Normalization','Probability')
            hold on
            %change ybins: adjusted for normalization
            histogram(tassm{iMov,1}{1,1}(plotIndices{iMov,1}{2*k}(:,1),j), ...
                xBins(j,:),'Normalization','Probability')
            title(plotAxis{3}(k))
            legend([num2str(sum(plotIndices{iMov,1}{2*k - 1,1}(:,1))) ' Matched SM'], ...
                [num2str(sum(plotIndices{iMov,1}{2*k,1}(:,1))), ' Unmatched SM'])
            xlabel(plotAxis{1}{j}(1,:))
            ylabel(['Normalized' 'Counts']);
            saveas(f(iMov).sm(j,1),[fn{2}{j} suff])
%{
            %histcounts returns vector of counts with (not) met threshold
            %values to be normalized
            xCounts = [histcounts(tracksAttributesActinSM{i,1}{1,1} ...
                (plotIndices{2*k - 1},j), xBins(j,:,1)); ...
                histcounts(tracksAttributesActinSM{i,1}{1,1} ...
                (plotIndices{2*k},j), xBins(j,:,1))];

            %vector of normalized counts with (not) met threshold values to
            %be used as bins
            xNorm = [xCounts(1,:) / sum(xCounts(1,:)); ...
                xCounts(2,:) / sum(xCounts(2,:))];

            %histogramming for (non-)thresholded sm property ... [1, 1 : 2]
            %stairs(xBins(j,:,1),[xNorm(1,:), 0])
            hold on
            %stairs(xBins(j,:,1),[xNorm(2,:), 0])
%}
        end
        colswp = 1;
        f(iMov).sm_speckle(j,1) = figure;
        for k = 1 : length(plotAxis{2})
%{
            yBins(k + 3,:) = linspace(min(tracksAttributesActinSM{i,1}{2,1} ...
                (plotIndices{1},2*k + 1)), max(tracksAttributesActinSM{i,1}{2,1} ...
                (plotIndices{1},2*k + 1)),19);
            yCounts = histcounts(tracksAttributesActinSM{i,1}{2,1} ...
                (plotIndices{1},2*k + 1), yBins(k + 3,:));
            yNorm = yCounts / sum(yCounts);
%}                                                
            %% histograms of speckle attributes vs. single molecule attribute
            yBins(k,:) = linspace(min(tassm{iMov,1}{2,1}(:,2*k - 1)), ...
                max(tassm{iMov,1}{2,1}(:,2*k - 1)),19);

            f(iMov).sm_speckle(j,k) = figure;
            histogram2(tassm{iMov,1}{1,1}(plotIndices{iMov,1}{colswp,1}(:,1),j), ...
                tassm{iMov,1}{2,1}(plotIndices{iMov,1}{colswp,1}(:,1),2*k - 1), ...
                xBins(j,:),yBins(k,:),'DisplayStyle','tile','ShowEmptyBins','on');
            colorbar;
            xlabel(plotAxis{1}{j}(1,:));
            ylabel(plotAxis{2}{k}(1,:));
            saveas(f(iMov).sm_speckle(j,k),[fn{2}{j} '_vs._' fn{3}{k} suff])

            testResults{iMov,1}.covar{j,k} = ...
                cov(tassm{iMov,1}{1,1}(plotIndices{iMov,1}{colswp,1}(:,1),j), ...
                tassm{iMov,1}{2,1}(plotIndices{iMov,1}{colswp,1}(:,1),2*k - 1));

            %% histogram for every speckle attribute (only want to loop once)
            if j == 1

                for l = 1 : 2

                    f(iMov).speckle(k,l) = figure;
                    histogram(tassm{iMov,1}{2,1} ...
                        (plotIndices{iMov,1}{colswp,1}(:,1),2*k - 1), ...
                        yBins(k,:), 'Normalization','Probability')
                    hold on
%%%%%%%%%%% MAKE SURE BINS ARE OF EQUAL WIDTH %%%%%%%%%%%%%%%
                    histogram(tassm{iMov,1}{2,1} ...
                        (plotIndices{iMov,1}{colswp + 1,1}(:,1),2*k - 1), ...
                        yBins(k,:),'Normalization','Probability')

                    title([plotAxis{3}(l), ' vs. ', plotAxis{2}{k}(1: end - 1)])
                    legend([num2str(sum(plotIndices{iMov,1}{2*k - 1,1}(:,1))) ...
                        ' matched with SM'],[num2str(sum(plotIndices ...
                        {iMov,1}{2*k,1}(:,1))), ' not matched with SM'])
                    xlabel(plotAxis{2}{k}(1,:))
                    ylabel('Counts');
                    saveas(f(iMov).speckle(k,l),[fn{3}{k} suff])
                end
                colswp = 3;
            end
        end
    end
end

end