% Code for the manuscript Clark et al, 2024 for crowding.
%   Each row represents a single trial.
%   "Condition" is -1 for uncrowded and +1 for crowded conditions.
%   "Offset" is the average euclidean distance from the center
%       of the target (0,0);
%   "DC" is the diffusion constant.
%   "Curv" is the curvature
%   "Size" is the stimulus width in arcminutes
%   "Span is the average span of the eye trace
%   "SameConeProb" is the individual trial probabilities of the
%       same cones being stimulated by both a target and a
%       flanker.
%   "Area" = is the area of the trace in x and y for .68
%   "Dtrend" = the D value but removing the bias in direction
%   "Bias" = estimate for directional offset
%   "prl" = location on the monitor subjects spend the most time with gaze
%
% TO DO - add code for D, Area, 
%       - make public version
% Code written by Ashley M. Clark, 2024 Active Perception Lab

load('dataAllSingleSize.mat');
temp = load('dataAllSizes.mat');
dataAll2 = temp.dataAll2;
%%% for default figs pick a trial you want to check analysis on
trial = 1889;
stimulusSpacingExperiment = 1.4; %all stimuli were spaced 1.4x C-C

%%

allProbsUnique = unique(dataAll.SameConeProb);
for u = 1:length(allProbsUnique)
    tempIdx = [];
    tempIdx = find(dataAll.SameConeProb == allProbsUnique(u));
    meanEachProb(u) = mean(dataAll.Area(tempIdx));
end
% mdlAll = fitlm((newDataAll.SameConeProb,(newDataAll.AreaP))
tempVal = [allProbsUnique'];
mdlAll = fitlm(tempVal,meanEachProb);

figure;
% hold on
for ii = 1:length(unique([dataAll.SubjectID(:)]))
    subplot(4,4,ii)
    scatter1 = scatter([dataAll.SameConeProb(find(dataAll.Condition == 1 & ...
        dataAll.SubjectID == ii))],...
        [(dataAll.Area(find(dataAll.Condition == 1 & ...
        dataAll.SubjectID == ii)))],...
        'MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor',[.8 .8 .8]);
    scatter1.MarkerFaceAlpha = 0;
    scatter1.MarkerEdgeAlpha = 0;
    if ii == 10
        mdlEach{ii} = fitlm([dataAll.SameConeProb(find(dataAll.Condition == 1 & ...
            dataAll.SubjectID == ii & dataAll.Area < 80))],...
            ([(dataAll.Area(find(dataAll.Condition == 1 & ...
            dataAll.SubjectID == ii & dataAll.Area < 80)))]));
    else
        mdlEach{ii} = fitlm([dataAll.SameConeProb(find(dataAll.Condition == 1 & ...
            dataAll.SubjectID == ii))],([(dataAll.Area(find(dataAll.Condition == 1 & ...
            dataAll.SubjectID == ii)))]));
    end
    hold on
    plot(mdlEach{ii});
    title(sprintf('p = %.2f, R-squared = %.2f',...
        mdlEach{ii}.Coefficients.pValue(2),...
        mdlEach{ii}.Rsquared.Ordinary));
    if ii == 1
        xlabel('Probability')
        ylabel('Ocular Drift Area (arcmin^2)')
    else
        xlabel('');
        ylabel('');
    end
    title(sprintf('S%i',ii))
        text(.10,125,'p < 0.001');
%     else
%         text(.10,125,sprintf('p > %.2f',mdlEach{ii}.Coefficients.pValue(2)));
%     end
    text(.10,110,sprintf('R-squared = %.2f',mdlEach{ii}.Rsquared.Ordinary));
    legend off
    xlim ([0 .8])
    ylim([0 150])
    axis square
end


%% Figure 1B
figure;
plot(dataAll.X{trial},dataAll.Y{trial},'-','Color','k');

%% Figure 1C
axisForCones = 15; %in arcminutes
yholder = ([-axisForCones:.5:axisForCones]); %setting the vernoi spacing to be a 1/2 arcmin
counter = 1;
for i = 1:(length(yholder))
    if counter == 1
        xV(i,:) = ones(1,length(yholder))*(yholder(i));
        yV(i,:) = yholder;
        counter = counter - 1;
    else
        xV(i,:) = (ones(1,length(yholder))*(yholder(i)));
        yV(i,:) = yholder-.25;
        counter = counter + 1;
    end
end

figure;
subplot(1,2,1)
voronoiAMC(xV,yV);
axis square
axis([-axisForCones axisForCones -axisForCones axisForCones])
% [Xe, Ye] = EM_Brownian2(30, Fs, nSamples, 1);
Xe = dataAll.X{trial};
Ye = dataAll.Y{trial};
hold on
AlphaVal = 1;
for i = 1
    rectangle('Position',[Xe(i)-2.5,Ye(i)-2.5,5,5],...
        'EdgeColor',[0 0 1 AlphaVal],'FaceColor', [0 0 1 AlphaVal])
    hold on
end
factor = 5*1.4;
for i = 1
    rectangle('Position',[(Xe(i)-2.5)+factor,Ye(i)-2.5,5,5],...
        'EdgeColor',[1 0 0 AlphaVal],'FaceColor', [1 0 0 AlphaVal])
    hold on
end
for i = 1
    rectangle('Position',[(Xe(i)-2.5)-factor,Ye(i)-2.5,5,5],...
        'EdgeColor',[1 0 0 AlphaVal],'FaceColor', [1 0 0 AlphaVal])
    hold on
end

ylabel('Y Position (arcmin')
xlabel('X Position (arcmin')

subplot(1,2,2)
voronoiAMC(xV,yV);
axis square
axis([-axisForCones axisForCones -axisForCones axisForCones])
hold on
AlphaVal = .01;
for i = 1:length(Xe)
    rectangle('Position',[Xe(i)-2.5,Ye(i)-2.5,5,5],...
        'EdgeColor',[0 0 1 AlphaVal],'FaceColor', [0 0 1 AlphaVal])
    hold on
end
factor = 5*1.4;
for i = 1:length(Xe)
    rectangle('Position',[(Xe(i)-2.5)+factor,Ye(i)-2.5,5,5],...
        'EdgeColor',[1 0 0 AlphaVal],'FaceColor', [1 0 0 AlphaVal])
    hold on
end
for i = 1:length(Xe)
    rectangle('Position',[(Xe(i)-2.5)-factor,Ye(i)-2.5,5,5],...
        'EdgeColor',[1 0 0 AlphaVal],'FaceColor', [1 0 0 AlphaVal])
    hold on
end

ylabel('Y Position (arcmin')
xlabel('X Position (arcmin')

%% Figure 1D
coneImage = Tiff('Z190_R_2023_07_25_coneTagImg.TIF');
imageData = read(coneImage);

% get cone density map, subject 9
tempData = load('Z190_R_2023_07_25_CD_data.mat');
centerImage = [tempData.PRL_X(1),tempData.PRL_Y(1)];

pixPerArc = size(imageData)/60;
pixPerArc = round(mean(pixPerArc));
rowsInclude = round(cell2mat(centerImage(1)))-(pixPerArc*10):...
    round(cell2mat(centerImage(1)))+(pixPerArc*10);
columnsInclude =  round(cell2mat(centerImage(2)))-(pixPerArc*10):...
    round(cell2mat(centerImage(2)))+(pixPerArc*10);
figure;
imshow(imageData(rowsInclude,columnsInclude)); %image is 1 degree width and height, crop so its pm 10

% get heatmap of eye movements, subject 9
figure;
idx = find([dataAll.SubjectID(:)] == 9 & ...
    round([dataAll.Size(:)],3) == 1.1230 & ... %threshold
    [dataAll.Condition(:)] == 1);

x = [];
y = [];

threshold = 2.1129; %threshold size for this subject

xAll = [dataAll.X{idx}];
yAll = [dataAll.Y{idx}];
stVector = [-2:1:2]; %width of square stim in arcminutes
for sty = 1:length(stVector)
    for stx = 1:length(stVector)
        xT = xAll(ismember(round(xAll),...
            [round(xAll+(threshold*stimulusSpacingExperiment)) ...
            round(xAll-(threshold*stimulusSpacingExperiment))]));
        yT = yAll(ismember(round(yAll),...
            [round(yAll+(threshold*stimulusSpacingExperiment)) ...
            round(yAll-(threshold*stimulusSpacingExperiment))]));
        
        x = [x (xT + stVector(stx))];
        y = [y (yT + stVector(sty))];
    end
end


figure;
temp = generateHeatMapSimple( ...
    x,... %flip along the vertical axis to match map
    -y, ...
    'Bins', 40,... %for a 20 window width, each bin is .5arcmin
    'StimulusSize', threshold,...
    'AxisValue', 10,... %plus or minus
    'Uncrowded', 0,...
    'Borders', 1);
axis square;

%% Figure 3 A&B

%load(dataAll2.mat); %larger set for every trial - not just at threshold
xl = [0.1, 11];
yl = [0, 1.05];
gamma = 0.25;

% load(thresh.mat);

for ii = 1:length(unique([dataAll2.SubjectID(:)]))
   for c = 1:2
       if c == 1
           cond = -1; %uncrowded
       else
           cond = 1; %crowded
       end
       idx = find([dataAll2.SubjectID(:)] == ii & ...
           [dataAll2.Perf(:)] < 3 & ...
           [dataAll2.Condition(:)] == cond);

       [thresh(ii,c),~,~,xi{ii,c},yi{ii,c}] = psyfitCrowding(...
          round([dataAll2.Size(idx)]*2,1), ... %rounding to reflect monitor resolution
          [dataAll2.Perf(idx)], 'DistType', 'Normal',...
           'Xlim', xl, 'Ylim', yl,...
           'Chance', gamma, 'Extra');

   end
end


%% do quick check to see what the spacing thresh is and what the spacing needed to have no overlap is

for ii = 1:length(unique([dataAll2.SubjectID(:)]))
    
    idx = find(yi{ii,2} < .95 & yi{ii,2} > .9);
    avSpacing(ii) = mean(xi{ii,2}(idx))*1.4;
    stimSize = mean(xi{ii,2}(idx));
    %what spacing and size is needed to reach plateu probability for this
    %dc
    idxTrialsDC = find(dataAll.Condition(:) == 1 &...
        dataAll.SubjectID(:) == ii);
    dc = mean(find(dataAll.DC(idxTrialsDC)));

    counter = 1;
    fc = 1000;
    N = 500; %(500ms)
    DC = dc;
    [Xe, Ye] = EM_Brownian(DC,  fc, N, 1000);
%     stimSize = .5:.5:5;
    spacing = 1:.15:20;
    for st = 1%:length(stimSize)
        for sp = 1:length(spacing)
            probSameCone(st,sp) = coneAndFlankerProbability_OG([Xe(:)'],stimSize,spacing(sp));
            %         spac(counter) = spacing;
            %         siz(counter) = stimSize;
            %         counter = counter + 1;
        end
    end
    
    [col,row]=find(probSameCone >.01 & probSameCone < .1)  % corresponding indices

    spacingTested(ii) = nanmean(spacing(row'));
    valUsed(ii) = value;

end

mdl_Sp = fitlm(spacingTested,avSpacing);
figure;
plot(mdl_Sp)
xlabel('Predicted Spacing Needed from D');
ylabel('Spacing at 92.5\% Performance');
title(sprintf('p = %.2f',table2array(mdl_Sp.Coefficients(2,4))));

%% Figure 3C
figure;
for ii = 1:length(unique([dataAll.SubjectID(:)]))
    plot([1 2],[thresh(ii,1),thresh(ii,2)],'-o');
    hold on
end
xlim([.5 2.5])
xticks([1 2])
xticklabels({'Uncrowded','Crowded'});

%% Figure 3D
figure;
for ii = 1:length(unique([dataAll.SubjectID(:)]))
    clear perf
    idx = (find([dataAll2.SubjectID(:)] == ii & ...
        [dataAll2.Perf(:)] < 3 & ...
        [dataAll2.Condition(:)] == -1));
    allSizesUncrowded = unique(dataAll2.Size(idx));
    for i = 1:length(allSizesUncrowded)
        idxPerf = find([dataAll2.SubjectID(:)] == ii & ...
            [dataAll2.Perf(:)] < 3 & ...
            [dataAll2.Condition(:)] == -1 &...
            [dataAll2.Size(:)] == allSizesUncrowded(i));
        perf(i,1) = mean([dataAll2.Perf(idxPerf)]);
    end
    if ii == 8 || ii == 4 %specifiying when there are not enough trials in one of the conditions
        [~,chooseSizeIdx] = min(abs(perf - 0.625));
    elseif  ii == 6
        [~,chooseSizeIdx] = min(abs(perf - 0.85));
    else
        [~,chooseSizeIdx] = min(abs(perf - 0.7));
    end
    
    idx2 = find([dataAll2.SubjectID(:)] == ii & ...
        [dataAll2.Perf(:)] < 3 & ...
        [dataAll2.Condition(:)] == 1);
    
    allSizesCrowded = unique([dataAll2.Size(idx2)]);
    [~,chooseSizeIdxC] = min(abs(allSizesCrowded - allSizesUncrowded(chooseSizeIdx)));
    
    idx3 = find([dataAll2.SubjectID(:)] == ii & ...
        [dataAll2.Perf(:)] < 3 & ...
        [dataAll2.Condition(:)] == 1 & ...
        [dataAll2.Size(:)] == allSizesCrowded(chooseSizeIdxC));
    
    crowdedSize(ii,1) = (allSizesUncrowded(chooseSizeIdx));
    crowdedSize(ii,2) = (allSizesCrowded(chooseSizeIdxC));
    
    plot([1 2], [perf(chooseSizeIdx) mean([dataAll2.Perf(idx3)])],...
        '-o');
    
    perfCom(ii,1) = perf(chooseSizeIdx);
    perfCom(ii,2) = mean([dataAll2.Perf(idx3)]);
    hold on
end
xlim([.5 2.5])
ylim([0 1])
xticks([1 2])
xticklabels({'Uncrowded','Crowded'});

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Figure 4A
figure;
subplot(1,2,1)

idxSmall = find([dataAll.Area(:)] < 7);
idxLarge = find([dataAll.Area(:)] > 20);

temp = generateHeatMapSimple( ...
    [dataAll.X{idxSmall}],... %flip along the vertical axis to match map
    [dataAll.Y{idxSmall}], ...
    'Bins', 40,... %for a 20 window width, each bin is .5arcmin
    'StimulusSize', thresh(1),...
    'AxisValue', 20,... %plus or minus
    'Uncrowded', 0,...
    'Borders', 1);
axis square;

subplot(1,2,2)
temp = generateHeatMapSimple( ...
    [dataAll.X{idxLarge}],... %flip along the vertical axis to match map
    [dataAll.Y{idxLarge}], ...
    'Bins', 40,... %for a 20 window width, each bin is .5arcmin
    'StimulusSize', thresh(2),...
    'AxisValue', 20,... %plus or minus
    'Uncrowded', 0,...
    'Borders', 1);
axis square;

%% Figure 4B

for ii = 1:length(unique([dataAll.SubjectID(:)]))
    idx = [];
    idx = find([dataAll.SubjectID(:)] == ii & ...
            [dataAll.Perf(:)] < 3 & ...
            [dataAll.Condition(:)] == -1);
    areaPolyA{1,ii} = [dataAll.Area(idx)];
    meanPerf(1,ii) = mean([dataAll.Perf(idx)]);
    idx = [];
    idx = find([dataAll.SubjectID(:)] == ii & ...
            [dataAll.Perf(:)] < 3 & ...
            [dataAll.Condition(:)] == 1);
    areaPolyA{2,ii} = [dataAll.Area(idx)];
    meanPerf(2,ii) = mean([dataAll.Perf(idx)]);
end



smallD = [];
largeD = [];
indSubBinUpper = ones(1,length(unique([dataAll.SubjectID(:)])))*30;
indSubBinLower = ones(1,length(unique([dataAll.SubjectID(:)])))*70;
[smallD,largeD] = calculatePerfProbAtThresh(indSubBinUpper,indSubBinLower,unique([dataAll.SubjectID(:)]),...
    dataAll,areaPolyA,'Area',1,smallD,largeD);

indSubBinUpper(1) = 70;
indSubBinLower(1) = 20;
indSubBinUpper(2) = 75;
indSubBinLower(2) = 20;
indSubBinUpper(3) = 70;
indSubBinLower(3) = 35;
indSubBinUpper(4) = 70;
indSubBinLower(4) = 20;
indSubBinUpper(5) = 70; 
indSubBinLower(5) = 40;
indSubBinUpper(6) = 60;
indSubBinLower(6) = 60;
indSubBinUpper(7) = 85;
indSubBinLower(7) = 41; 
indSubBinUpper(8) = 60;
indSubBinLower(8) = 50;
indSubBinUpper(9) = 70;
indSubBinLower(9) = 20;
indSubBinUpper(10) = 50;
indSubBinLower(10) = 50;
indSubBinUpper(11) = 45;
indSubBinLower(11) = 45;
indSubBinUpper(12) = 50;
indSubBinLower(12) = 30;
indSubBinUpper(13) = 70;
indSubBinLower(13) = 30;

[smallD,largeD] = calculatePerfProbAtThresh(indSubBinUpper,indSubBinLower,unique([dataAll.SubjectID(:)]),...
    dataAll,areaPolyA,'Area',2,smallD,largeD);

[~,pU] = ttest(smallD.Unc.Performance,largeD.Unc.Performance);
[~,pC] = ttest(smallD.Cro.Performance,largeD.Cro.Performance);

figure;
makePlotPerfDiffUnityLine([smallD.Unc.Performance;smallD.Cro.Performance],...
    [largeD.Unc.Performance; largeD.Cro.Performance], pU, pC);
suptitle('Area')
saveas(gca,'AreaUnityLine.epsc');


mean([smallD.Unc.Performance largeD.Unc.Performance])
mean([smallD.Cro.Performance largeD.Cro.Performance])
std([smallD.Unc.Performance largeD.Unc.Performance])
std([smallD.Cro.Performance largeD.Cro.Performance])

[h,p,ci,stats] = ttest([smallD.Unc.Performance largeD.Unc.Performance],...
    [smallD.Cro.Performance largeD.Cro.Performance])

%% Fit the different thresholds based on small and large AREA

xl = [0.1, 11];
yl = [0, 1.05];
gamma = 0.25;


for ii = 1:length(unique([dataAll2.SubjectID(:)]))
    for c = 1:2
        if c == 1
            cond = -1; %uncrowded
        else
            cond = 1; %crowded
        end
        idx = find([dataAll2.SubjectID(:)] == ii & ...
            [dataAll2.Perf(:)] < 3 & ...
            [dataAll2.Condition(:)] == cond);
        
        meanArea = median(dataAll2.Area(idx));
        
        idxS = find([dataAll2.SubjectID(:)] == ii & ...
            [dataAll2.Perf(:)] < 3 & ...
            [dataAll2.Condition(:)] == cond & ...
            [dataAll2.Area(:)] < meanArea - (std(dataAll2.Area(idx))/3));
        
        [smallAreaThresh(ii,c)] = psyfitCrowding(...
            round([dataAll2.Size(idxS)]*2,1), ... %rounding to reflect monitor resolution
            [dataAll2.Perf(idxS)], 'DistType', 'Normal',...
            'Xlim', xl, 'Ylim', yl,...
            'Chance', gamma, 'Extra');
        
        idxL = find([dataAll2.SubjectID(:)] == ii & ...
            [dataAll2.Perf(:)] < 3 & ...
            [dataAll2.Condition(:)] == cond & ...
            [dataAll2.Area(:)] > meanArea + (std(dataAll2.Area(idx))/3));
        
        [largeAreaThresh(ii,c)] = psyfitCrowding(...
            round([dataAll2.Size(idxL)]*2,1), ... %rounding to reflect monitor resolution
            [dataAll2.Perf(idxL)], 'DistType', 'Normal',...
            'Xlim', xl, 'Ylim', yl,...
            'Chance', gamma, 'Extra');
        
    end
end


 figure;plot([1 2],[smallAreaThresh(:,2) largeAreaThresh(:,2)],'-o')
[h,p] = ttest([smallAreaThresh(:,2)'], [largeAreaThresh(:,2)'])
xlim([0 3])
xticks([1 2])
xticklabels('Small Area','Larger Area')
xticklabels({'Small Area','Larger Area'})
ylabel('Threshold')
ylim([1 4])
title(sprintf('p = %.2f',p))

%% Theoretical Figures: 
%   3 different DC and different sizes and different spacings (color is
%   probability of overlap where spacing is on x axis and y axis speed, the
%   color of the bin is the probability o shared stimulation over the
%   stimulated cones color is the probability


% idx = find(dataAll2.DC < 100 & dataAll2.Condition ==1)
% figure;
% plot(dataAll2.DC(idx),dataAll2.Size(idx),'.');
dc = [5 15 30 60];
for d = 1:4
    counter = 1;
    fc = 1000;
    N = 500; %(500ms)
    DC = dc(d);
    [Xe, Ye] = EM_Brownian(DC,  fc, N, 1000);
    stimSize = .5:.5:5;
    spacing = 0:.5:4;
    for st = 1:length(stimSize)
        for sp = 1:length(spacing)
            probSameCone(st,sp) = coneAndFlankerProbability_OG([Xe(:)'],stimSize(st),spacing(sp));
            %         spac(counter) = spacing;
            %         siz(counter) = stimSize;
            %         counter = counter + 1;
        end
    end
    
    subplot(2,2,d)
    image(probSameCone,'CDataMapping','scaled')
    yticklabels(stimSize)
    xticks(1:length(spacing))
    xticklabels(spacing)
    cb = colorbar();
    ylabel(cb,'Probability of Shared Cone Stimulation')
    ylabel('Stimulus Width (arcmin)');
    xlabel('Stimulus Spacing (arcmin)');
    title(sprintf('D = %.2f arcmin^2/sec',DC));
end
%% Figure 4C
figure;

errorbar([1 2],([mean(smallD.Unc.Performance)-mean(meanPerf(1,:)),...
    mean(largeD.Unc.Performance)-mean(meanPerf(1,:))]),...
    [sem(smallD.Unc.Performance),...
    sem(largeD.Unc.Performance)],'-d','Color','b',...
    'MarkerFaceColor','b','MarkerSize',12);
hold on
errorbar([1 2],([mean(smallD.Cro.Performance)-mean(meanPerf(2,:)),...
    mean(largeD.Cro.Performance)-mean(meanPerf(2,:))]),[sem(smallD.Cro.Performance),...
    sem(largeD.Cro.Performance)],'-d','Color','r',...
    'MarkerFaceColor','r','MarkerSize',12);
xlim([0.5 2.5])
xticks([1 2])
xticklabels({'Small','Large'})
xlabel('Fixation Area')
ylabel('Change in Performance from Threshold')
line([0.5 2.5],[0 0],'LineStyle','--');
% yticks([.8 .9 1])
% ylim([.78 1.05])
saveas(gca,'AreaThreshChange.epsc');


%% Get Cone Probability for Trial 1

xTemp = dataAll.X{1,1};
stimSize = dataAll.Size(1);

probSameCone(1) = coneAndFlankerProbability_OG(xTemp,stimSize,1.4);

%% Figure 4D

for ii = 1:length(unique([dataAll.SubjectID(:)]))
    idxTrials = find(dataAll.SubjectID == ii & ...
        dataAll.Condition == 1 & ...
        dataAll.Area < 200);
    halfwayPoint = mean(dataAll.Area(idxTrials));
    halfwayPointSTD = std(dataAll.Area(idxTrials));

    perfSmall(ii) = mean(dataAll.SameConeProb(find(dataAll.SubjectID == ii & ...
        dataAll.Condition == 1 & ...
        dataAll.Area < 200 & ...
        dataAll.Area < halfwayPoint)));
    perfLarge(ii) = mean(dataAll.SameConeProb(find(dataAll.SubjectID == ii & ...
        dataAll.Condition == 1 & ...
        dataAll.Area < 200 & ...
        dataAll.Area > halfwayPoint)));
end

figure;
errorbar([1 2],[mean(perfSmall),...
    mean(perfLarge)],[sem(perfSmall),...
    sem(perfLarge)],'-d','Color','g',...
    'MarkerFaceColor','g','MarkerSize',12);
xlim([0.5 2.5])
xticks([1 2])
xticklabels({'Small','Large'})
xlabel('Fixation Area')
ylabel({'Probability of Cone Being Stimulated','by both Target and Flanker'})
[~,pC] = ttest(perfSmall,perfLarge);
saveas(gca,'ConeProbThreshChange.epsc');

%% Figure 5B (7 subjects)
clear thresh
load('stabThresh.mat');

figure;
for ii = 1:length(thresh)
    forLegs(ii) = plot([1.1 2.2],[thresh{1,ii},...
        thresh{2,ii}],'-','Color',[.7 .7 .7]);
    hold on
end
threshMat = cell2mat(thresh);

errorbar([1 2],[mean(threshMat(1,:)) mean(threshMat(2,:))],...
    [sem(threshMat(1,:)) sem(threshMat(2,:))],'-o',...
    'MarkerFaceColor','k','Color','k');
xticks([1 2])
xlim([0.5 2.5])
xticklabels({'Stabilized','Unstabilized'});
ylabel('Nominal Critical Spacing')
[h,p] = ttest(threshMat(1,:),threshMat(2,:));
title(sprintf('p = %.2f',p));
ylim([0.8 2.2]);

%% Looking at DC instead of Area
for ii = 1:length(unique([dataAll.SubjectID(:)]))
    idx = [];
    idx = find([dataAll.SubjectID(:)] == ii & ...
            [dataAll.Perf(:)] < 3 & ...
            [dataAll.Condition(:)] == -1);
    areaPolyD{1,ii} = [dataAll.DC(idx)];
    meanPerf(1,ii) = mean([dataAll.Perf(idx)]);
    idx = [];
    idx = find([dataAll.SubjectID(:)] == ii & ...
            [dataAll.Perf(:)] < 3 & ...
            [dataAll.Condition(:)] == 1);
    areaPolyD{2,ii} = [dataAll.DC(idx)];
    meanPerf(2,ii) = mean([dataAll.Perf(idx)]);
end



smallD = [];
largeD = [];
indSubBinUpper = ones(1,length(unique([dataAll.SubjectID(:)])))*75;
indSubBinLower = ones(1,length(unique([dataAll.SubjectID(:)])))*20;
indSubBinUpper(7) = 85;
indSubBinUpper(13) = 30;
indSubBinLower(13) = 20;
indSubBinUpper(6) = 95;
indSubBinLower(6) = 60;
indSubBinUpper(9) = 60;
indSubBinLower(9) = 30;
indSubBinUpper(1) = 70;
indSubBinLower(1) = 50;

[smallD,largeD] = calculatePerfProbAtThresh(indSubBinUpper,indSubBinLower,unique([dataAll.SubjectID(:)]),...
    dataAll,areaPolyD,'DC',1,smallD,largeD);

indSubBinUpper = ones(1,length(unique([dataAll.SubjectID(:)])))*75;
indSubBinLower = ones(1,length(unique([dataAll.SubjectID(:)])))*30;

[smallD,largeD] = calculatePerfProbAtThresh(indSubBinUpper,indSubBinLower,unique([dataAll.SubjectID(:)]),...
    dataAll,areaPolyD,'DC',2,smallD,largeD);

[~,pU] = ttest(smallD.Unc.Performance,largeD.Unc.Performance);
[~,pC] = ttest(smallD.Cro.Performance,largeD.Cro.Performance);



figure;
makePlotPerfDiffUnityLine([smallD.Unc.Performance;smallD.Cro.Performance],...
    [largeD.Unc.Performance; largeD.Cro.Performance], pU, pC);
suptitle('DC')

figure;

errorbar([1 2],([mean(smallD.Unc.Performance)-mean(meanPerf(1,:)),...
    mean(largeD.Unc.Performance)-mean(meanPerf(1,:))]),...
    [sem(smallD.Unc.Performance),...
    sem(largeD.Unc.Performance)],'-d','Color','b',...
    'MarkerFaceColor','b','MarkerSize',12);
hold on
errorbar([1 2],([mean(smallD.Cro.Performance)-mean(meanPerf(2,:)),...
    mean(largeD.Cro.Performance)-mean(meanPerf(2,:))]),[sem(smallD.Cro.Performance),...
    sem(largeD.Cro.Performance)],'-d','Color','r',...
    'MarkerFaceColor','r','MarkerSize',12);
xlim([0.5 2.5])
xticks([1 2])
xticklabels({'Small','Large'})
xlabel('Fixation DC')
ylabel('Change in Performance from Threshold')
line([0.5 2.5],[0 0],'LineStyle','--');
% yticks([.8 .9 1])
% ylim([.78 1.05])

%% PRL
for ii = 1:length(unique([dataAll.SubjectID(:)]))
    idx = [];
    idx = find([dataAll.SubjectID(:)] == ii & ...
            [dataAll.Perf(:)] < 3 & ...
            [dataAll.Condition(:)] == -1);
    areaPolyD{1,ii} = [dataAll.Offset(idx)];
    meanPerf(1,ii) = mean([dataAll.Perf(idx)]);
    idx = [];
    idx = find([dataAll.SubjectID(:)] == ii & ...
            [dataAll.Perf(:)] < 3 & ...
            [dataAll.Condition(:)] == 1);
    areaPolyD{2,ii} = [dataAll.Offset(idx)];
    meanPerf(2,ii) = mean([dataAll.Perf(idx)]);
end



smallD = [];
largeD = [];
indSubBinUpper = ones(1,length(unique([dataAll.SubjectID(:)])))*75;
indSubBinLower = ones(1,length(unique([dataAll.SubjectID(:)])))*30;

[smallD,largeD] = calculatePerfProbAtThresh(indSubBinUpper,indSubBinLower,unique([dataAll.SubjectID(:)]),...
    dataAll,areaPolyD,'Offset',1,smallD,largeD);

indSubBinUpper = ones(1,length(unique([dataAll.SubjectID(:)])))*75;
indSubBinLower = ones(1,length(unique([dataAll.SubjectID(:)])))*30;

[smallD,largeD] = calculatePerfProbAtThresh(indSubBinUpper,indSubBinLower,unique([dataAll.SubjectID(:)]),...
    dataAll,areaPolyD,'Offset',2,smallD,largeD);

[~,pU] = ttest(smallD.Unc.Performance,largeD.Unc.Performance);
[~,pC] = ttest(smallD.Cro.Performance,largeD.Cro.Performance);



figure;
makePlotPerfDiffUnityLine([smallD.Unc.Performance;smallD.Cro.Performance],...
    [largeD.Unc.Performance; largeD.Cro.Performance], pU, pC);
suptitle('Offset')

figure;

errorbar([1 2],([mean(smallD.Unc.Performance)-mean(meanPerf(1,:)),...
    mean(largeD.Unc.Performance)-mean(meanPerf(1,:))]),...
    [sem(smallD.Unc.Performance),...
    sem(largeD.Unc.Performance)],'-d','Color','b',...
    'MarkerFaceColor','b','MarkerSize',12);
hold on
errorbar([1 2],([mean(smallD.Cro.Performance)-mean(meanPerf(2,:)),...
    mean(largeD.Cro.Performance)-mean(meanPerf(2,:))]),[sem(smallD.Cro.Performance),...
    sem(largeD.Cro.Performance)],'-d','Color','r',...
    'MarkerFaceColor','r','MarkerSize',12);
xlim([0.5 2.5])
xticks([1 2])
xticklabels({'Small','Large'})
xlabel('Fixation DC')
ylabel('Change in Performance from Threshold')
line([0.5 2.5],[0 0],'LineStyle','--');
% yticks([.8 .9 1])
% ylim([.78 1.05])

%% Look at Speed instead of Area
counter = 1;
for ii = 1:size(dataAll)
    if length(dataAll.X{counter})<499
        Fz = 341;
    else
        Fz = 1000;
    end
    tmp_x = sgfilt(dataAll.X{counter},3,41,1);
    tmp_y = sgfilt(dataAll.Y{counter},3,41,1);
    x = tmp_x(floor((41/2)+10):end-floor(41/2))*Fz;
    y = tmp_y(floor((41/2)+10):end-floor(41/2))*Fz;
    vel_tmp = sqrt(x .^2 + y .^ 2) ;
    speed_tmp = sqrt(x .^ 2+ y .^ 2);
    idxSpd = speed_tmp < 100;
    dataAll.Speed(counter) = nanmean(speed_tmp);
    counter = counter + 1;
end

for ii = 1:length(unique([dataAll.SubjectID(:)]))
    idx = [];
    idx = find([dataAll.SubjectID(:)] == ii & ...
            [dataAll.Perf(:)] < 3 & ...
            [dataAll.Condition(:)] == -1);
    areaPolyD{1,ii} = [dataAll.Speed(idx)];
    meanPerf(1,ii) = mean([dataAll.Perf(idx)]);
    idx = [];
    idx = find([dataAll.SubjectID(:)] == ii & ...
            [dataAll.Perf(:)] < 3 & ...
            [dataAll.Condition(:)] == 1);
    areaPolyD{2,ii} = [dataAll.Speed(idx)];
    meanPerf(2,ii) = mean([dataAll.Perf(idx)]);
end



smallD = [];
largeD = [];
indSubBinUpper = ones(1,length(unique([dataAll.SubjectID(:)])))*30;
indSubBinLower = ones(1,length(unique([dataAll.SubjectID(:)])))*70;

[smallD,largeD] = calculatePerfProbAtThresh(indSubBinUpper,indSubBinLower,unique([dataAll.SubjectID(:)]),...
    dataAll,areaPolyD,'Speed',1,smallD,largeD);

indSubBinUpper = ones(1,length(unique([dataAll.SubjectID(:)])))*30;
indSubBinLower = ones(1,length(unique([dataAll.SubjectID(:)])))*70;

[smallD,largeD] = calculatePerfProbAtThresh(indSubBinUpper,indSubBinLower,unique([dataAll.SubjectID(:)]),...
    dataAll,areaPolyD,'Speed',2,smallD,largeD);

[~,pU,~,statsU] = ttest(smallD.Unc.Performance,largeD.Unc.Performance);
[~,pC,~,statsC] = ttest(smallD.Cro.Performance,largeD.Cro.Performance);



figure;
makePlotPerfDiffUnityLine([smallD.Unc.Performance;smallD.Cro.Performance],...
    [largeD.Unc.Performance; largeD.Cro.Performance], pU, pC);
suptitle('Speed')

figure;

errorbar([1 2],([mean(smallD.Unc.Performance)-mean(meanPerf(1,:)),...
    mean(largeD.Unc.Performance)-mean(meanPerf(1,:))]),...
    [sem(smallD.Unc.Performance),...
    sem(largeD.Unc.Performance)],'-d','Color','b',...
    'MarkerFaceColor','b','MarkerSize',12);
hold on
errorbar([1 2],([mean(smallD.Cro.Performance)-mean(meanPerf(2,:)),...
    mean(largeD.Cro.Performance)-mean(meanPerf(2,:))]),[sem(smallD.Cro.Performance),...
    sem(largeD.Cro.Performance)],'-d','Color','r',...
    'MarkerFaceColor','r','MarkerSize',12);
xlim([0.5 2.5])
xticks([1 2])
xticklabels({'Small','Large'})
xlabel('Fixation Speed')
ylabel('Change in Performance from Threshold')
line([0.5 2.5],[0 0],'LineStyle','--');
% yticks([.8 .9 1])
% ylim([.78 1.05])

%% Area in each condition
for ii = 1:length(unique([dataAll.SubjectID(:)]))
    idxTrials = find(dataAll.SubjectID == ii & ...
        dataAll.Condition == 1 & ...
        dataAll.Area < 200);
    
    meanAreaCr(ii) = mean(dataAll.Area(idxTrials));
    
    idxTrials = find(dataAll.SubjectID == ii & ...
        dataAll.Condition == -1 & ...
        dataAll.Area < 200);
    
    meanAreaUn(ii) = mean(dataAll.Area(idxTrials));
end

nanmean(meanAreaCr)
nanmean(meanAreaUn)
nanstd(meanAreaCr)
nanstd(meanAreaUn)

[h,p,ci,stats] = ttest(meanAreaCr, meanAreaUn);
stats
p


% %% Velocity as a function of crowding
% % for ii = 1:length(unique([dataAll.SubjectID(:)]))
% %     figure;
% %     subplot(1,2,1)
% %     
%     idxTrials = find(dataAll.Condition == -1 & ...
%         dataAll.Area < 200);
% %    mdlu = fitlm([dataAll.Speed(idxTrials)],[dataAll.Perf(idxTrials)])
%       mdlu = fitglm(dataAll.Speed(idxTrials),...
%           [dataAll.Perf(idxTrials) dataAll.SubjectID(idxTrials)],...
%            "linear","Distribution","binomial","link","logit")
% 
% %    plot(mdlu);
% %    
% %    
% %    
% %    subplot(1,2,2)
% %    idxTrials = find(dataAll.SubjectID == ii & ...
% %         dataAll.Condition == 1 & ...
% %         dataAll.Area < 200);
% %    mdlc = fitlm([dataAll.Speed(idxTrials)],[dataAll.Perf(idxTrials)])
% %    plot(mdlc);
% % end
% 
% nanmean(meanAreaCr)
% nanmean(meanAreaUn)
% 
% [h,p,ci,stats] = ttest(meanAreaCr, meanAreaUn);
% stats
% p

%% Testing Crowding between the two systems

crowdingEffect = thresh(:,2)-thresh(:,1)

[h,p,ci,stats] = ttest2(crowdingEffect(1:5),crowdingEffect(6:end));

%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function prob = coneAndFlankerProbability_OG(x,stimSize,multiplier)


stimSize = round(stimSize);%;2;%(arcmin)
spacingCones = 0.5;% (arcmin)
stimSizeInCones = stimSize / spacingCones;

targetActivation = abs(round(x) + round( (x-round(x))/0.5) * 0.5);
targetActRight = targetActivation + (stimSizeInCones/2);
targetActLeft = targetActivation - (stimSizeInCones/2);

xFlanker = targetActivation + (stimSize * (multiplier/.5));
flankerActivation = abs(round(xFlanker) + ...
    round( (xFlanker-round(xFlanker))/0.5) * 0.5);
flankerActRight= flankerActivation + (stimSizeInCones/2);
flankerActLeft = flankerActivation - (stimSizeInCones/2);

[targetAll,~, numTimesT] = unique([targetActivation targetActRight targetActLeft]);
[flankerAll,~, numTimesF] = unique([flankerActivation flankerActRight flankerActLeft]);

a_countsT = accumarray(numTimesT,1);
targetAllThresh = targetAll(find(a_countsT'>4));

a_countsF = accumarray(numTimesF,1);
flankerAllThresh = flankerAll(find(a_countsF'>5));

sameCones = intersect(targetAllThresh, flankerAllThresh);
prob = length(sameCones)/length(targetAllThresh);
end

function [result,limit] = generateHeatMapSimple( xValues, yValues, varargin )
%Creates 2D Distribution Map of Traces
n_bins = 40;
axisValue = 40;
stimuliSize = 1;
offset = 0;
doubleTarget = 0;
condition = 0;
borders = 1;
rotateIm = 0;
% poiEndName = 'TimeTargetOFF';
% axisWindow = 60;
% filepath = pwd;
% newFileName = 'pptrials';
% conversionFactor = 1; %machine sampling rate conversion factor (ie DPI = 1, dDPI = 1000/330)
% trialId = 1:length(pptrials);

k = 1;
Properties = varargin;
while k <= length(Properties) && ischar(Properties{k})
    switch (Properties{k})
        case 'Bins'
            n_bins =  Properties{k + 1};
            Properties(k:k+1) = [];
        case 'StimulusSize'
            stimuliSize = Properties{k + 1};
            Properties(k:k+1) = [];
        case 'AxisValue'
            axisValue = Properties{k + 1};
            Properties(k:k+1) = [];
        case 'Offset'
            offset = Properties{k + 1};
            Properties(k:k+1) = [];
        case 'DoubleTargets'
            doubleTarget = Properties{k + 1};
            Properties(k:k+1) = [];
        case 'Uncrowded' %(should be 1 for U, 2 for Crowded, 4 = Fixation)
            condition = Properties{k + 1};
            Properties(k:k+1) = [];
        case 'Borders'
            borders = Properties{k + 1};
            Properties(k:k+1) = [];
        case 'Rotate'
            rotateIm = Properties{k + 1};
            Properties(k:k+1) = [];
        otherwise
            k = k + 1;
    end
end

idx = find(xValues > -axisValue & ...
    xValues < axisValue &...
    yValues > -axisValue & ...
    yValues < axisValue);

xValues = xValues(idx);
yValues = yValues(idx);

limit.xmin = floor(min(xValues));
limit.xmax = ceil(max(xValues));
limit.ymin = floor(min(yValues));
limit.ymax = ceil(max(yValues));

%  limit.xmin = -60;
%         limit.xmax = 60;
%         limit.ymin = -60;
%         limit.ymax = 60;

result = MyHistogram2(xValues, yValues, [limit.xmin,limit.xmax,n_bins;limit.ymin,limit.ymax,n_bins]);
result = result./(max(max(result)));

if rotateIm ~= 0
    result = imrotate((result),rotateIm);
end

load('./MyColormaps.mat');
mycmap(1,:) = [1 1 1];
set(gcf, 'Colormap', mycmap)

hold on
temp = pcolor(linspace(limit.xmin, limit.xmax, size(result, 1)),...
    linspace(limit.ymin, limit.ymax, size(result, 1)),...
    result');

if stimuliSize == 0
    if condition == 0
    else
        centerX = 8;
        centerY = -8;
        width = 16;
        height = 16;
        rectangle('Position',[-centerX+stimuliSize, centerY, width, height],'LineWidth',3,'LineStyle','-')
    end
else
    stimuliSize = stimuliSize/2;
    width = 2 * stimuliSize;
    height = width * 5;
    centerX = (-stimuliSize+offset);
    centerY = (-stimuliSize*5);
    %     rectangle('Position',[centerX, centerY, width, height],'LineWidth',2)
    
    if doubleTarget
        % if uEcc(numEcc)> 0
        rectangle('Position',[-centerX+stimuliSize, centerY, width, height],'LineWidth',2,'LineStyle','--')
        % end
    end
    
    %     spacing = 1.4;
    %     arcminEcc = uEcc * params.pixelAngle;
    if condition == 2
        
        rectangle('Position',[-(width + (width * 1.4)) + offset, centerY, width, height],'LineWidth',1) %Right
        
        rectangle('Position',[(width * 1.4) + offset, centerY, width, height], 'LineWidth',1) %Left
        
        rectangle('Position',[centerX, -(height + (height * 1.4)), width, height], 'LineWidth',1) % Bottom
        
        rectangle('Position',[centerX, (height * 1.4), width, height], 'LineWidth',1) %Top
        
    end
    if condition == 1 || condition == 2
        rectangle('Position',[centerX, centerY, width, height],'LineWidth',3)
    end
    if condition == 4
        p = plot(0,0,'--ks','MarkerSize',stimuliSize);
        p(1).LineWidth = 3;
    end
end

set(gca, 'FontSize', 12)

if borders == 0
    set(gca,'xtick',[], 'ytick', []);
end

caxis([.1 ceil(max(max(result)))])
shading interp%/flat

axis([-axisValue axisValue -axisValue axisValue]);


end

function [smallD,largeD] = calculatePerfProbAtThresh(indSubBinUpper,indSubBinLower,subjectsAll,...
    dataAll,var,var2,c,smallD,largeD)

if c == 1
    cVal = -1;
    name = 'Unc';
elseif c ==2
    cVal = 1;
    name = 'Cro';
end

for ii = 1:length(subjectsAll)
    
    temp1 = find([dataAll.SubjectID] == ii &...
        [dataAll.Condition] == cVal);
    
    newIdx = find([dataAll.Size(temp1)] == mode(dataAll.Size(temp1)));
   
    sizePick = mode(dataAll.Size(temp1));
    
    forPercent(1,ii) = indSubBinUpper(ii);
    forPercent(2,ii) = indSubBinLower(ii);
    bin(ii,1) = prctile(var{c, ii},forPercent(1,ii),"all");
    bin(ii,2) = prctile(var{c, ii},forPercent(2,ii),"all");
    
    smallD.(name).TrialsIdx = find([dataAll.SubjectID] == ii &...
        [dataAll.Condition] == cVal & ...
        [dataAll.(var2)] < bin(ii,2)&...
        [dataAll.Size] == sizePick);
    
    largeD.(name).TrialsIdx = find([dataAll.SubjectID] == ii &...
        [dataAll.Condition] == cVal & ...
        [dataAll.(var2)] > bin(ii,1)&...
        [dataAll.Size] == sizePick);
    
    smallD.(name).ConeProb(ii) = mean(dataAll.SameConeProb(smallD.(name).TrialsIdx));
    largeD.(name).ConeProb(ii) = mean(dataAll.SameConeProb(largeD.(name).TrialsIdx));
    smallD.(name).ConeProbSEM(ii) = sem(dataAll.SameConeProb(smallD.(name).TrialsIdx));
    largeD.(name).ConeProbSEM(ii) = sem(dataAll.SameConeProb(largeD.(name).TrialsIdx));
    
    smallD.(name).stimSize(ii) = sizePick;
    largeD.(name).stimSize(ii) = sizePick;
    
    smallD.(name).trialNum(ii) = length(smallD.(name).TrialsIdx);
    largeD.(name).trialNum(ii) = length(largeD.(name).TrialsIdx);
    smallD.(name).Performance(ii) = mean(dataAll.Perf(smallD.(name).TrialsIdx));
    largeD.(name).Performance(ii) = mean(dataAll.Perf(largeD.(name).TrialsIdx));
    
    smallD.(name).NumTrials(ii) = length(dataAll.Perf(smallD.(name).TrialsIdx));
    largeD.(name).NumTrials(ii) = length(dataAll.Perf(largeD.(name).TrialsIdx));
    
    smallD.(name).varz(ii) = mean(dataAll.(var2)(smallD.(name).TrialsIdx));
    largeD.(name).varz(ii) = mean(dataAll.(var2)(largeD.(name).TrialsIdx));
    
end
end

function makePlotPerfDiffUnityLine(smallPerf,largePerf, pU, pC)

figure;
plot(smallPerf(2,:),largePerf(2,:),'o','Color','b');
line([.4 1],[.4 1])
subplot(1,2,1)
plot(smallPerf(1,:),largePerf(1,:),'o','Color','b');
axis([0 1 0 1])
xlabel({'Performance in Trials','for Smaller Fixation Area'});
ylabel({'Performance in Trials','for Larger Fixation Area'});
line([0 1],[0 1],'Color',[.2 .2 .2])
hold on
errorbar(mean(smallPerf(1,:)),...
mean(largePerf(1,:)),sem(smallPerf(1,:)),...
sem(largePerf(1,:)),sem(smallPerf(1,:)),...
sem(largePerf(1,:)),'o','MarkerSize',10,'Color','k');
title(sprintf('Uncrowded, p = %.3f', pU))
text(.05,.9,{'better performance with','large fixation area'},...
    'FontSize',8)
text(.4,.1,{'better performance with','small fixation area'},...
    'FontSize',8)
axis square

subplot(1,2,2)
plot(smallPerf(2,:),largePerf(2,:),'o','Color','r');
axis([0 1 0 1])
xlabel({'Performance in Trials','for Smaller Fixation Area'});
ylabel({'Performance in Trials','for Larger Fixation Area'});
line([0 1],[0 1],'Color',[.2 .2 .2])
hold on
errorbar(mean(smallPerf(2,:)),...
mean(largePerf(2,:)),sem(smallPerf(2,:)),sem(smallPerf(2,:)),...
sem(largePerf(2,:)),sem(largePerf(2,:)),'o','MarkerSize',10,'Color','k')
title(sprintf('Crowded, p = %.3f',pC))
text(.05,.9,{'better performance with','large fixation area'},...
    'FontSize',8)
text(.4,.1,{'better performance with','small fixation area'},...
    'FontSize',8)
axis square
end

