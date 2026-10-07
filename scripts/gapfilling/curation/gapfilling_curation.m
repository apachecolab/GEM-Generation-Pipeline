% Two-stage iterative false-positive/false-negative curation of gap-filled
% models, run once per scoring method against that method's own gap-filled
% models. Stage 1 runs until no organism's score improves for `patience`
% iterations; stage 2 re-seeds each organism from its own best and stops on
% the first regression.
%
% Inputs:
%   data/experimental/medium/min_med_CSourceScreen.mat                    - minMed, nutrients, vitamins
%   data/experimental/carbon_source_screening/binarizerd_CSourceScreen_Jun2024.xlsx - growth/no-growth ground truth
%   data/models/gapfilled/<method>/<organism>_gapfilled.mat               - gap-filled COBRA models, from gapfilling_selection.m
%   <donorDir>/*.mat                                                      - external donor COBRA models, only if donorDir is set
%   <donorDir's parent>/essential/*.mat                                   - from donor_essentiality_sweep.m, only if donorReactionMode = 'essential'
%   scripts/gapfilling/core/                                             - growModelInCSources, normalizeMetIDs
%
% Outputs:
%   data/models/curated/<method>/<organism>_curated.mat                 - curated COBRA model
%   data/gapfilling/curation_results/<method>/modelAccuracy.mat          - accData: per-organism confusion matrix and metrics, latest iteration
%   data/gapfilling/curation_results/<method>/iterationTracking.mat      - tracking: per-metric score, organism x iteration, both stages concatenated
%   data/gapfilling/donor_cache/<donorDir basename>/growth.mat           - cached donor predicted-growth matrix, only if donorDir is set
%
% Author: Alan R. Pacheco, Clement Lefebvre
% NICEgame: Evangelia Vayena
% Modified by: Sudharshan Ravi

clear all; clc;

%% Parameters
scoringMethod = 'mcc'; % 'mcc' | 'auc' | 'pearson' | 'nrmse' | 'wnrmse' | 'log-nrmse' | 'smape' | 'all'.
                       % 'all' expands to {'mcc','auc'} when growthDataMode = 'binary', or every method when 'continuous'.
growthDataMode = 'binary'; % 'binary' or 'continuous'. 'continuous' requires both growth data files below; 'binary' requires only binaryGrowthDataFile.
expMetric = 'growthRate'; % Only used when growthDataMode = 'continuous': 'growthRate' | 'maxOD' | 'endpointOD'.

donorDir = fullfile('data','donors','fba'); % relative to the repository root. '' = internal donors (this run's own gap-filled models). Otherwise a directory of external .mat donor models.
donorReactionMode = 'essential'; % 'fba' | 'essential'. Only used when donorDir is set.

maxK = 3; % Maximum number of new reactions to combine
maxNumCombos = 10000; % Maximum number of combinations to test (in order to limit computational time). Set to Inf for unlimited.
maxIter = 10; % Maximum number of curation iterations, per stage
patience = 2; % consecutive non-improving iterations (global) before stage 1 stops

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
rootDir = fullfile(scriptDir,'..','..','..');
if ~isempty(donorDir); donorDir = fullfile(rootDir,donorDir); end
dataDir = fullfile(rootDir,'data');
coreDir = fullfile(rootDir,'scripts','gapfilling','core');
addpath(coreDir);

mediaFile = fullfile(dataDir,'experimental','medium','min_med_CSourceScreen.mat');

% Growth data. binaryGrowthDataFile is always required. continuousGrowthDataFile
% is only read when growthDataMode = 'continuous', and must be an .xlsx file with the
% same layout as binaryGrowthDataFile (same carbon source columns, same order),
% with continuous growth rate instead of binary.
binaryGrowthDataFile = fullfile(dataDir,'experimental','carbon_source_screening','binarizerd_CSourceScreen_Jun2024.xlsx');
continuousGrowthDataFile = '';

isExternalDonor = ~isempty(donorDir);
if isExternalDonor
    donorDirClean = regexprep(donorDir,[filesep '+$'],'');
    [~,donorDirName] = fileparts(donorDirClean);
    donorCacheDir = fullfile(dataDir,'gapfilling','donor_cache',donorDirName,filesep);
    if ~exist(donorCacheDir,'dir'); mkdir(donorCacheDir); end
    donorGrowthCacheFile = fullfile(donorCacheDir,'growth.mat');
    donorEssentialDir = fullfile(fileparts(donorDirClean),'essential',filesep); % sibling of donorDir, written by donor_essentiality_sweep.m
end

allMethods = {'mcc','auc','pearson','nrmse','wnrmse','log-nrmse','smape'};
binaryMethods = {'mcc','auc'};
if strcmp(scoringMethod,'all')
    if strcmp(growthDataMode,'continuous')
        methods = allMethods;
    else
        methods = binaryMethods;
    end
else
    if strcmp(growthDataMode,'binary') && ~ismember(scoringMethod,binaryMethods)
        error('scoringMethod ''%s'' requires growthDataMode = ''continuous''.',scoringMethod);
    end
    methods = {scoringMethod};
end
nM = length(methods);

outModelBase = fullfile(dataDir,'models','curated');
outResultBase = fullfile(dataDir,'gapfilling','curation_results');
for mi = 1:nM
    if ~exist(fullfile(outModelBase,methods{mi}),'dir'); mkdir(fullfile(outModelBase,methods{mi})); end
    if ~exist(fullfile(outResultBase,methods{mi}),'dir'); mkdir(fullfile(outResultBase,methods{mi})); end
end

%% Load the medium and experimental growth data
disp('Loading files...')

disp('    Loading medium file...')
load(mediaFile)
minMed = replace(minMed,'[e]','_e');
nutrients = replace(nutrients,'[e]','_e');
vitamins = replace(vitamins,'[e]','_e');

disp('    Loading experimental data file...')
if ~isfile(binaryGrowthDataFile)
    error(['Binarized growth data not found: ' binaryGrowthDataFile]);
end
[~, ~, raw] = xlsread(binaryGrowthDataFile);
rawHeader = raw(1,2:end);
strainNames = raw(2:end,1);
nutrients = raw(1,2:end);
growthData = raw(2:end,2:end);
growthData = reshape([growthData{:}],size(growthData));
waterIndex = find(ismember(nutrients,'h2o'));
nutrients(waterIndex) = [];
nutrientNames(waterIndex) = [];
growthData(:,waterIndex) = [];
growthData(isnan(growthData)) = 0;
nutrients = strcat(nutrients,'_e');

contStrainNames = {};
contMatrix = [];
if strcmp(growthDataMode,'continuous')
    if ~isfile(continuousGrowthDataFile)
        error('growthDataMode is ''continuous'' but continuousGrowthDataFile was not found.');
    end
    [~, ~, rawCont] = xlsread(continuousGrowthDataFile);
    if ~isequal(rawCont(1,2:end),rawHeader)
        error('continuousGrowthDataFile''s carbon source columns must match binaryGrowthDataFile exactly, same names, same order.');
    end
    contStrainNames = rawCont(2:end,1);
    contMatrix = rawCont(2:end,2:end);
    contMatrix = reshape([contMatrix{:}],size(contMatrix));
    contMatrix(:,waterIndex) = []; % Keep columns aligned with nutrients, which has water removed
end

%% Check if COBRA toolbox is loaded
if ~exist('optimizeCbModel','file'); initCobraToolbox(false,'agent'); end
changeCobraSolver('gurobi');

%% External donor pool: discover, precompute/cache predicted growth, snapshot essentiality data
donorIDs = {};
donorPredGrowth = [];
donorPredGrowthCount = [];
essentialData = containers.Map();
if isExternalDonor
    donorFiles = dir(fullfile(donorDir,'*.mat'));
    donorIDs = erase({donorFiles.name},'.mat');

    if exist(donorGrowthCacheFile,'file')
        disp('Loading cached donor predicted-growth matrix...')
        cached = load(donorGrowthCacheFile); % donorIDs, donorPredGrowth
        if isequal(cached.donorIDs,donorIDs)
            donorPredGrowth = cached.donorPredGrowth;
        else
            disp('donorDir listing changed since the cache was written. Recomputing...')
        end
    end
    if isempty(donorPredGrowth)
        fprintf('Precomputing predicted growth for %d donor(s)...\n', length(donorIDs));
        donorPredGrowth = zeros(length(donorIDs),length(nutrients));
        for i = 1:length(donorIDs)
            fprintf('    %d/%d  %s\n', i, length(donorIDs), donorIDs{i});
            d = load(fullfile(donorDir,donorFiles(i).name));
            donorModel = normalizeMetIDs(d.model);
            donorPredGrowth(i,:) = growModelInCSources(donorModel,minMed,vitamins,[],nutrients,10);
        end
        save(donorGrowthCacheFile,'donorIDs','donorPredGrowth','-v7.3')
    end
    donorPredGrowthCount = sum(donorPredGrowth,2);

    if strcmp(donorReactionMode,'essential')
        essentialFiles = dir(fullfile(donorEssentialDir,'*.mat'));
        sweptDonorIDs = erase({essentialFiles.name},'.mat');
        for i = 1:length(sweptDonorIDs)
            d = load(fullfile(donorEssentialDir,[sweptDonorIDs{i} '.mat']),'results');
            essentialData(sweptDonorIDs{i}) = d.results;
        end
        fprintf('%d/%d donor(s) have essentiality data on disk.\n', length(sweptDonorIDs), length(donorIDs));
    end
end

donorCtx = struct('isExternalDonor',isExternalDonor,'donorDir',donorDir,'donorIDs',{donorIDs}, ...
    'donorPredGrowth',donorPredGrowth,'donorPredGrowthCount',donorPredGrowthCount, ...
    'donorReactionMode',donorReactionMode,'essentialData',essentialData);

%% ==================== Per scoring method ====================
for mi = 1:nM
    thisMethod = methods{mi};
    fprintf(['\n\n############################\n' thisMethod '\n############################\n'])

    modelDir = fullfile(dataDir,'models','gapfilled',thisMethod,filesep);
    workingDir = fullfile(dataDir,'gapfilling','curation_working',thisMethod,filesep); % scratch: deleted once final output is written
    bestSnapshotDir = fullfile(workingDir,'best_snapshot',filesep);
    bestPerStrainDir = fullfile(workingDir,'best_per_strain',filesep);
    finalDir = fullfile(outModelBase,thisMethod,filesep);
    resultsDir = fullfile(outResultBase,thisMethod,filesep);
    if ~exist(workingDir,'dir'); mkdir(workingDir); end
    if ~exist(bestSnapshotDir,'dir'); mkdir(bestSnapshotDir); end
    if ~exist(bestPerStrainDir,'dir'); mkdir(bestPerStrainDir); end
    accuracySaveFileName = fullfile(resultsDir,'modelAccuracy');
    trackingSaveFileName = fullfile(resultsDir,'iterationTracking');

    %% Order the organisms by strain name
    allFiles = dir(fullfile(modelDir,'*_gapfilled.mat'));
    organismIDs = erase({allFiles.name},'_gapfilled.mat');

    %% Initialise the working dir from the gap-filled models
    disp('Initialising working dir from gap-filled models...')
    for i = 1:length(organismIDs)
        src = [modelDir organismIDs{i} '_gapfilled.mat'];
        if exist(src,'file'); copyfile(src,[workingDir organismIDs{i} '.mat']); end
    end

    bestScore_prev = [];
    stagnantCount = 0;
    lastIter = 1;
    maxIterTotal = maxIter + maxIter; % stage 1 + stage 2, tracked as one continuous trajectory
    tracking = struct();

    %% ==================== STAGE 1 ====================
    for iter = 1:maxIter

        fprintf(['\n\n============================\nIteration ' num2str(iter) '\n============================\n'])

        if iter == 1
            sourceDir = modelDir;
            sourceSuffix = '_gapfilled.mat';
        else
            sourceDir = workingDir;
            sourceSuffix = '.mat';
        end

        [organismIDsCurr,predGrowthAll,trueGrowthAll,growthRatesAll,scores] = computeAccuracy( ...
            organismIDs,sourceDir,sourceSuffix,nutrients,minMed,vitamins,growthData,strainNames, ...
            growthDataMode,contMatrix,contStrainNames);

        driverScore_All = scores.(driverField(thisMethod));

        accData = scores; accData.nutrients = nutrients; accData.organismIDs = organismIDsCurr;
        accData.predGrowthAll = predGrowthAll; accData.trueGrowthAll = trueGrowthAll; accData.growthRatesAll = growthRatesAll;
        save(accuracySaveFileName,'accData','-v7.3')

        if iter == 1
            tracking.organismIDs = organismIDsCurr;
            tracking.nutrients = nutrients;
            for f = fieldnames(rmfield(scores,'confusion'))'
                tracking.(f{1}) = NaN(length(organismIDsCurr),maxIterTotal);
            end
            bestScore_prev = -Inf(length(organismIDsCurr),1);
        end
        [~,idxMap] = ismember(organismIDsCurr,tracking.organismIDs);
        for f = fieldnames(rmfield(scores,'confusion'))'
            tracking.(f{1})(idxMap,iter) = scores.(f{1});
        end
        save(trackingSaveFileName,'tracking','-v7.3')

        %% Convergence check -- patience-based, global stop
        improvedNow = driverScore_All > bestScore_prev(idxMap);
        bestScore_prev(idxMap(improvedNow)) = driverScore_All(improvedNow);

        improvedOrgIdx = find(improvedNow);
        for k = 1:length(improvedOrgIdx)
            orgName = organismIDsCurr{improvedOrgIdx(k)};
            copyfile([sourceDir orgName sourceSuffix],[bestPerStrainDir orgName '.mat'])
        end

        if any(improvedNow)
            stagnantCount = 0;
            copyfile([workingDir '*.mat'],bestSnapshotDir)
        else
            stagnantCount = stagnantCount + 1;
        end
        lastIter = iter;
        if stagnantCount >= patience
            disp(['No organism improved for ' num2str(patience) ' iteration(s). Converged after ' num2str(iter) ' iteration(s).'])
            break
        end

        troubleshootModels(organismIDsCurr,predGrowthAll,trueGrowthAll,growthRatesAll,sourceDir,sourceSuffix,workingDir, ...
            nutrients,minMed,vitamins,growthData,strainNames,growthDataMode,contMatrix,contStrainNames, ...
            maxK,maxNumCombos,donorCtx);

    end % stage 1

    %% ==================== STAGE 2: refine from each organism's own best ====================
    stage1LastIter = lastIter;
    disp(['Stage 1 complete after ' num2str(stage1LastIter) ' iteration(s). Starting stage 2 from each organism''s own best...'])

    for i = 1:length(organismIDs)
        src = [bestPerStrainDir organismIDs{i} '.mat'];
        if exist(src,'file'); copyfile(src,[workingDir organismIDs{i} '.mat']); end
    end

    sourceDir = workingDir;
    sourceSuffix = '.mat';
    stage2Baseline = [];

    for iter2 = 1:maxIter
        iter = stage1LastIter + iter2;

        fprintf(['\n\n============================\nStage 2, Iteration ' num2str(iter2) ' (global iteration ' num2str(iter) ')\n============================\n'])

        [organismIDsCurr,predGrowthAll,trueGrowthAll,growthRatesAll,scores] = computeAccuracy( ...
            organismIDs,sourceDir,sourceSuffix,nutrients,minMed,vitamins,growthData,strainNames, ...
            growthDataMode,contMatrix,contStrainNames);

        driverScore_All = scores.(driverField(thisMethod));

        accData = scores; accData.nutrients = nutrients; accData.organismIDs = organismIDsCurr;
        accData.predGrowthAll = predGrowthAll; accData.trueGrowthAll = trueGrowthAll; accData.growthRatesAll = growthRatesAll;
        save(accuracySaveFileName,'accData','-v7.3')

        [~,idxMap] = ismember(organismIDsCurr,tracking.organismIDs);
        for f = fieldnames(rmfield(scores,'confusion'))'
            tracking.(f{1})(idxMap,iter) = scores.(f{1});
        end
        save(trackingSaveFileName,'tracking','-v7.3')

        if iter2 == 1
            stage2Baseline = -Inf(length(tracking.organismIDs),1);
        end
        regressedNow = driverScore_All < stage2Baseline(idxMap);
        if any(regressedNow)
            disp('Stage 2: an organism dropped below its own best reached so far in stage 2. Stopping.')
            break
        end
        improvedNow2 = driverScore_All > stage2Baseline(idxMap);
        if iter2 > 1 && ~any(improvedNow2)
            disp('Stage 2: no organism improved this iteration. Stopping.')
            break
        end
        stage2Baseline(idxMap) = max(stage2Baseline(idxMap),driverScore_All);
        copyfile([workingDir '*.mat'],bestSnapshotDir)
        lastIter = iter;

        modelsToTroubleshoot = organismIDsCurr(find(any(predGrowthAll~=trueGrowthAll,2)));
        if isempty(modelsToTroubleshoot)
            disp('Stage 2: nothing left to troubleshoot. Stopping.')
            break
        end

        troubleshootModels(organismIDsCurr,predGrowthAll,trueGrowthAll,growthRatesAll,sourceDir,sourceSuffix,workingDir, ...
            nutrients,minMed,vitamins,growthData,strainNames,growthDataMode,contMatrix,contStrainNames, ...
            maxK,maxNumCombos,donorCtx);

    end % stage 2

    %% Final save: every organism from the last coherent snapshot
    for i = 1:length(tracking.organismIDs)
        orgID = tracking.organismIDs{i};
        d = load([bestSnapshotDir orgID '.mat'],'model');
        model = bracketizeMetIDs(d.model); %#ok<NASGU>
        save([finalDir orgID '_curated.mat'],'model','-v7.3')
    end
    if isfolder(workingDir); rmdir(workingDir,'s'); end

    for f = fieldnames(tracking)'
        if ~ismember(f{1},{'organismIDs','nutrients'})
            tracking.(f{1}) = tracking.(f{1})(:,1:lastIter);
        end
    end
    save(trackingSaveFileName,'tracking','-v7.3')

    disp([thisMethod ' curation complete.'])

end % method loop

curationWorkingRoot = fullfile(dataDir,'gapfilling','curation_working');
if isfolder(curationWorkingRoot) && numel(dir(curationWorkingRoot)) <= 2 % dir() always lists '.' and '..'
    rmdir(curationWorkingRoot)
end

%% Local functions

function fieldName = driverField(thisMethod)
    switch thisMethod
        case 'mcc';       fieldName = 'MCC_All';
        case 'auc';       fieldName = 'AUC_All';
        case 'pearson';   fieldName = 'Pearson_All';
        case 'nrmse';     fieldName = 'NRMSE_All';
        case 'wnrmse';    fieldName = 'wNRMSE_All';
        case 'log-nrmse'; fieldName = 'LogNRMSE_All';
        case 'smape';     fieldName = 'SMAPE_All';
    end
end

function [organismIDsCurr,predGrowthAll,trueGrowthAll,growthRatesAll,scores] = computeAccuracy( ...
    organismIDs,sourceDir,sourceSuffix,nutrients,minMed,vitamins,growthData,strainNames, ...
    growthDataMode,contMatrix,contStrainNames)

organismIDsCurr = organismIDs;
[Acc_All,AccBal_All,TPR_All,FPR_All,MCC_All,AUC_All] = deal(zeros(length(organismIDsCurr),1));
[Pearson_All,NRMSE_All,wNRMSE_All,LogNRMSE_All,SMAPE_All] = deal(zeros(length(organismIDsCurr),1));
confusion = zeros(length(organismIDsCurr),length(nutrients));
organismsNotPresent = [];
[predGrowthAll,trueGrowthAll,growthRatesAll] = deal(zeros(length(organismIDsCurr),length(nutrients)));

for i = 1:length(organismIDsCurr)
    try
        d = load([sourceDir organismIDsCurr{i} sourceSuffix]);
    catch
        organismsNotPresent = [organismsNotPresent i]; %#ok<AGROW>
        warning(['No file found for organism ' organismIDsCurr{i} '. Skipping...'])
        continue
    end
    disp(['Calculating accuracy for ' organismIDsCurr{i} '...'])
    d.model = normalizeMetIDs(d.model);

    [predGrowth,growthRates] = growModelInCSources(d.model,minMed,vitamins,[],nutrients,10);
    trueGrowth = growthData(find(ismember(strainNames,organismIDsCurr{i})),:);

    [tp,tn,fp,fn] = confusionCounts(predGrowth,trueGrowth);
    Acc_All(i) = length(find(predGrowth == trueGrowth))/length(nutrients);
    TPR_All(i) = tp/max(tp+fn,1);
    FPR_All(i) = fp/max(fp+tn,1);
    AccBal_All(i) = (tp/max(tp+fn,1) + tn/max(fp+tn,1))/2;
    MCC_All(i) = mccScore(predGrowth,trueGrowth);
    AUC_All(i) = scoreVec(growthRates,[],[],trueGrowth,'auc');

    confusion(i,intersect(find(predGrowth==1),find(trueGrowth==1))) = 1;
    confusion(i,intersect(find(predGrowth==0),find(trueGrowth==0))) = -1;
    confusion(i,intersect(find(predGrowth==1),find(trueGrowth==0))) = 0.5;
    confusion(i,intersect(find(predGrowth==0),find(trueGrowth==1))) = -0.5;

    predGrowthAll(i,:) = predGrowth;
    trueGrowthAll(i,:) = trueGrowth;
    growthRatesAll(i,:) = growthRates;

    if strcmp(growthDataMode,'continuous')
        contIdx = find(ismember(contStrainNames,organismIDsCurr{i}));
        if ~isempty(contIdx); expGrowth = contMatrix(contIdx,:); else; expGrowth = trueGrowth; end
        expGrowth(isnan(expGrowth)) = 0;
        maxExp = max(expGrowth);
        if maxExp > 0; expGrowth_norm = expGrowth/maxExp; else; expGrowth_norm = zeros(size(expGrowth)); end

        Pearson_All(i) = scoreVec(growthRates,expGrowth,expGrowth_norm,trueGrowth,'pearson');
        NRMSE_All(i) = scoreVec(growthRates,expGrowth,expGrowth_norm,trueGrowth,'nrmse');
        wNRMSE_All(i) = scoreVec(growthRates,expGrowth,expGrowth_norm,trueGrowth,'wnrmse');
        LogNRMSE_All(i) = scoreVec(growthRates,expGrowth,expGrowth_norm,trueGrowth,'log-nrmse');
        SMAPE_All(i) = scoreVec(growthRates,expGrowth,expGrowth_norm,trueGrowth,'smape');
    end
end

organismIDsCurr(organismsNotPresent) = [];
Acc_All(organismsNotPresent) = []; AccBal_All(organismsNotPresent) = [];
TPR_All(organismsNotPresent) = []; FPR_All(organismsNotPresent) = [];
MCC_All(organismsNotPresent) = []; AUC_All(organismsNotPresent) = [];
Pearson_All(organismsNotPresent) = []; NRMSE_All(organismsNotPresent) = [];
wNRMSE_All(organismsNotPresent) = []; LogNRMSE_All(organismsNotPresent) = [];
SMAPE_All(organismsNotPresent) = [];
confusion(organismsNotPresent,:) = [];
predGrowthAll(organismsNotPresent,:) = [];
trueGrowthAll(organismsNotPresent,:) = [];
growthRatesAll(organismsNotPresent,:) = [];

scores = struct('Acc_All',Acc_All,'AccBal_All',AccBal_All,'TPR_All',TPR_All,'FPR_All',FPR_All, ...
    'MCC_All',MCC_All,'AUC_All',AUC_All,'Pearson_All',Pearson_All,'NRMSE_All',NRMSE_All, ...
    'wNRMSE_All',wNRMSE_All,'LogNRMSE_All',LogNRMSE_All,'SMAPE_All',SMAPE_All,'confusion',confusion);
end

function troubleshootModels(organismIDsCurr,predGrowthAll,trueGrowthAll,growthRatesAll,sourceDir,sourceSuffix,workingDir, ...
    nutrients,minMed,vitamins,growthData,strainNames,growthDataMode,contMatrix,contStrainNames, ...
    maxK,maxNumCombos,donorCtx)

modelsToTroubleshoot = organismIDsCurr(find(any(predGrowthAll~=trueGrowthAll,2)));

for mmm = 1:length(modelsToTroubleshoot)
    ticStart = tic;
    troubleshootOrganismID = modelsToTroubleshoot{mmm};
    fprintf(['\nTroubleshooting model: ' troubleshootOrganismID '\n\n'])
    d = load([sourceDir troubleshootOrganismID sourceSuffix]);
    modelProblem = normalizeMetIDs(d.model);
    troubleshootRow = find(ismember(organismIDsCurr,troubleshootOrganismID));
    trueGrowth = growthData(find(ismember(strainNames,troubleshootOrganismID)),:);
    expGrowth = trueGrowth;
    if strcmp(growthDataMode,'continuous')
        contIdx = find(ismember(contStrainNames,troubleshootOrganismID));
        if ~isempty(contIdx); expGrowth = contMatrix(contIdx,:); end
        expGrowth(isnan(expGrowth)) = 0;
    end
    fixedFlag = 0;

    falseNegCSources = nutrients(intersect(find(predGrowthAll(troubleshootRow,:)==0),find(trueGrowthAll(troubleshootRow,:)==1)));

    for qqq = 1:length(falseNegCSources)
        problemCSource = falseNegCSources{qqq};
        csi = find(ismember(nutrients,problemCSource));

        %% Donor selection
        if ~donorCtx.isExternalDonor
            candidateIndices = intersect(find(predGrowthAll(:,csi)==1),find(trueGrowthAll(:,csi)==1));
            if isempty(candidateIndices)
                fprintf(['    Could not troubleshoot ' problemCSource '. No organisms to reference!\n'])
                continue
            end
            specialistIndices = candidateIndices(sum(predGrowthAll(candidateIndices,:),2) == min(sum(predGrowthAll(candidateIndices,:),2)));
            if length(specialistIndices) > 1
                if strcmp(growthDataMode,'continuous')
                    gaps = zeros(length(specialistIndices),1);
                    for k = 1:length(specialistIndices)
                        predRanksD = rankDescend(growthRatesAll(specialistIndices(k),:));
                        donorContIdx = find(ismember(contStrainNames,organismIDsCurr{specialistIndices(k)}));
                        if ~isempty(donorContIdx)
                            donorExpGrowth = contMatrix(donorContIdx,:);
                        else
                            donorExpGrowth = trueGrowthAll(specialistIndices(k),:);
                        end
                        trueRanksD = rankDescend(donorExpGrowth);
                        gaps(k) = abs(predRanksD(csi) - trueRanksD(csi));
                    end
                    [~,bestSpecialist] = min(gaps);
                else
                    numRxnsNeeded = zeros(length(specialistIndices),1);
                    for k = 1:length(specialistIndices)
                        dK = load([sourceDir organismIDsCurr{specialistIndices(k)} sourceSuffix]);
                        modelK = applyMedium(normalizeMetIDs(dK.model),minMed,vitamins,problemCSource);
                        FBAsolnK = optimizeCbModel(modelK);
                        numRxnsNeeded(k) = length(setdiff(modelK.rxns(find(FBAsolnK.x ~= 0)),modelProblem.rxns));
                    end
                    [~,bestSpecialist] = min(numRxnsNeeded);
                end
                referenceOrganismID = organismIDsCurr{specialistIndices(bestSpecialist)};
            else
                referenceOrganismID = organismIDsCurr{specialistIndices};
            end
            dRef = load([sourceDir referenceOrganismID sourceSuffix]);
            modelRef = normalizeMetIDs(dRef.model);
        else
            candidateDonorIndices = find(donorCtx.donorPredGrowth(:,csi)==1);
            if strcmp(donorCtx.donorReactionMode,'essential')
                candidateDonorIndices = candidateDonorIndices(arrayfun(@(idx) ...
                    ~isempty(essentialRxnsFor(donorCtx.essentialData,donorCtx.donorIDs{idx},problemCSource)),candidateDonorIndices));
            end
            if isempty(candidateDonorIndices)
                fprintf(['    Could not troubleshoot ' problemCSource '. No external donors predicted to grow!\n'])
                continue
            end
            specialistDonorIndices = candidateDonorIndices(donorCtx.donorPredGrowthCount(candidateDonorIndices) == min(donorCtx.donorPredGrowthCount(candidateDonorIndices)));
            if length(specialistDonorIndices) > 1
                numRxnsNeeded = zeros(length(specialistDonorIndices),1);
                for k = 1:length(specialistDonorIndices)
                    donorID = donorCtx.donorIDs{specialistDonorIndices(k)};
                    if strcmp(donorCtx.donorReactionMode,'essential')
                        numRxnsNeeded(k) = length(setdiff(essentialRxnsFor(donorCtx.essentialData,donorID,problemCSource),modelProblem.rxns));
                    else
                        dK = load(fullfile(donorCtx.donorDir,[donorID '.mat']));
                        donorModelK = normalizeMetIDs(dK.model);
                        donorModelK = applyMedium(donorModelK,minMed,vitamins,problemCSource);
                        FBAsolnK = optimizeCbModel(donorModelK);
                        numRxnsNeeded(k) = length(setdiff(donorModelK.rxns(find(FBAsolnK.x ~= 0)),modelProblem.rxns));
                    end
                end
                [~,bestSpecialist] = min(numRxnsNeeded);
                referenceOrganismID = donorCtx.donorIDs{specialistDonorIndices(bestSpecialist)};
            else
                referenceOrganismID = donorCtx.donorIDs{specialistDonorIndices};
            end
            dRef = load(fullfile(donorCtx.donorDir,[referenceOrganismID '.mat']));
            modelRef = normalizeMetIDs(dRef.model);
        end

        problemCSourceName = modelRef.metNames{find(ismember(modelRef.mets,problemCSource))};
        problemCSourceName = lower(problemCSourceName);
        problemCSourceName = replace(problemCSourceName,'l-','L-');
        problemCSourceName = replace(problemCSourceName,'d-','D-');
        fprintf(['    Troubleshooting ' problemCSourceName '. Reference model: ' referenceOrganismID '\n'])

        %% Candidate reactions from the donor
        if donorCtx.isExternalDonor && strcmp(donorCtx.donorReactionMode,'essential')
            testReactions = setdiff(essentialRxnsFor(donorCtx.essentialData,referenceOrganismID,problemCSource),modelProblem.rxns);
        else
            modelRefBounded = applyMedium(modelRef,minMed,vitamins,problemCSource);
            FBAsoln = optimizeCbModel(modelRefBounded);
            testReactions = setdiff(modelRefBounded.rxns(find(FBAsoln.x ~= 0)),modelProblem.rxns);
        end

        if isempty(testReactions); continue; end

        %% Exhaustive single-reaction test
        workedSingleReactions = {};
        for r = 1:length(testReactions)
            modelTest = modelProblem;
            modelTest = addReactionFromRef(modelTest,modelRef,testReactions{r});
            modelTest = applyMedium(modelTest,minMed,vitamins,problemCSource);
            FBAsoln = optimizeCbModel(modelTest);
            if FBAsoln.f > 1e-3
                workedSingleReactions{end+1} = testReactions{r}; %#ok<AGROW>
            end
        end

        if ~isempty(workedSingleReactions)
            if length(workedSingleReactions) > 1
                if strcmp(growthDataMode,'continuous')
                    trueRanksR = rankDescend(expGrowth);
                    gapsR = zeros(length(workedSingleReactions),1);
                    for r = 1:length(workedSingleReactions)
                        modelTest = modelProblem;
                        modelTest = addReactionFromRef(modelTest,modelRef,workedSingleReactions{r});
                        [~,growthRatesR] = growModelInCSources(modelTest,minMed,vitamins,[],nutrients,10);
                        predRanksR = rankDescend(growthRatesR);
                        gapsR(r) = abs(predRanksR(csi) - trueRanksR(csi));
                    end
                    [~,bestSingle] = min(gapsR);
                else
                    newFalsePosCounts = zeros(length(workedSingleReactions),1);
                    for r = 1:length(workedSingleReactions)
                        modelTest = modelProblem;
                        modelTest = addReactionFromRef(modelTest,modelRef,workedSingleReactions{r});
                        predGrowthR = growModelInCSources(modelTest,minMed,vitamins,[],nutrients,10);
                        newFalsePosCounts(r) = length(intersect(find(predGrowthR==1),find(trueGrowth==0)));
                    end
                    [~,bestSingle] = min(newFalsePosCounts);
                end
            else
                bestSingle = 1;
            end
            reactionsToAdd = workedSingleReactions(bestSingle);
            fixedFlag = 1;
            modelProblem = addReactionFromRef(modelProblem,modelRef,reactionsToAdd{1});

            predGrowth = growModelInCSources(modelProblem,minMed,vitamins,[],nutrients,10);
            falsePositivesNew = nutrients(intersect(find(predGrowth==1),find(trueGrowth==0)));
            fprintf(['        Adding reaction ' reactionsToAdd{1} ' enables growth on ' problemCSourceName ' with the fewest new false positives (' num2str(length(falsePositivesNew)) ').\n\n'])
            continue
        end

        %% Escalate to combination search
        maxKCurr = min(maxK,length(testReactions));
        if maxKCurr < 2; continue; end
        numCombos = size(nchoosek(1:length(testReactions),maxKCurr),1);
        while numCombos > maxNumCombos
            maxKCurr = maxKCurr - 1;
            if maxKCurr < 2; break; end
            numCombos = size(nchoosek(1:length(testReactions),maxKCurr),1);
        end
        if maxKCurr < 2; continue; end
        if maxKCurr < min(maxK,length(testReactions))
            fprintf(['        Warning: maximum number of combinations to test exceeded. Reduced to combinations of ' num2str(maxKCurr) ' reactions.\n'])
        end

        testReactionCombos = zeros(0,maxKCurr);
        for KCurr = 2:maxKCurr
            combosK = nchoosek(1:length(testReactions),KCurr);
            testReactionCombos = [testReactionCombos; combosK, zeros(size(combosK,1),maxKCurr-KCurr)]; %#ok<AGROW>
        end

        workedSelectReactions = cell(size(testReactionCombos,1),maxKCurr);
        for r = 1:size(testReactionCombos,1)
            if (r > 1) && any(~cellfun(@isempty,workedSelectReactions(1:r-1,1))); break; end
            modelTest = modelProblem;
            idxCombo = testReactionCombos(r,testReactionCombos(r,:) > 0);
            reactionsToTest = testReactions(idxCombo);
            for rr = 1:length(reactionsToTest)
                modelTest = addReactionFromRef(modelTest,modelRef,reactionsToTest{rr});
            end
            modelTest = applyMedium(modelTest,minMed,vitamins,problemCSource);
            FBAsoln = optimizeCbModel(modelTest);
            if FBAsoln.f > 1e-3
                workedSelectReactions(r,1:length(reactionsToTest)) = reactionsToTest;
            end
        end
        workedSelectReactions(cellfun(@isempty,workedSelectReactions(:,1)),:) = [];

        if isempty(workedSelectReactions)
            fprintf(['        No combinations of selected reactions enable growth on ' problemCSourceName '. Adding all new reactions from reference model...\n'])
            modelTest = modelProblem;
            for r = 1:length(testReactions)
                modelTest = addReactionFromRef(modelTest,modelRef,testReactions{r});
            end
            modelTest = applyMedium(modelTest,minMed,vitamins,problemCSource);
            FBAsoln = optimizeCbModel(modelTest);
            if FBAsoln.f > 1e-3
                fixedFlag = 1;
                reactionsToAdd = testReactions;
            else
                fprintf(['        Adding all new reactions from reference model (' num2str(length(testReactions)) ') did not enable growth on ' problemCSourceName '.\n'])
                continue
            end
        else
            disp(['        Found ' num2str(size(workedSelectReactions,1)) ' reaction combination(s) that enable growth on ' problemCSourceName '.'])
            newFalsePosCounts = zeros(size(workedSelectReactions,1),1);
            for r = 1:size(workedSelectReactions,1)
                modelTest = modelProblem;
                reactionsToTest = workedSelectReactions(r,find(cellfun(@isempty,workedSelectReactions(r,:))==0));
                for rr = 1:length(reactionsToTest)
                    modelTest = addReactionFromRef(modelTest,modelRef,reactionsToTest{rr});
                end
                modelTest = applyMedium(modelTest,minMed,vitamins,problemCSource);
                FBAsoln = optimizeCbModel(modelTest);
                if FBAsoln.f > 1e-3
                    predGrowthR = growModelInCSources(modelTest,minMed,vitamins,[],nutrients,10);
                    newFalsePosCounts(r) = length(intersect(find(predGrowthR==1),find(trueGrowth==0)));
                end
            end
            [~,bestSoln] = min(newFalsePosCounts);
            fixedFlag = 1;
            reactionsToAdd = workedSelectReactions(bestSoln,find(cellfun(@isempty,workedSelectReactions(bestSoln,:))==0));
        end

        for r = 1:length(reactionsToAdd)
            modelProblem = addReactionFromRef(modelProblem,modelRef,reactionsToAdd{r});
        end
        predGrowth = growModelInCSources(modelProblem,minMed,vitamins,[],nutrients,10);
        falsePositivesNew = nutrients(intersect(find(predGrowth==1),find(trueGrowth==0)));
        fprintf(['        Adding ' num2str(length(reactionsToAdd)) ' reaction(s) enables growth on ' problemCSourceName ' with ' num2str(length(falsePositivesNew)) ' new false positives.\n\n'])
    end

    %% Correct false positives: remove reactions only active during incorrect growth
    predGrowth = growModelInCSources(modelProblem,minMed,vitamins,[],nutrients,10);
    falsePositives = nutrients(intersect(find(predGrowth==1),find(trueGrowth==0)));
    if ~isempty(falsePositives)
        fprintf('\n    Attempting to correct false positives...')

        activeRxnsTP = {};
        nutrientsTP = nutrients(find(trueGrowth));
        for nnn = 1:length(nutrientsTP)
            modelTest = modelProblem;
            modelTest = applyMedium(modelTest,minMed,vitamins,nutrientsTP(nnn));
            FBAsoln = optimizeCbModel(modelTest);
            if FBAsoln.f > 1e-3
                activeRxnsTP = [activeRxnsTP;modelTest.rxns(find(FBAsoln.x))]; %#ok<AGROW>
            end
        end
        activeRxnsTP = unique(activeRxnsTP);

        for nnn = 1:length(falsePositives)
            modelTest = modelProblem;
            modelTest = applyMedium(modelTest,minMed,vitamins,falsePositives(nnn));
            FBAsoln = optimizeCbModel(modelTest);
            if FBAsoln.f > 1e-3
                activeRxnsFP = modelTest.rxns(find(FBAsoln.x));
                activeRxnsFP = [activeRxnsFP;findRxnsFromMets(modelTest,falsePositives{nnn})];

                if isempty(find(ismember({'malt_e','arab__L_e','xyl__D_e','fru_e','glc__D_e','tre_e','cellb_e','sucr_e','glyc_e'},falsePositives{nnn}),1))
                    activeRxnsFP = setdiff(activeRxnsFP,modelTest.rxns(find(findExcRxns(modelTest))));
                    transportRxnsCurr = {};
                    for rrr = 1:length(activeRxnsFP)
                        if length(unique(modelTest.metCompSymbol(find(ismember(modelTest.mets,findMetsFromRxns(modelTest,activeRxnsFP(rrr))))))) > 1
                            transportRxnsCurr = [transportRxnsCurr;activeRxnsFP(rrr)]; %#ok<AGROW>
                        end
                    end
                    activeRxnsFP = setdiff(activeRxnsFP,transportRxnsCurr);
                end
                rxnsToRemove = setdiff(activeRxnsFP,activeRxnsTP);
                if isfield(modelTest,'manuallyAddedRxns')
                    rxnsToRemove = setdiff(rxnsToRemove,strcat('R_',modelTest.manuallyAddedRxns));
                end
                rxnsToRemove = setdiff(rxnsToRemove,'R_ATPM');

                modelTest = removeRxns(modelTest,rxnsToRemove);
                predGrowth = growModelInCSources(modelTest,minMed,vitamins,[],nutrients,10);
                falseNegativesNew = nutrients(intersect(find(predGrowth==0),find(trueGrowth==1)));

                if (length(intersect(find(predGrowth==1),find(trueGrowth==0))) < length(falsePositives)) && (length(intersect(find(predGrowth==0),find(trueGrowth==1))) <= length(falseNegativesNew))
                    modelProblem = removeRxns(modelProblem,rxnsToRemove);
                end
            end
        end
        fixedFlag = 1;
    end

    %% Save working copy
    if fixedFlag
        model = modelProblem;

        uniqueRxns = unique(model.rxns);
        for r = 1:length(uniqueRxns)
            dups = find(ismember(model.rxns,uniqueRxns{r}));
            if length(dups) > 1
                sVecNew = zeros(length(model.mets),1);
                for d2 = 1:length(dups)
                    sVecNew(find(model.S(:,dups(d2)))) = model.S(find(model.S(:,dups(d2))),dups(d2));
                end
                model.S(:,dups(1)) = sVecNew;
                selRxns = ones(length(model.rxns),1);
                selRxns(dups(2:end)) = 0;
                model = removeFieldEntriesForType(model, ~selRxns, 'rxns', numel(model.rxns));
            end
        end

        uniqueMets = unique(model.mets);
        for m = 1:length(uniqueMets)
            dups = find(ismember(model.mets,uniqueMets{m}));
            if length(dups) > 1
                sVecNew = zeros(1,length(model.rxns));
                for d2 = 1:length(dups)
                    sVecNew(find(model.S(dups(d2),:))) = model.S(dups(d2),find(model.S(dups(d2),:)));
                end
                model.S(dups(1),:) = sVecNew;
                selMets = ones(length(model.mets),1);
                selMets(dups(2:end)) = 0;
                model = removeFieldEntriesForType(model, ~selMets, 'mets', numel(model.mets));
            end
        end

        save([workingDir troubleshootOrganismID '.mat'],'model','-v7.3')

        elapsedT = toc(ticStart);
        fprintf(['\nModel ' troubleshootOrganismID ' complete. Time elapsed: ' num2str(round(elapsedT/60,2)) ' minutes.\n'])
    end
end
end

function [tp,tn,fp,fn] = confusionCounts(predBin,trueGrowth)
    tp = sum((predBin==1) & (trueGrowth==1));
    tn = sum((predBin==0) & (trueGrowth==0));
    fp = sum((predBin==1) & (trueGrowth==0));
    fn = sum((predBin==0) & (trueGrowth==1));
end

function s = mccScore(predBin,trueGrowth)
    [tp,tn,fp,fn] = confusionCounts(predBin,trueGrowth);
    denom = sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn));
    if denom > 0; s = ((tp*tn)-(fp*fn))/denom; else; s = 0; end
end

function s = scoreVec(predRates,expGrowth,expGrowth_norm,trueGrowth,method)
    s = 0;
    switch method
        case 'pearson'
            if std(predRates) > 0 && std(expGrowth) > 0
                s = corr(predRates',expGrowth','type','Pearson');
            end
        case 'nrmse'
            if max(predRates) > 0
                pred_norm = predRates / max(predRates);
                s = 1 - sqrt(mean((expGrowth_norm - pred_norm).^2));
            end
        case 'wnrmse'
            if max(predRates) > 0
                pred_norm = predRates / max(predRates);
                w = max(expGrowth_norm,pred_norm);
                w_sum = sum(w);
                if w_sum > 0
                    s = 1 - sqrt(sum(w .* (expGrowth_norm - pred_norm).^2) / w_sum);
                end
            end
        case 'log-nrmse'
            epsilon = 0.005;
            if max(predRates) > 0
                pred_norm = predRates / max(predRates);
                s = 1 - sqrt(mean((log(expGrowth_norm + epsilon) - log(pred_norm + epsilon)).^2));
            end
        case 'smape'
            epsilon = 1e-6;
            if max(predRates) > 0
                pred_norm = predRates / max(predRates);
                r = abs(pred_norm - expGrowth_norm) ./ (max(pred_norm,expGrowth_norm) + epsilon);
                s = 1 - mean(r);
            end
        case 'auc'
            pos_idx = find(trueGrowth == 1);
            n_pos = length(pos_idx);
            n_neg = sum(trueGrowth == 0);
            if max(predRates) > 0 && n_pos > 0 && n_neg > 0
                ranks = tiedrank(predRates);
                s = (sum(ranks(pos_idx)) - n_pos*(n_pos+1)/2) / (n_pos * n_neg);
            end
    end
end

function mdl = applyMedium(mdl,mm,vits,csource)
    mdl.lb(find(findExcRxns(mdl))) = 0;
    mdl.lb(find(ismember(mdl.rxns,intersect(mdl.rxns(find(findExcRxns(mdl))),findRxnsFromMets(mdl,mm))))) = -1000;
    for vi = 1:length(vits)
        vitRxns = intersect(mdl.rxns(find(findExcRxns(mdl))),findRxnsFromMets(mdl,vits(vi)));
        if ~isempty(vitRxns); mdl.lb(ismember(mdl.rxns,vitRxns)) = -1000; end
    end
    mdl.lb(find(ismember(mdl.rxns,intersect(mdl.rxns(find(findExcRxns(mdl))),findRxnsFromMets(mdl,csource))))) = -10;
end

function mdl = addReactionFromRef(mdl,refModel,rxnID)
% refModel may have no .rev field (modern COBRA schema) -- reversibility is
% then derived from bounds instead: model.rev(intersect(find(model.lb < 0),
% find(model.ub > 0))) = 1.
    rxnIdx = find(ismember(refModel.rxns,rxnID));
    if isfield(refModel,'rev')
        isRev = refModel.rev(rxnIdx) > 0;
    else
        isRev = refModel.lb(rxnIdx) < 0 && refModel.ub(rxnIdx) > 0;
    end
    mdl = addReaction(mdl,rxnID,'metaboliteList',refModel.mets(find(refModel.S(:,rxnIdx)))', ...
        'stoichCoeffList',refModel.S(find(refModel.S(:,rxnIdx)),rxnIdx), ...
        'reversible',isRev);
end

function ranks = rankDescend(vec)
% Rank 1 = highest value. Ties resolved by array order (stable sort).
    [~,sortIdx] = sort(vec,'descend');
    ranks = zeros(size(vec));
    ranks(sortIdx) = 1:length(vec);
end

function rxns = essentialRxnsFor(essentialData,donorID,carbonSource)
    rxns = {};
    if ~isKey(essentialData,donorID); return; end
    results = essentialData(donorID);
    eci = find(ismember({results.carbon_source},carbonSource));
    if isempty(eci) || isempty(results(eci).essential); return; end
    rxns = results(eci).rxns(results(eci).essential);
end
