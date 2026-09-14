% Selects the best combination of gap-filling reaction sets for every
% organism with a completed AltRxns file, using NICEgame and matTFA (Salvy
% et al., 2019; Vayena et al., 2022). Tests combinations of alternative
% reaction sets across every carbon source and builds the final gapfilled
% model from the winning combination.
%
% Inputs:
%   data/genomes/reference_key/AtLSPHERE_RefSeq.mat                       - StrainRefSeqKey (organism -> RefSeq accession)
%   data/experimental/medium/min_med_CSourceScreen.mat                    - minMed, nutrients, vitamins
%   data/experimental/carbon_source_screening/binarizerd_CSourceScreen_Jun2024.xlsx - growth/no-growth ground truth
%   data/models/draft/cobrapy/mat/<accession>.mat                         - draft models
%   data/gapfilling/alt_rxnsets/AltRxns_<organism>.mat                    - alternative reaction sets, from gapfilling_alt_rxnsets.m
%   scripts/gapfilling/NICEgame/databases/BiGG_universal.mat
%   scripts/gapfilling/NICEgame/matTFA-master/    - matTFA toolbox
%   scripts/gapfilling/NICEgame/gapfilling/       - NICEgame core functions
%   scripts/gapfilling/core/                      - growAltSolnModel, growModelInCSources
%
% Outputs:
%   data/models/gapfilled/<method>/<organism>_gapfilled.mat            - gapfilled COBRA model
%   data/gapfilling/selection_results/<method>/<organism>_selection.mat - scoreCombo, combosToTest, predGrowth_raw/_norm/_bin, gapfilled_bin
%
% Author: Alan R. Pacheco, Clement Lefebvre
% NICEgame: Evangelia Vayena
% Modified by: Sudharshan Ravi

%%
clear all; clc;

%% Parameters
mandatoryGrowthNutrients = {}; % Limit solutions to those that enable growth on these carbon sources (if the organism grows on them experimentally and if solutions that enable growth on them have an FPR of at most mandatoryFPRCutoff times that of best solutions). Only applies to scoringMethod = 'mcc'.
mandatoryFPRCutoff = 1.5; % Upper limit of false positive rate (relative to that of best overall solution) to accept for solutions that allow growth on all mandatory nutrients

maxChooseK = 3; % Maximum number of alternative solutions to combine
maxBestCombos = 10; % Maximum number of combined alternative solutions to test via FBA

growthDataMode = 'binary'; % 'binary' or 'continuous'. 'continuous' requires both growth data files below; 'binary' requires only binaryGrowthDataFile.
scoringMethod = 'mcc'; % 'mcc' | 'auc' | 'pearson' | 'nrmse' | 'wnrmse' | 'log-nrmse' | 'smape' | 'all'.
                       % 'all' expands to {'mcc','auc'} when growthDataMode = 'binary', or every method when 'continuous'.
expMetric = 'growthRate'; % Only used when growthDataMode = 'continuous': 'growthRate' | 'maxOD' | 'endpointOD'.

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
rootDir = fullfile(scriptDir,'..','..','..');
dataDir = fullfile(rootDir,'data');
ngDir = fullfile(rootDir,'scripts','gapfilling','NICEgame');
coreDir = fullfile(rootDir,'scripts','gapfilling','core');

refSeqKeyFile = fullfile(dataDir,'genomes','reference_key','AtLSPHERE_RefSeq.mat');
modelDir = fullfile(dataDir,'models','draft','cobrapy','mat');
mediaFile = fullfile(dataDir,'experimental','medium','min_med_CSourceScreen.mat');
altRxnsDir = fullfile(dataDir,'gapfilling','alt_rxnsets');

% Growth data. binaryGrowthDataFile is always required. continuousGrowthDataFile
% is only read when growthDataMode = 'continuous', and must be an .xlsx file with the
% same layout as binaryGrowthDataFile (same carbon source columns, same order),
% with continuous growth rate instead of binary.
binaryGrowthDataFile = fullfile(dataDir,'experimental','carbon_source_screening','binarizerd_CSourceScreen_Jun2024.xlsx');
continuousGrowthDataFile = '';

modelDBFile = fullfile(ngDir,'databases','BiGG_universal.mat');
thermoDBFile = fullfile(ngDir,'matTFA-master','thermoDatabases','thermo_data.mat');
compartmentDataFile = fullfile(ngDir,'matTFA-master','pytfa','models','CompartmentData.mat');
compoundsFile = fullfile(ngDir,'databases','Compounds.mat');
pathMatTFA = fullfile(ngDir,'matTFA-master','matTFA');

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

outModelBase = fullfile(dataDir,'models','gapfilled');
outResultBase = fullfile(dataDir,'gapfilling','selection_results');
for mi = 1:nM
    if ~exist(fullfile(outModelBase,methods{mi}),'dir'); mkdir(fullfile(outModelBase,methods{mi})); end
    if ~exist(fullfile(outResultBase,methods{mi}),'dir'); mkdir(fullfile(outResultBase,methods{mi})); end
end

%% Load files
disp('Loading files...')

addpath(genpath(pathMatTFA));
addpath(genpath(fullfile(ngDir,'gapfilling')));
addpath(coreDir);

load(refSeqKeyFile)

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

disp('    Loading model database file...')
load(modelDBFile)
modelDB = model; clear model;
if ~isfield(modelDB,'metCompSymbol')
    metCompSymbol = cell(length(modelDB.mets),1);
    metComps = zeros(length(modelDB.mets),1);
    for m = 1:length(modelDB.mets)
        met = modelDB.mets{m};
        symbol = met(end);
        metCompSymbol{m} = symbol;
        if strcmp(symbol,'c')
            metComps(m) = 1;
        elseif strcmp(symbol,'p')
            metComps(m) = 2;
        elseif strcmp(symbol,'e')
            metComps(m) = 3;
        end
    end
    modelDB.metComps = metComps;
    modelDB.metCompSymbol = metCompSymbol;
end
if isfield(modelDB,'metCharges')
    modelDB.metCharge = modelDB.metCharges;
    modelDB = rmfield(modelDB,'metCharges');
end
modelDB = removeRxns(modelDB,modelDB.rxns(find(contains(lower(modelDB.rxnNames),'mevalonate')))); % Remove mevalonate-containing reactions

disp('    Loading thermodynamics database...')
flagTFA = 1;
load(thermoDBFile);
thermoData = DB_AlbertyUpdate;

disp('    Loading compartment files...')
load(compartmentDataFile)
load(compoundsFile)

%% Check if COBRA toolbox is loaded
if ~exist('optimizeCbModel','file'); initCobraToolbox(false,'agent'); end
changeCobraSolver('gurobi');

%% Select the best combination for every organism with a completed AltRxns file
altFiles = dir(fullfile(altRxnsDir,'AltRxns_*.mat'));

for MF = 1:length(altFiles)

[~,fname] = fileparts(altFiles(MF).name);
organismID = fname(9:end);

fprintf(['\n-------------------------------------------------------------------------------------\n\nModel: ' organismID '\n\n'])

organismIndex = find(strcmp(StrainRefSeqKey.Strain,organismID));
if isempty(organismIndex)
    warning(['Could not find a RefSeq ID for organism ' organismID '. Skipping...'])
    continue
end
modelFile = fullfile(modelDir,[StrainRefSeqKey.RefSeqAssembly{organismIndex} '.mat']);
if ~isfile(modelFile)
    warning(['Draft model not found for organism ' organismID '. Skipping...'])
    continue
end

trueGrowth = growthData(ismember(strainNames,organismID),:);

% Load the model
disp('    Loading draft model (note: COBRA may flag errors)...')
model = readCbModel(modelFile);
excRxns = find(contains(model.rxns,'EX_'));

% Change metabolite compartment suffix
for i = 1:length(model.mets)
    model.mets{i} = replace(model.mets{i},'[C_','_');
    model.mets{i} = replace(model.mets{i},']','');
end

% Change charge fieldname (matTFA uses singular 'metCharge')
disp('    Formatting model...')
model.metCharge = model.metCharges;
model = rmfield(model,'metCharges');

% Add missing exchange reactions to the model
excRxnMets = regexprep(model.rxns(excRxns), 'R_EX_|EX_', '');
missingExcMetsNutrients = nutrients(find(ismember(nutrients,excRxnMets) == 0));
missingExcMetsMinMed = minMed(find(ismember(minMed,excRxnMets) == 0));
missingExcMets = [missingExcMetsNutrients';missingExcMetsMinMed];
if isempty(find(ismember(excRxnMets,'co2_e'), 1))
    missingExcMets = [missingExcMets;'co2_e'];
end
for m = 1:length(missingExcMets)

    if isempty(find(ismember(model.mets,missingExcMets{m}), 1))
        metIndexModelDB = find(ismember(modelDB.mets,missingExcMets{m}));
        if isempty(metIndexModelDB)
            warning(['Metabolite ' missingExcMets{m} ' not found in modelDB. Cannot add its exchange reaction, skipping.'])
            continue
        end
        model.mets = [model.mets;missingExcMets{m}];
        model.metNames = [model.metNames;modelDB.metNames(metIndexModelDB)];
        model.metFormulas = [model.metFormulas;modelDB.metFormulas(metIndexModelDB)];
        model.metCharge = [model.metCharge;modelDB.metCharge(metIndexModelDB)];
        model.b = [model.b;modelDB.b(metIndexModelDB)];
        model.S = [model.S;zeros(1,length(model.rxns))];
        model.csense = [model.csense;'E'];
    end

    model = addReaction(model,['EX_' missingExcMets{m}],'reactionFormula',strcat(missingExcMets{m},' <=>'));
end
pl = '.'; if length(missingExcMets) > 1; pl = 's.'; end
if ~isempty(missingExcMets); disp(['        Added ' num2str(length(missingExcMets)) ' missing exchange reaction' pl]); end

% Add missing fields to model
model.CompartmentData = CompartmentData;
[model.metSEEDID,model.metCompSymbol] = deal(cell(length(model.mets),1));
for i = 1:length(model.mets)

    s = {model.mets{i}(1:end-2)};
    s1 = split(s{1},'__');
    if length(s1) > 1
        metname = strcat(s1{1},'-',s1{2});
    else
        metname = s{1};
    end

    f = find(ismember(Compounds.abbreviation,metname));
    if length(f) > 0
        model.metSEEDID(i) = Compounds.id(f(1));
    else
        f2 = find(contains(Compounds.aliases,metname));
        if length(f2) > 0
          model.metSEEDID(i) = Compounds.id(f2(1));
        else
            model.metSEEDID{i} = metname;
        end
    end
    model.metCompSymbol{i} = model.mets{i}(end);
end
model.rev = zeros(length(model.rxns),1);
model.rev(intersect(find(model.lb < 0), find(model.ub > 0))) = 1;

modelOrig = model;

% Nothing to test against: skip gap-filling, pass the draft model through.
% Still records the draft model's own FBA predictions, so every organism
% has an entry in every method's output regardless of branch.
if ~any(trueGrowth == 1)
    fprintf('    No carbon source this strain grows on -- nothing to gap-fill against. Passing through.\n');
    predGrowth_bin = growModelInCSources(modelOrig,minMed,vitamins,[],nutrients,10);
    predGrowth_raw = zeros(size(predGrowth_bin));
    predGrowth_norm = predGrowth_raw;
    gapfilled_bin = 0;
    modelSave = modelOrig;
    if isfield(modelSave,'modelNotes')
        modelSave.modelNotes = [modelSave.modelNotes(1:end-17) '<p>no carbon source with positive experimental growth -- not gap-filled</p> </html> </notes>'];
    else
        modelSave.modelNotes = '<notes><p>no carbon source with positive experimental growth -- not gap-filled</p> </html> </notes>';
    end
    modelSave.modelID = organismID;
    for mi = 1:nM
        model = modelSave;
        save(fullfile(outModelBase,methods{mi},[organismID '_gapfilled.mat']),'model','-v7.3');
        save(fullfile(outResultBase,methods{mi},[organismID '_selection.mat']), ...
            'gapfilled_bin','predGrowth_raw','predGrowth_norm','predGrowth_bin','-v7.3');
    end
    continue
end

% Load AltRxns
d = load(fullfile(altRxnsDir,altFiles(MF).name));
ActRxnsAll = d.ActRxnsAll;
csFields = fieldnames(ActRxnsAll);
numGFCS = length(csFields);
maxNumAlt = 0;
for i = 1:numGFCS
    ar = ActRxnsAll.(csFields{i});
    if ~isempty(ar) && ~isempty(ar{1,1})
        maxNumAlt = max(maxNumAlt,size(ar{1,1},1));
    end
end
if maxNumAlt == 0; maxNumAlt = 1; end

% Experimental reference for this strain
if strcmp(growthDataMode,'continuous')
    contIdx = find(strcmp(contStrainNames,organismID));
    if ~isempty(contIdx)
        expGrowth = contMatrix(contIdx,:);
    else
        expGrowth = trueGrowth;
    end
else
    expGrowth = trueGrowth;
end
expGrowth(isnan(expGrowth)) = 0;
maxExp = max(expGrowth);
if maxExp > 0; expGrowth_norm = expGrowth / maxExp;
else;          expGrowth_norm = zeros(size(expGrowth));
end

numTestCS = length(nutrients);

% Build growthRates3D (numGFCS x maxNumAlt x numTestCS) -- once per strain, shared across methods
growthRates3D = zeros(numGFCS,maxNumAlt,numTestCS);
predGrowth3D = zeros(numGFCS,maxNumAlt,numTestCS);

todisp = '';
for i = 1:numGFCS
    ar = ActRxnsAll.(csFields{i});
    if isempty(ar) || isempty(ar{1,1}); continue; end
    numAltCurr = size(ar{1,1},1);
    fprintf(repmat('\b',1,length(todisp)));
    todisp = sprintf('    %d/%d  %s',i,numGFCS,csFields{i});
    fprintf(todisp);
    [growthRates3D(i,1:numAltCurr,:),predGrowth3D(i,1:numAltCurr,:)] = ...
        growAltSolnModel(ar,modelDB,[1,numAltCurr],minMed,vitamins,[],nutrients,modelOrig,0);
end
fprintf('\n');

growthRates2D = reshape(growthRates3D,[numGFCS*maxNumAlt,numTestCS]);
predGrowth2D = reshape(predGrowth3D,[numGFCS*maxNumAlt,numTestCS]);
% Row r maps to: i_cs = mod(r-1,numGFCS)+1, n_alt = ceil(r/numGFCS)
validRows = find(any(growthRates2D > 1e-3,2));

% No valid gap-filling solutions anywhere: pass draft model through for every method
if isempty(validRows)
    warning([organismID ': no valid gap-filling solutions found.']);
    predGrowth_bin = growModelInCSources(modelOrig,minMed,vitamins,[],nutrients,10);
    predGrowth_raw = zeros(size(predGrowth_bin));
    predGrowth_norm = predGrowth_raw;
    gapfilled_bin = 0;
    modelSave = modelOrig;
    if isfield(modelSave,'modelNotes')
        modelSave.modelNotes = [modelSave.modelNotes(1:end-17) '<p>no gap-filling solutions found</p> </html> </notes>'];
    else
        modelSave.modelNotes = '<notes><p>no gap-filling solutions found</p> </html> </notes>';
    end
    modelSave.modelID = organismID;
    for mi = 1:nM
        model = modelSave;
        save(fullfile(outModelBase,methods{mi},[organismID '_gapfilled.mat']),'model','-v7.3');
        save(fullfile(outResultBase,methods{mi},[organismID '_selection.mat']), ...
            'gapfilled_bin','predGrowth_raw','predGrowth_norm','predGrowth_bin','-v7.3');
    end
    continue
end

mandatoryIdx = find(ismember(nutrients,mandatoryGrowthNutrients));

% Steps 2-5: run per scoring method (tensor above is shared)
for mi = 1:nM
    thisMethod = methods{mi};
    outModelDir = fullfile(outModelBase,thisMethod);
    outResultDir = fullfile(outResultBase,thisMethod);

    fprintf('  [%s]\n',thisMethod);

    allCombos = zeros(0,maxChooseK);
    for k = 1:maxChooseK
        if length(validRows) >= k
            ck = nchoosek(validRows,k);
            allCombos = [allCombos;[ck,zeros(size(ck,1),maxChooseK-k)]];
        end
    end

    fprintf('    %d valid rows, %d combos\n',length(validRows),size(allCombos,1));
    approxScores = zeros(size(allCombos,1),1);
    mandatoryOK = false(size(allCombos,1),1);

    for c = 1:size(allCombos,1)
        rows = allCombos(c,allCombos(c,:) > 0);
        if strcmp(thisMethod,'mcc')
            combinedPred = double(any(predGrowth2D(rows,:) > 0,1));
            approxScores(c) = mccScore(combinedPred,trueGrowth);
            if ~isempty(mandatoryIdx)
                mandatoryOK(c) = sum(combinedPred(mandatoryIdx)) > 0;
            end
        else
            approxPred = max(growthRates2D(rows,:),[],1);
            approxScores(c) = scoreVec(approxPred,expGrowth,expGrowth_norm,trueGrowth,thisMethod);
        end
    end

    % Optional mandatory-nutrient override (mcc only)
    if strcmp(thisMethod,'mcc') && ~isempty(mandatoryIdx) && any(mandatoryOK)
        bestFPR = fprAtBestScore(allCombos,approxScores,predGrowth2D,trueGrowth);
        mandatoryFPRs = zeros(sum(mandatoryOK),1);
        mIdx = find(mandatoryOK);
        for j = 1:length(mIdx)
            rows = allCombos(mIdx(j),allCombos(mIdx(j),:) > 0);
            combinedPred = double(any(predGrowth2D(rows,:) > 0,1));
            [~,~,fp,~] = confusionCounts(combinedPred,trueGrowth);
            tn = sum((combinedPred==0) & (trueGrowth==0));
            mandatoryFPRs(j) = fp / max(fp+tn,1);
        end
        if min(mandatoryFPRs) <= bestFPR * mandatoryFPRCutoff
            keep = false(size(allCombos,1),1);
            keep(mIdx) = true;
            allCombos = allCombos(keep,:);
            approxScores = approxScores(keep);
        end
    end

    % Select top maxBestCombos, deduplicated by reaction union
    [~,sortIdx] = sort(approxScores,'descend');
    selectedRows = zeros(maxBestCombos,maxChooseK);
    selectedRxnSets = {};
    numSelected = 0;

    for s = 1:length(sortIdx)
        c = sortIdx(s);
        rows = allCombos(c,:);
        rxnSet = rowsToRxnSet(rows(rows > 0),numGFCS,csFields,ActRxnsAll);

        isDup = false;
        for prev = 1:numSelected
            if isequal(rxnSet,selectedRxnSets{prev}); isDup = true; break; end
        end

        if ~isDup
            numSelected = numSelected + 1;
            selectedRows(numSelected,:) = rows;
            selectedRxnSets{numSelected} = rxnSet;
        end
        if numSelected >= maxBestCombos; break; end
    end
    selectedRows = selectedRows(1:numSelected,:);
    fprintf('    %d combos for FBA\n',numSelected);

    combosToTest = zeros(numSelected,maxChooseK*2);
    for p = 1:numSelected
        for j = 1:maxChooseK
            r = selectedRows(p,j);
            if r == 0; continue; end
            combosToTest(p,j*2-1) = mod(r-1,numGFCS) + 1;
            combosToTest(p,j*2) = ceil(r/numGFCS);
        end
    end

    [growthRatesCombos,predGrowthCombos] = deal(zeros(numSelected,numTestCS));
    for p = 1:numSelected
        fprintf('    combo %d/%d\n',p,numSelected);
        coords = combosToTest(p,find(combosToTest(p,:)));
        ActRxnsNew = buildActRxns(coords,csFields,ActRxnsAll);
        [growthRatesCombos(p,:),predGrowthCombos(p,:)] = ...
            growAltSolnModel(ActRxnsNew,modelDB,[1,1],minMed,vitamins,[],nutrients,modelOrig,0);
    end

    scoreCombo = zeros(numSelected,1);
    mccCombo = zeros(numSelected,1);

    for p = 1:numSelected
        predRates = growthRatesCombos(p,:);
        predBin = predGrowthCombos(p,:);

        if strcmp(thisMethod,'mcc')
            scoreCombo(p) = mccScore(predBin,trueGrowth);
        else
            scoreCombo(p) = scoreVec(predRates,expGrowth,expGrowth_norm,trueGrowth,thisMethod);
        end
        mccCombo(p) = mccScore(predBin,trueGrowth);
    end

    [~,bestIdx] = max(scoreCombo);
    bestIdx = bestIdx(1);
    bestCoords = combosToTest(bestIdx,find(combosToTest(bestIdx,:)));

    predGrowth_raw = growthRatesCombos(bestIdx,:);
    predGrowth_bin = predGrowthCombos(bestIdx,:);
    if max(predGrowth_raw) > 0
        predGrowth_norm = predGrowth_raw / max(predGrowth_raw);
    else
        predGrowth_norm = zeros(size(predGrowth_raw));
    end
    gapfilled_bin = 1;

    % Build gapfilled model from winning combination
    ActRxnsWinner = buildActRxns(bestCoords,csFields,ActRxnsAll);
    uniqueNewRxns = ActRxnsWinner{1,1}{1,1};
    uniqueNewFormulas = ActRxnsWinner{1,1}{1,2};

    model = modelOrig;
    for r = 1:length(uniqueNewRxns)
        s = split(uniqueNewFormulas{r},' ');
        metsToAdd = setdiff(s(contains(s,'_')),model.mets);
        for m = 1:length(metsToAdd)
            metIdx = find(ismember(modelDB.mets,metsToAdd{m}));
            if isempty(metIdx); continue; end
            model.mets = [model.mets;metsToAdd{m}];
            model.metNames = [model.metNames;modelDB.metNames(metIdx)];
            model.metFormulas = [model.metFormulas;modelDB.metFormulas(metIdx)];
            model.metCharge = [model.metCharge;modelDB.metCharge(metIdx)];
            model.b = [model.b;modelDB.b(metIdx)];
            model.S = [model.S;zeros(1,length(model.rxns))];
            model.csense = [model.csense;'E'];
            if isfield(model,'metCompSymbol')
                model.metCompSymbol = [model.metCompSymbol;metsToAdd{m}(end)];
            end
        end
        model = addReaction(model,['R_' uniqueNewRxns{r}],'reactionFormula',uniqueNewFormulas{r});
    end

    for i = 1:length(model.mets)
        met = model.mets{i};
        if length(met) < 2; continue; end
        sfx = met(end-1:end);
        nm = met(1:end-2);
        sfx = strrep(sfx,'_c','[c]');
        sfx = strrep(sfx,'_p','[p]');
        sfx = strrep(sfx,'_e','[e]');
        model.mets{i} = [nm sfx];
    end

    winnerDesc = sprintf('gapfilled on %s (alt %d)',csFields{bestCoords(1)},bestCoords(2));
    for q = 2:length(bestCoords)/2
        winnerDesc = [winnerDesc sprintf(' + %s (alt %d)',csFields{bestCoords(q*2-1)},bestCoords(q*2))];
    end
    winnerDesc = [winnerDesc sprintf('; %s = %.4f; MCC = %.4f',thisMethod,scoreCombo(bestIdx),mccCombo(bestIdx))];

    if isfield(model,'modelNotes')
        model.modelNotes = [model.modelNotes(1:end-17) '<p>' winnerDesc '</p> </html> </notes>'];
    else
        model.modelNotes = ['<notes><p>' winnerDesc '</p> </html> </notes>'];
    end
    model.modelID = organismID;

    save(fullfile(outModelDir,[organismID '_gapfilled.mat']),'model','-v7.3');
    save(fullfile(outResultDir,[organismID '_selection.mat']), ...
        'scoreCombo','mccCombo','growthRatesCombos','predGrowthCombos', ...
        'combosToTest','bestIdx','csFields','thisMethod','expMetric', ...
        'predGrowth_raw','predGrowth_norm','predGrowth_bin','gapfilled_bin','-v7.3');
end

end

%% Remove path
rmpath(genpath(pathMatTFA));

%% Local functions

function rxnSet = rowsToRxnSet(rows,numGFCS,csFields,ActRxnsAll)
    rxnSet = {};
    for idx = 1:length(rows)
        r = rows(idx);
        if r == 0; continue; end
        i_cs = mod(r-1,numGFCS) + 1;
        n_alt = ceil(r/numGFCS);
        ar = ActRxnsAll.(csFields{i_cs});
        if isempty(ar) || isempty(ar{1,1}); continue; end
        if n_alt > size(ar{1,1},1); continue; end
        rxnSet = [rxnSet;ar{1,1}{n_alt,1}];
    end
    rxnSet = sort(unique(rxnSet));
end

function ActRxnsNew = buildActRxns(coords,csFields,ActRxnsAll)
    allNewRxns = ActRxnsAll.(csFields{coords(1)}){1,1}{coords(2),1};
    allNewFormulas = ActRxnsAll.(csFields{coords(1)}){1,1}{coords(2),2};
    for q = 2:length(coords)/2
        allNewRxns = [allNewRxns;ActRxnsAll.(csFields{coords(q*2-1)}){1,1}{coords(q*2),1}];
        allNewFormulas = [allNewFormulas;ActRxnsAll.(csFields{coords(q*2-1)}){1,1}{coords(q*2),2}];
    end
    uniqueRxns = unique(allNewRxns);
    uniqueFormulas = cell(length(uniqueRxns),1);
    for r = 1:length(uniqueRxns)
        match = allNewFormulas(ismember(allNewRxns,uniqueRxns(r)));
        uniqueFormulas(r) = match(1);
    end
    ActRxnsNew = {};
    ActRxnsNew{1,1}{1,1} = uniqueRxns;
    ActRxnsNew{1,1}{1,2} = uniqueFormulas;
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

function fpr = fprAtBestScore(allCombos,approxScores,predGrowth2D,trueGrowth)
    [~,bestC] = max(approxScores);
    rows = allCombos(bestC,allCombos(bestC,:) > 0);
    combinedPred = double(any(predGrowth2D(rows,:) > 0,1));
    [~,tn,fp,~] = confusionCounts(combinedPred,trueGrowth);
    fpr = fp / max(fp+tn,1);
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
