% Enumerates alternative gap-filling reaction sets for every draft model
% found, across every carbon source, using NICEgame and matTFA (Salvy et
% al., 2019; Vayena et al., 2022). Draft models are merged with the BiGG
% universal database under thermodynamic constraints, then gfbiomass finds
% candidate reactions that would enable growth. Selection of the best
% combination happens in a separate script.
%
% Inputs:
%   data/genomes/reference_key/AtLSPHERE_RefSeq.mat                       - StrainRefSeqKey (organism -> RefSeq accession)
%   data/experimental/medium/min_med_CSourceScreen.mat                    - minMed, nutrients, vitamins
%   data/experimental/carbon_source_screening/binarizerd_CSourceScreen_Jun2024.xlsx - growth/no-growth ground truth
%   data/models/draft/cobrapy/mat/<accession>.mat                         - draft models
%   scripts/gapfilling/NICEgame/databases/BiGG_universal.mat
%   scripts/gapfilling/NICEgame/databases/Compounds.mat
%   scripts/gapfilling/NICEgame/matTFA-master/    - matTFA toolbox
%   scripts/gapfilling/NICEgame/gapfilling/       - NICEgame core functions
%
% Outputs:
%   data/gapfilling/alt_rxnsets/AltRxns_<organism>.mat - ActRxnsAll, nutrients
%
% Author: Alan R. Pacheco, Clement Lefebvre
% NICEgame: Evangelia Vayena
% Modified by: Sudharshan Ravi

%%
clear all; clc;

%% Parameters
forceO2Uptake = 1; % Sets a mandatory uptake flux for oxygen (for aerobic gapfilling)
maxNumAlt = 25; % Maximum number of alternative solutions per carbon source
maxSolveTime = 500; % Seconds per MILP solve; [] = no limit. Set e.g. 500 to bound runtime.
retryMode = 'binary'; % Order of carbon sources tried when gapfilling fails and retries with one forced into the medium. 'binary': true-positive list order. 'continuous': ranked by growth rate.

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
rootDir = fullfile(scriptDir,'..','..','..');
dataDir = fullfile(rootDir,'data');
ngDir = fullfile(rootDir,'scripts','gapfilling','NICEgame');

refSeqKeyFile = fullfile(dataDir,'genomes','reference_key','AtLSPHERE_RefSeq.mat');
modelDir = fullfile(dataDir,'models','draft','cobrapy','mat');
mediaFile = fullfile(dataDir,'experimental','medium','min_med_CSourceScreen.mat');

% Growth data. growthDataFile (binary growth/no-growth) is always required.
% growthRateFile is only read when retryMode = 'continuous', and must be an
% .xlsx file with the same layout as growthDataFile (same carbon source
% columns, same order), with continuous growth rate instead of binary.
growthDataFile = fullfile(dataDir,'experimental','carbon_source_screening','binarizerd_CSourceScreen_Jun2024.xlsx');
growthRateFile = '';
modelDBFile = fullfile(ngDir,'databases','BiGG_universal.mat');
thermoDBFile = fullfile(ngDir,'matTFA-master','thermoDatabases','thermo_data.mat');
compartmentDataFile = fullfile(ngDir,'matTFA-master','pytfa','models','CompartmentData.mat');
compoundsFile = fullfile(ngDir,'databases','Compounds.mat');
pathMatTFA = fullfile(ngDir,'matTFA-master','matTFA');

saveDirGFSoln = fullfile(dataDir,'gapfilling','alt_rxnsets');
if ~exist(saveDirGFSoln,'dir'); mkdir(saveDirGFSoln); end

%% Find draft models and map to organism names
disp('Finding draft models...')
load(refSeqKeyFile)
modelFileList = dir(fullfile(modelDir,'*.mat'));
modelFiles = cell(length(modelFileList),1);
organismIDs = cell(length(modelFileList),1);
notFound = [];
for i = 1:length(modelFileList)
    accession = erase(modelFileList(i).name,'.mat');
    organismIndex = find(strcmp(StrainRefSeqKey.RefSeqAssembly,accession));
    if ~isempty(organismIndex)
        modelFiles{i} = accession;
        organismIDs{i} = StrainRefSeqKey.Strain{organismIndex};
    else
        warning(['Could not find a strain name for RefSeq accession ' accession '. Skipping...'])
        notFound = [notFound;i];
    end
end
modelFiles(notFound) = [];
organismIDs(notFound) = [];

%% Load files
disp('Loading files...')

addpath(genpath(pathMatTFA));
addpath(genpath(fullfile(ngDir,'gapfilling')));

disp('    Loading medium file...')
load(mediaFile)
minMed = replace(minMed,'[e]','_e');
nutrients = replace(nutrients,'[e]','_e');
vitamins = replace(vitamins,'[e]','_e');

disp('    Loading experimental data file...')
if ~isfile(growthDataFile)
    error(['Binarized growth data not found: ' growthDataFile]);
end
[~, ~, raw] = xlsread(growthDataFile);
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

% If retryMode = 'continuous', growthRateFile must be the same layout as
% growthDataFile: row 1 = carbon source names, column 1 = strain name,
% cells = continuous growth rate instead of binary growth/no-growth.
growthRateData = [];
if strcmp(retryMode,'continuous')
    if isfile(growthRateFile)
        [~, ~, rawRate] = xlsread(growthRateFile);
        if ~isequal(rawRate(1,2:end),rawHeader)
            error('growthRateFile''s carbon source columns must match growthDataFile exactly, same names, same order.');
        end
        growthRateData.strainNames = rawRate(2:end,1);
        growthRateMatrix = rawRate(2:end,2:end);
        growthRateMatrix = reshape([growthRateMatrix{:}],size(growthRateMatrix));
        growthRateMatrix(:,waterIndex) = []; % Keep columns aligned with nutrients, which has water removed
        growthRateData.growthRate = growthRateMatrix;
    else
        fprintf('retryMode is ''continuous'' but growthRateFile was not found. Falling back to binary.\n');
    end
end

disp('    Loading model database file...')
load(modelDBFile)
modelDB = model;
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
global TFA_MILP_SOLVER; TFA_MILP_SOLVER = 'gurobi_direct';

%% Enumerate alternative gap-filling reaction sets for each model
for MF = 1:length(modelFiles)

tic

modelFile = fullfile(modelDir,[modelFiles{MF} '.mat']);
organismID = organismIDs{MF};

fprintf(['\n-------------------------------------------------------------------------------------\n\nModel: ' organismID '\n\n'])

saveFileNameGFSoln = fullfile(saveDirGFSoln,['AltRxns_' organismID]);

growthNutrients = find(growthData(ismember(strainNames,organismID),:));

% Build the ranked list of carbon sources to retry with if gapfilling
% fails. Uses continuous growth rate if growthRateFile was given and has
% this organism, otherwise falls back to the binary true-positive list.
retryNutrients = nutrients(growthNutrients);
if ~isempty(growthRateData)
    strainIdx = find(strcmp(growthRateData.strainNames,organismID));
    if ~isempty(strainIdx)
        strainRates = growthRateData.growthRate(strainIdx,:);
        [sortedRates,rateOrder] = sort(strainRates,'descend');
        retryNutrients = nutrients(rateOrder(sortedRates > 0));
    end
end

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

% Set model
[modelOrig,sourceModel] = deal(model);

% Constrain model exchange reactions (but keep added exchange reactions open, as well as methanol)
sourceModel.lb(excRxns) = 0;
sourceModel.lb(intersect(find(contains(model.rxns,'EX_')),find(ismember(sourceModel.rxns,findRxnsFromMets(sourceModel,[minMed;vitamins]))))) = -50;

% Force uptake of oxygen
if forceO2Uptake
    sourceModel.ub(intersect(excRxns,find(ismember(sourceModel.rxns,findRxnsFromMets(sourceModel,'o2_e'))))) = -1;
end

alternateTruePositives = {};

% Merge model with database
fprintf('\n    Preparing for gapfilling...\n')
try [GFModelMaster, conflict] = PrepareForGapFilling(sourceModel,{modelDB},'',0,flagTFA,{},[],thermoData);
catch % If there is no feasible solution, add a true positive carbon source into the medium
    fails = 1;
    tp = 0;
    alternateTruePositives = [alternateTruePositives,setdiff(retryNutrients,alternateTruePositives)];
    alternateTruePositives(find(ismember(alternateTruePositives,setdiff(alternateTruePositives,retryNutrients)))) = [];
    while fails == 1
        tp = tp+1;
        if tp > length(alternateTruePositives) % Exit if no alternate carbon sources are left to try
            disp('Model could not be merged with the database. Skipping.')
            break
        end
        truePos = alternateTruePositives{tp};
        sourceModel.lb(excRxns) = 0;
        sourceModel.lb(intersect(excRxns,find(ismember(sourceModel.rxns,findRxnsFromMets(sourceModel,[minMed;vitamins;truePos]))))) = -50;
        disp(['    Re-trying with ' truePos ' in medium...'])
        try
            [GFModelMaster, conflict] = PrepareForGapFilling(sourceModel,{modelDB},'',0,flagTFA,{},[],thermoData);
            fails = 0;
        catch
            fails = 1;
        end
    end
    if fails == 1; continue; end
end

%% Gap-fill across all carbon sources
fprintf('\nGapfilling...')
[ActRxnsAll,foundSolution] = performGapfillingTFA(nutrients,GFModelMaster,excRxns,minMed,vitamins,modelOrig,maxNumAlt,maxSolveTime);

% If no solution was found, add a true positive carbon source into the medium and try gapfilling again
if foundSolution == 0
    sourceModelOrig = sourceModel; % Save sourceModel since exchange reactions will be opened in next step
    alternateTruePositives = [alternateTruePositives,setdiff(retryNutrients,alternateTruePositives)];
    alternateTruePositives(find(ismember(alternateTruePositives,setdiff(alternateTruePositives,retryNutrients)))) = [];
    tp = 0;
    while foundSolution == 0
        tp = tp+1;
        if tp > length(alternateTruePositives) % Exit if no gapfilling solutions were found after cycling through all alternate carbon sources
            disp('Model gapfilling failed.')
            break
        end
        truePos = alternateTruePositives{tp};
        disp(['No gapfilling solutions found. Re-trying with ' truePos ' in medium...'])
        sourceModel = sourceModelOrig;
        sourceModel.lb(find(ismember(sourceModel.rxns,intersect(sourceModel.rxns(excRxns),findRxnsFromMets(sourceModel,truePos))))) = -50;

        try
            [GFModelMaster, conflict] = PrepareForGapFilling(sourceModel,{modelDB},'',0,flagTFA,{},[],thermoData);
            [ActRxnsAll,foundSolution] = performGapfillingTFA(nutrients,GFModelMaster,excRxns,minMed,vitamins,modelOrig,maxNumAlt,maxSolveTime);
        catch
            % This carbon source also fails to merge with the database; try the next one.
        end
    end
end

%% Save gapfilling solutions
fprintf(['\nGapfilling for ' organismID ' complete. Saving results...'])
save(saveFileNameGFSoln,'ActRxnsAll','nutrients','-v7.3')
fprintf('\nDone.\n')

t = toc;
fprintf(['\nModel ' organismID ' complete. Time elapsed: ' num2str(round(t/60,2)) ' minutes.\n'])

end

%% Remove path
rmpath(genpath(pathMatTFA));

%% Main gapfilling function
function [ActRxnsAll,foundSolution] = performGapfillingTFA(nts,gfmodel,exrxns,mm,vits,origmodel,maxalt,maxtime)
    foundSolution = 0; % Toggled when at least one solution is found that leads to model growth
    for i = 1:length(nts)

        if i > 1
            fprintf(repmat('\b',1,length(todisp)-1));
        end

        todisp = ['\n    Carbon source ' num2str(i) ' of ' num2str(length(nts)) ': ' nts{i} ' '];
        fprintf(todisp);

        % Define the media for the model
        GFmodel = gfmodel;
        media = GFmodel.rxns(intersect(exrxns,find(ismember(GFmodel.rxns,findRxnsFromMets(GFmodel,[mm;vits;nts{i}])))));
        f = find(ismember(GFmodel.varNames,strcat('R_',media)));
        GFmodel.var_ub(f) = 50;
        GFmodel.var_lb(f) = -50;

        % Apply a lower bound on growth to gapfill for biomass production
        GFmodel.var_lb(find(ismember(GFmodel.varNames,strcat('F_',origmodel.rxns(find(origmodel.c)))))) = 0.01;

        % Gapfill, generating the alternative solutions. A single carbon
        % source failing here should not lose the other 45.
        try
            [resultStat,ActRxns,DPsAll] = gfbiomass(GFmodel,GFmodel.indUSE,maxalt,maxtime,1,0,0,'');
        catch
            ActRxns = cell(1,1); % Same shape gfbiomass itself uses for "cannot be gap-filled"
        end

        % Record alternative solutions
        ActRxnsAll.(matlab.lang.makeValidName(nts{i})) = ActRxns;

        % Determine if a solution was found
        if ~isempty(ActRxns{1,1})
            foundSolution = 1;
        else
            todisp = ['\n']; % Don't erase the warning message
        end
    end
    fprintf(repmat('\b',1,length(todisp)-1));
    fprintf('\n')
end
