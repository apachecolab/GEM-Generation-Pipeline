% Precomputes reaction-level essentiality (findEssentialReactions.m) for
% every model in donorDir, against every carbon source in the pipeline's
% panel. Required before running gapfilling_curation.m with
% donorReactionMode = 'essential'; not needed for donorReactionMode = 'fba'.
%
% Each donor's result is written to disk as soon as it finishes and is
% skipped on a rerun, so an interrupted sweep can just be restarted.
%
% Inputs:
%   <donorDir>/*.mat                                   - external donor COBRA models (variable name 'model')
%   data/experimental/medium/min_med_CSourceScreen.mat - minMed, nutrients, vitamins
%   scripts/gapfilling/core/                           - growModelInCSources, findEssentialReactions, normalizeMetIDs
%
% Outputs:
%   <donorDir's parent>/essential/<donorModelID>.mat
%       - results: 1 x numCS struct array (carbon_source, mu_star, rxns, essential)
%
% Author: Sudharshan Ravi

clear all; clc;

%% Parameters
donorDir = fullfile('data','donors','fba'); % required: directory of external donor .mat files to sweep, relative to the repository root
nWorkers = 4; % used only if Parallel Computing Toolbox is available
growthThreshold = 5e-3; % essentiality only computed when mu_star exceeds this
essentialThreshold = 1e-6;
vmax = 10;

if isempty(donorDir); error('Set donorDir to the external donor directory to sweep.'); end

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
rootDir = fullfile(scriptDir,'..','..','..');
donorDir = fullfile(rootDir,donorDir);
dataDir = fullfile(rootDir,'data');
coreDir = fullfile(rootDir,'scripts','gapfilling','core');
addpath(coreDir);

mediaFile = fullfile(dataDir,'experimental','medium','min_med_CSourceScreen.mat');

donorDirClean = regexprep(donorDir,[filesep '+$'],'');
essentialDir = fullfile(fileparts(donorDirClean),'essential',filesep); % sibling of donorDir
if ~exist(essentialDir,'dir'); mkdir(essentialDir); end

%% Load medium
load(mediaFile)
minMed = replace(minMed,'[e]','_e');
nutrients = replace(nutrients,'[e]','_e');
vitamins = replace(vitamins,'[e]','_e');
numCS = length(nutrients);

%% Check if COBRA toolbox is loaded
if ~exist('optimizeCbModel','file'); initCobraToolbox(false,'agent'); end
changeCobraSolver('gurobi');

%% Discover donor models, skip any already swept
donorFiles = dir(fullfile(donorDir,'*.mat'));
donorIDs = erase({donorFiles.name},'.mat');
alreadyDone = cellfun(@(id) exist(fullfile(essentialDir,[id '.mat']),'file') > 0, donorIDs);
pending = find(~alreadyDone);
fprintf('%d/%d donor(s) already swept. %d remaining.\n', sum(alreadyDone), length(donorIDs), length(pending));

%% Parallel pool, only if the Parallel Computing Toolbox is licensed and installed
useParallel = license('test','Distrib_Computing_Toolbox') && ~isempty(ver('parallel'));
if useParallel
    pool = gcp('nocreate');
    if isempty(pool) || pool.NumWorkers ~= nWorkers
        if ~isempty(pool); delete(pool); end
        parpool(nWorkers);
    end
else
    disp('Parallel Computing Toolbox not available. Running sequentially.')
end

%% Sweep
parfor pi = 1:length(pending)
    di = pending(pi);
    donorID = donorIDs{di};

    % parfor workers are separate processes and don't share the client
    % session's globals -- each one needs its own solver configured before
    % it can call optimizeCbModel, regardless of whether the function is
    % already visible on its path.
    initCobraToolbox(false,'agent');
    changeCobraSolver('gurobi');

    d = load(fullfile(donorDir,donorFiles(di).name));
    donorModel = normalizeMetIDs(d.model);

    results = repmat(struct('carbon_source','','mu_star',NaN,'rxns',{{}},'essential',[]),1,numCS);

    [~,growthRates] = growModelInCSources(donorModel,minMed,vitamins,[],nutrients,vmax);

    for ci = 1:numCS
        results(ci).carbon_source = nutrients{ci};
        results(ci).mu_star = growthRates(ci);

        if growthRates(ci) > growthThreshold
            model = donorModel;
            model.lb(find(findExcRxns(model))) = 0;
            model.lb(find(ismember(model.rxns,intersect(model.rxns(find(findExcRxns(model))),findRxnsFromMets(model,minMed))))) = -1000;
            for vi = 1:length(vitamins)
                vitRxns = intersect(model.rxns(find(findExcRxns(model))),findRxnsFromMets(model,vitamins(vi)));
                if ~isempty(vitRxns); model.lb(ismember(model.rxns,vitRxns)) = -1000; end
            end
            model.lb(find(ismember(model.rxns,intersect(model.rxns(find(findExcRxns(model))),findRxnsFromMets(model,nutrients{ci}))))) = -vmax;

            results(ci).rxns = model.rxns;
            results(ci).essential = findEssentialReactions(model, essentialThreshold);
        end
    end

    parsaveEssentialResults(fullfile(essentialDir,[donorID '.mat']), results);
    fprintf('%s done (%d/%d).\n', donorID, pi, length(pending));
end

disp('Essentiality sweep complete.')

function parsaveEssentialResults(outFile, results)
    save(outFile,'results','-v7.3')
end
