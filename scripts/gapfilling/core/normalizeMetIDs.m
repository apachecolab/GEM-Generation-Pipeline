function model = normalizeMetIDs(model)
% Normalizes a model's metabolite IDs to this project's working convention
% (compartment suffix '_c'/'_e'/'_p'). Detects bracket convention (e.g.
% 'atp[c]') and converts it; leaves an already-underscore-suffixed model
% (e.g. 'atp_c') untouched. Needed both for external donor models and for
% this pipeline's own gap-filled/curated models, which gapfilling_selection.m
% saves in bracket notation despite working in underscore notation
% internally. Also derives .rev from bounds when the field is absent, which
% many external COBRA exports omit.
%
% Inputs:
%   model - COBRA model as loaded from disk
%
% Outputs:
%   model - same model, met IDs normalized, .rev guaranteed present
%
% Author: Sudharshan Ravi

if any(endsWith(model.mets,']'))
    model.mets = regexprep(model.mets,'\[([cep])\]$','_$1');
end

if ~isfield(model,'rev')
    model.rev = zeros(length(model.rxns),1);
    model.rev(intersect(find(model.lb < 0), find(model.ub > 0))) = 1;
end
end
