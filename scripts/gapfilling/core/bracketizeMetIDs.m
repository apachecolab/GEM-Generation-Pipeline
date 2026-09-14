function model = bracketizeMetIDs(model)
% Converts a model's metabolite IDs from this pipeline's working convention
% (compartment suffix '_c'/'_e'/'_p') to bracket notation ('atp[c]'),
% matching the convention gapfilling_selection.m saves its own gap-filled
% models in. Applied on save so every file under data/models/ uses the
% same convention regardless of which script produced it.
%
% Inputs:
%   model - COBRA model with underscore-suffixed metabolite IDs
%
% Outputs:
%   model - same model, metabolite IDs in bracket notation
%
% Author: Sudharshan Ravi

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
end
