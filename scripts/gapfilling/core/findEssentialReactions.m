function essential = findEssentialReactions(model, threshold)
% Reaction-level essentiality under the model's current bounds (medium and
% carbon source already applied by the caller). For every reaction, force
% its bounds to [0,0] and re-solve FBA; a reaction is essential if growth
% collapses below threshold * wild-type growth without it. Operates purely
% on reaction bounds, independent of any gene-reaction mapping.
%
% Inputs:
%   model     - COBRA model with medium/carbon-source bounds already applied
%   threshold - essential if grRateKO/grRateWT < threshold (default 1e-6)
%
% Outputs:
%   essential - n_rxns x 1 logical vector, true where that reaction is
%               essential for growth under the model's current bounds
%
% Author: Sudharshan Ravi

if nargin < 2 || isempty(threshold); threshold = 1e-6; end

solWT = optimizeCbModel(model, 'max');
if solWT.stat ~= 1 || solWT.f <= 0
    error('findEssentialReactions: wild-type model does not grow (stat = %d, f = %.4g).', solWT.stat, solWT.f);
end
grRateWT = solWT.f;

n = length(model.rxns);
essential = false(n,1);

for i = 1:n
    modelDel = changeRxnBounds(model, model.rxns{i}, 0, 'b');
    solKO = optimizeCbModel(modelDel, 'max');
    if solKO.stat == 1
        grRateKO = solKO.f;
    else
        grRateKO = 0;
    end
    essential(i) = (grRateKO / grRateWT) < threshold;
end
end
