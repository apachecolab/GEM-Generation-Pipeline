function [pG,growthRatesAlt] = growModelInCSources(mdl,mm,vits,vitaminBounds,nts,vmax)

    if isempty(vitaminBounds); vitaminBounds = repmat(-1000,length(vits),1); end

    [growthRatesAlt,pG] = deal(zeros(1,length(nts)));
    for j = 1:length(nts)

        modelTest = mdl;

        % Define the model's medium
        modelTest.lb(find(findExcRxns(modelTest))) = 0;
        modelTest.lb(find(ismember(modelTest.rxns,intersect(modelTest.rxns(find(findExcRxns(modelTest))),findRxnsFromMets(modelTest,mm))))) = -1000;
        for v = 1:length(vits)
            vitRxns = intersect(modelTest.rxns(find(findExcRxns(modelTest))),findRxnsFromMets(modelTest,vits(v)));
            if ~isempty(vitRxns)
                modelTest.lb(find(ismember(modelTest.rxns,vitRxns))) = vitaminBounds(v);
            end
        end
        modelTest.lb(find(ismember(modelTest.rxns,intersect(modelTest.rxns(find(findExcRxns(modelTest))),findRxnsFromMets(modelTest,nts{j}))))) = -vmax;

        FBAsoln = optimizeCbModel(modelTest);
        growthRatesAlt(j) = FBAsoln.f;
    end
    pG(find(growthRatesAlt > 1e-3)) = 1;
end