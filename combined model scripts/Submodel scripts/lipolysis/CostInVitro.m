function [cost, responses, costHSL]=CostInVitro(model, param, stimulus, data, diab, doHSL)
if nargin <6, doHSL=0; end
try
    
    [responses, steady, ~, basal]=simulateInVitro(model, param, diab, stimulus);
    
    if any(abs(basal{:,:})<1e-3)
        penalty=sum(1./abs(basal{:,:}))+50;
    else
        penalty = 0;
    end
    
    residuals=(responses.Normal{:,'Glycerol'}-data.InVitro.Glycerol.Mean).^2./data.InVitro.Glycerol.SEM.^2;
    costGlycerol=sum(residuals(~any(isnan(data.InVitro.Glycerol{:,3:4}),2))); %only som residuals where

    residuals=(responses.Normal{:,'FA'}-data.InVitro.FA.Mean).^2./data.InVitro.FA.SEM.^2;
    costFA=sum(residuals(~any(isnan(data.InVitro.FA{:,3:4}),2))); %only som residuals where

    ind=ismember(responses.Normal{:,{'Ins'}},data.InVitro.PKB473{:,{'Ins'}},'rows');
    residuals=(responses.Normal{ind,'PKB473'}-data.InVitro.PKB473.Mean).^2./data.InVitro.PKB473.SEM.^2;
    costPKB=sum(residuals(~any(isnan(data.InVitro.PKB473{:,3:4}),2))); %only som residuals where%
  
    ind=ismember(responses.Normal{:,{'Ins'}},data.InVitro_diabetes.Reesterification{:,{'Ins'}},'rows');
    residuals=(responses.Diabetes{ind,'Reesterification'}-data.InVitro_diabetes.Reesterification.Mean).^2./data.InVitro_diabetes.Reesterification.SEM.^2;
    costReesterificationDiab=sum(residuals(~any(isnan(data.InVitro_diabetes.Reesterification{:,3:4}),2))); %only som re

    residuals=(responses.Normal{:,'HSL'}-data.InVitro.HSL.Mean).^2./data.InVitro.HSL.SEM.^2;
    costHSL=sum(residuals(~any(isnan(data.InVitro.HSL{:,3:4}),2))); %only som residuals where%
    
    cost = costGlycerol+costFA+costPKB;
    cost = cost + costReesterificationDiab;

    if doHSL
        cost=cost+costHSL;
    end

    cost=cost+penalty+steady;

catch err
    cost=1e24; %returns a "large" error.
    costHSL=inf;
    responses=[];
end



end

