function [cost] = CostCombined(params,model, expInd, data, stimLipo, stimAdi, limit, doHSL)

if (nargin < 8), doHSL = 1; end

if iscolumn(params); params=params'; end

params(expInd)=exp(params(expInd));
nParams=find(strcmp(IQMparameters(model),'pip'))-3;

if length(params)==nParams+2 %If the parameter vector contains two more (diabetes + diab_reest) parameter than the free parameters
    diabetes=params(end-1:end);
    params(end-1:end)=[];
elseif length(params)==nParams+1 %If the parameter vector contains one more (diabetes) parameter than the free parameters
    diabetes=[params(end) 1];
    params(end)=[];
else
    diabetes=[1 0.33]; % not used
end

costNyman = CostNyman([params diabetes], model, [], data.Rajan, limit);
costInVivo = CostInVivo(params, model, data.InVivo, 0 );%% InVivo data
costInVitro = CostInVitro(model, params, stimLipo, data,diabetes, doHSL);%% InVitro data
costAdi = CostAdi(params, [], model, data.Adi, stimAdi, limit);

cost=costNyman+costInVivo+costInVitro+costAdi;

end