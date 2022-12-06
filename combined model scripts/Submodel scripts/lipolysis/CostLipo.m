function [cost] = CostLipo(params,model, expInd, data, stimulus, doHSL)

if (nargin < 6), doHSL = 0; end

if iscolumn(params); params=params'; end

params(expInd)=exp(params(expInd));
nParams=find(strcmp(IQMparameters(model),'pip'))-3;

if length(params)==nParams+2 %If the parameter vector contains two more (diabetes + diab_reest) parameter than the free parameters
    diabetes=params(end-1:end);
    params(end-1:end)=[];
elseif length(params)==nParams+1 %If the parameter vector contains one more (diab_reest) parameter than the free parameters
    diabetes=[1 params(end)];
    params(end)=[];
else
    diabetes=[1 0.33]; % not used
end

costInVivo = CostInVivo(params, model, data.InVivo, 0);%% InVivo data
costInVitro = CostInVitro(model, params, stimulus, data,diabetes, doHSL);%% InVitro data

cost=costInVivo+costInVitro;

end
