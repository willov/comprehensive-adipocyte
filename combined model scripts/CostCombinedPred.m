function [p, cost] = CostCombinedPred(params,model, expInd, data, stimLipo, stimAdi, limit, doHSL, stimAdiVal, toEstimate, pOpt)


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

[costNyman, yGlu] = CostNyman([params diabetes], model, [], data.Rajan, limit);
[costInVivo, yInVivo] = CostInVivo(params, model, data.InVivo, 0);%% InVivo data
[costInVitro,yInVitro,costHSL] = CostInVitro(model, params, stimLipo, data,diabetes, doHSL);%% InVitro data
[costAdi, costCL_ATP,yAdi] = CostAdi(params, [], model, data.Adi, stimAdi, limit, stimAdiVal);

cost=costNyman+costInVivo+costInVitro+costAdi;
if cost<1e20
    %% Select the right prediction to estimate the uncertainty of.
    if strcmp(toEstimate.submodel, 'lipoadi')
        p = costHSL+costCL_ATP;
    elseif strcmp(toEstimate.submodel, 'glu')
        tIdx=yGlu.(toEstimate.experiment){"Point",toEstimate.variable}==str2double(toEstimate.point);
        p=yGlu.(toEstimate.experiment){"Sim",toEstimate.variable}(tIdx);
    elseif strcmp(toEstimate.submodel, 'lipo')
        if contains(toEstimate.experiment, 'Fig')
            tIdx=yInVivo.(toEstimate.variable).Time==str2double(toEstimate.point);
            p = yInVivo.(toEstimate.variable).(toEstimate.experiment)(tIdx);
        else
            x=str2double(toEstimate.point);
            if x==0
                doseInd=yInVitro.(toEstimate.experiment).Ins==x & yInVitro.(toEstimate.experiment).Iso==0;
            else
                doseInd=yInVitro.(toEstimate.experiment).Ins==x;
            end
            p = yInVitro.(toEstimate.experiment).(toEstimate.variable)(doseInd);
        end
    elseif strcmp(toEstimate.submodel, 'adi')
        tIdx=yAdi{"Time", "Measures"}==str2double(toEstimate.point);
        p=yAdi{toEstimate.experiment, "Measures"}(tIdx);
    elseif strcmp(toEstimate.submodel, 'param')
        paramIdx = strcmp(IQMparameters(model), toEstimate.variable);
        params=[params diabetes];
        p = params(paramIdx);
    end

    %% Set the value based on searched direction (min or max)
    if strcmp(toEstimate.polarity, 'max')
        polarity=-1;
    else
        polarity=1;
    end
    p = p*polarity;

    if cost> limit
        p=p+abs(p)+abs(pOpt)*(1+(cost-limit));
    end

else
    p=cost;
end