function [model,data, lb, ub, nParams,expInd, stimLipo, stimAdi, stimAdiVal, dgf, opts] = Init_combined(modelName, doHSL, useCL_ATP)
if nargin<1; modelName='lipolysis'; end
if nargin<2, doHSL=1; end
if nargin<3, useCL_ATP=1; end
format short e

Setup_combined()

%% Setup lipolysis
[model,data, lb, ub, nParams,expInd, stimLipo, dgfLipo, opts] = Init_lipo(modelName, doHSL);

% Setup glucose
lb = [lb 0];
ub = [ub 1];

[~, expData, ~, ~, ~, ~, dgfRajan, ~] = Init_glu(modelName);
data.Rajan = expData;

%% Setup adiponectin
if useCL_ATP
    toEstimateOn={'Ca_ATP','Ca_noATP','noCa_ATP','noCa_noATP','highCa_ATP','EPI_ATP','CL_Ca', 'CL_ATP'};
    toValidateOn={};
else
    toEstimateOn={'Ca_ATP','Ca_noATP','noCa_ATP','noCa_noATP','highCa_ATP','EPI_ATP','CL_Ca'};
    toValidateOn={'CL_ATP'};
end
[~, dataAdi, stimAdi, stimAdiVal, dgfAdi]=Init_adi(modelName, toEstimateOn, toValidateOn);

data.Adi = dataAdi;
dgf = dgfLipo+dgfAdi+dgfRajan;

end