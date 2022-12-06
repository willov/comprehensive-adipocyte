function [simdata_60_n, simdata_60_d]=CoreSimulation(coreModelName, inhibition)
%% Load data
% close all
if nargin<2, inhibition=1; end

useCL_ATP=1;
doDiabetes=1;

warning off
[coreModel,~, ~, ~, ~, expInd] = Init_combined(coreModelName, doDiabetes, useCL_ATP);
warning on

load('../combined model scripts/Results/With HSL and CL_ATP/combined, opt(654.812129) 220214-162036.mat', 'optParam')
optParam(expInd)=exp(optParam(expInd));

pNames = IQMparameters(coreModel);

%% Defining parameter values
inhibIdx = ismember(pNames, {'kG4a','kG4e'});
optParam(inhibIdx)=optParam(inhibIdx).*inhibition;

sim = SimulateGlucose(optParam, coreModel);
simdata_60_n = sim.n_60;
simdata_60_d = sim.d_60;
end