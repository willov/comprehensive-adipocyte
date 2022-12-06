function [model, expData, estimation, validation, dgf, pNames, expInd, nParams, lb, ub, opts]=Init_adi(modelName, toEstimateOn, toValidateOn)
%% Loading data-sets, to be used in objective function/plotting

load('expDataAdiponectin','expData')
%% Model
design=matlab.lang.makeValidName({'Ca+ATP';'Ca+noATP';'noCa+ATP';'noCa+noATP';'highCa+ATP';...
    'EPI+ATP';'CL+ATP';'CL+Ca'});
%   Ca,           ATP,         cAMP,   Epi, CL
stimulus=[...
    0.0015         3           0.1     0    0;   %'Ca+ATP'
    0.0015         0           0.1     0    0;   %'Ca+noATP'
    0              3           0.1     0    0;   %'noCa+ATP'
    0              0           0.1     0    0;   %'noCa+noATP'
    0.015          3           0.1     0    0;   %'highCa+ATP'
    0              3           0       5    0;   %'EPI+ATP' 
    0              3           0       0    1;   %'CL+ATP'  
    0.0015         3           0       0    1];   %'CL+Ca'  

experimentalSetup=table(stimulus(:,1:3), stimulus(:,4:end), 'VariableNames',{'Pipette','Agonist'},'RowNames',design);

if nargin<1 || isempty(modelName)
    modelName='adr_endo'; % Name of model file.
end
if nargin<2
toEstimateOn={'Ca_ATP','Ca_noATP','noCa_ATP','noCa_noATP','highCa_ATP','EPI_ATP','CL_Ca'};
end
scaleExp={'Ca_ATP','Ca_noATP','noCa_ATP','noCa_noATP'};

if nargin<3
toValidateOn={'CL_ATP'};
end

expData(:,~ismember(expData.Properties.VariableNames,[toEstimateOn toValidateOn scaleExp(~ismember(scaleExp,toEstimateOn))]))=[];
estimation=experimentalSetup(toEstimateOn,:);
validation=experimentalSetup(toValidateOn,:);

model=str2func(modelName); % Sets the model as a function.

pNames=IQMparameters(model);
nParams = find(strcmp(pNames,'pip'))-3; %-3 to ignore pip and diabetes parameters

dgf=numel(expData{'Mean',toEstimateOn})-1;
dgf(2)=numel(expData{'Mean',toValidateOn});
if dgf(1)<0, dgf(1) =numel(expData{'Mean',toEstimateOn})-1; end 
ub=1e6*ones(1,nParams-2); % removes VApip and VAcell and basal cAMp/ATP/Ca2+. Sets the upper bounds of the parameter values. 
lb=1e-6*ones(1,nParams-2); % Sets the lower bounds.

[lb,ub] = CustomBounds(pNames, lb, ub);
expInd=~cellfun(@isempty,(regexp(pNames(1:nParams),'^k.+'))); 

lb(expInd)=log(lb(expInd));
ub(expInd)=log(ub(expInd));

if contains(modelName,'_hill')
    hillInd=~cellfun(@isempty,regexp(pNames,'^n[0-9]+$')); %Finds all names that begins with n, followed with numbers and then ends.
    ub(hillInd)=log(3);
    lb(hillInd)=log(1);
end

%% Setup optimization
opts.ndiverse     = 50; %'auto'; %100; %500; %5; %
opts.maxtime      = 750; % In cess this option will be overwritten
opts.maxeval      = 1e4;
opts.log_var      = [];

opts.local.solver = 'dhc'; %'dhc'; %'fmincon'; %'nl2sol'; %'mix'; %
opts.local.finish = opts.local.solver;
opts.local.bestx = 0;
opts.local.tol = 1;
opts.local.iterprint = 1;

opts.dim_refset   = 'auto'; %

end