function [model,data, lb, ub, nParams,expInd, stimulus, dgf, opts] = Init_lipo(modelName, doHSL)
if nargin<1; modelName='lipolysis'; end
if nargin<2, doHSL=0; end

format short e

model=str2func(modelName);

%%  Setup in vivo data
load('stichfig1.mat', 'stichfig1'); stichfig1.ins = stichfig1.ins * 6.945*1e-3; % converting from mU/L to nM
load('stichfig2.mat', 'stichfig2'); stichfig2.ins = stichfig2.ins * 6.945*1e-3;
load('stichfig3.mat', 'stichfig3'); stichfig3.ins = stichfig3.ins * 6.945*1e-3;

data.InVivo.Fig1.Mean=stichfig1.val(1:30);
data.InVivo.Fig1.SEM=stichfig1.se(1:30);
data.InVivo.Fig1.Time=stichfig1.time(1:30);
data.InVivo.Fig1.Ins=stichfig1.ins;

data.InVivo.Fig2Epi.Mean=stichfig2.epi(1:30);
data.InVivo.Fig2Epi.SEM=stichfig2.epi_se(1:30);
data.InVivo.Fig2Epi.Time=stichfig2.time(1:30);
data.InVivo.Fig2Epi.Ins=stichfig2.ins(1:30);

data.InVivo.Fig2Iso.Mean=stichfig2.iso(1:30);
data.InVivo.Fig2Iso.SEM=stichfig2.iso_se(1:30);
data.InVivo.Fig2Iso.Time=stichfig2.time(1:30);
data.InVivo.Fig2Iso.Ins=stichfig2.ins(1:30);

data.InVivo.Fig3Epi.Mean=stichfig3.epi(1:30);
data.InVivo.Fig3Epi.SEM=stichfig3.epi_se(1:30);
data.InVivo.Fig3Epi.Time=stichfig3.time(1:30);
data.InVivo.Fig3Epi.Ins=stichfig3.ins(1:30);

data.InVivo.Fig3EpiPhe.Mean=stichfig3.epiphe(1:30);
data.InVivo.Fig3EpiPhe.SEM=stichfig3.epiphe_se(1:30);
data.InVivo.Fig3EpiPhe.Time=stichfig3.time(1:30);
data.InVivo.Fig3EpiPhe.Ins=stichfig3.ins(1:30);

%% Setup in vitro data
load('expDataLipolysis','expData')
expData.FA.SEM(end-2:end-1)= mean(expData.FA.SEM([2:end-3 end]));

expData.FA.SEM(expData.FA.SEM==0)=nan;
expData.Glycerol.SEM(expData.Glycerol.SEM==0)=nan;
expData.HSL.SEM(expData.HSL.SEM==0)=nan;
expData.PKB473.SEM(expData.PKB473.SEM==0)=nan;

data.InVitro.FA=expData.FA;
data.InVitro.Glycerol=expData.Glycerol;
data.InVitro.HSL=expData.HSL;
data.InVitro.PKB473=expData.PKB473;

%% Setup diabetes data
load('expDataDiabetes','expDataDiabetes')
data.InVitro_diabetes.FA=expDataDiabetes.FA;
data.InVitro_diabetes.Glycerol=expDataDiabetes.Glycerol;
data.InVitro_diabetes.HSL=expDataDiabetes.HSL;
data.InVitro_diabetes.Reesterification=expDataDiabetes.reesterificationDiabetes;
data.InVitro.Reesterification=expDataDiabetes.reesterification;

%% setup degrees of freedom
dgfGlycerol=sum(~isnan(data.InVitro.Glycerol.SEM));
dgfFA=sum(~isnan(data.InVitro.FA.SEM));
dgfHSL=sum(~isnan(data.InVitro.HSL.SEM));
dgfPKB=sum(~isnan(data.InVitro.PKB473.SEM));
dgfReesterificationDiab=sum(~isnan(data.InVitro_diabetes.Reesterification.SEM));

dgfFig1Epi=sum(~isnan(data.InVivo.Fig1.SEM));
dgfFig2Epi=sum(~isnan(data.InVivo.Fig2Epi.SEM));
dgfFig2Iso=sum(~isnan(data.InVivo.Fig2Iso.SEM));
dgfFig3Phe=sum(~isnan(data.InVivo.Fig3EpiPhe.SEM));

dgfInVitro=dgfGlycerol+dgfFA+dgfPKB+dgfReesterificationDiab;
dgfInVivo=dgfFig1Epi+dgfFig2Iso+dgfFig3Phe+dgfFig2Epi;

dgf=dgfInVivo+dgfInVitro-1-6; %-1 for drift parameter, -6 for scaling

if doHSL
    dgf=dgf+dgfHSL;
    dgf(2)=0;
else
    dgf(2)=dgfHSL;
end

%% Parameter bound setup
[pNames, ~] = IQMparameters(model);
nParams=find(strcmp(pNames,'pip'))-1;
expInd=~cellfun(@isempty,(regexp(pNames(1:nParams-2),'^k.+'))); % -2 for ignoring the diabetes parameters

lb=repmat(1e-8, 1, nParams-2);
ub=repmat(1e6, 1, nParams-2);

[lb,ub] = CustomBounds(pNames, lb, ub);

lb       = [lb 0]; % Add one extra value for the effect of diabetes
ub       = [ub 1]; % Add one extra value for the effect of diabetes

lb(expInd) = log(lb((expInd)));
ub(expInd) = log(ub((expInd)));

%% setup stimulus
stimulus=table();
ins=unique([expData.FA.Ins; expData.Glycerol.Ins]);
stimulus.Ins=[ins; 0];
stimulus.Iso=[0.01*ones(height(stimulus)-1,1); 0]; % 10nm = 0.01 µM

%% Setup optimization
opts.ndiverse     = 500; %'auto'; %100; %500; %5; %
opts.maxtime      = 750; % In cess this option will be overwritten
opts.maxeval      = 1e7;
opts.log_var      = [];

opts.local.solver = 'dhc'; %'dhc'; %'fmincon'; %'nl2sol'; %'mix'; %
opts.local.finish = opts.local.solver;
opts.local.bestx = 0;
opts.local.tol = 2;
opts.local.iterprint = 1;

opts.dim_refset   = 'auto'; %

if(strcmp(opts.local.solver,'fmincon'))
    opts.local.use_gradient_for_finish = 1; %DW: provide gradient to fmincon
else
    opts.local.use_gradient_for_finish = 0; %DW: provide gradient to fmincon
end
opts.local.check_gradient_for_finish = 0; %DW: gradient checker

end

