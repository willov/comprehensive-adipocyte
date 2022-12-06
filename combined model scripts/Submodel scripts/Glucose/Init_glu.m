function [model,expData, lb, ub, nParams,expInd, dgf, opts] = Init_glu(modelName)
format short e

model=str2func(modelName);
variables = IQMvariables(model);

%%  Setup data
load('expDataRajan.mat','expData'); 

% Average the model uncertainty
expData.IR_dr.sem_n(2:end)=mean(expData.IR_dr.sem_n(2:end));
expData.IR_dr.sem_d(2:end)=mean(expData.IR_dr.sem_d(2:end));
expData.IRS1_dr.sem_n(2:6) = mean(expData.IRS1_dr.sem_n(2:5));
expData.IRS1_dr.sem_d(2:6)=mean(expData.IRS1_dr.sem_d(2:6));
expData.IRS1307_dr.sem_n(2:6) = mean(expData.IRS1307_dr.sem_n(2:5));
expData.IRS1307_dr.sem_d(2:6)=mean(expData.IRS1307_dr.sem_d(2:6));
expData.PKB308_dr.sem_n(2:7)=mean(expData.PKB308_dr.sem_n(2:7));
expData.PKB308_dr.sem_d(2:7)=mean(expData.PKB308_dr.sem_d(2:7));
expData.PKB473_dr.sem_n(2:7)=mean(expData.PKB473_dr.sem_n(2:7));
expData.PKB473_dr.sem_d(2:7)=mean(expData.PKB473_dr.sem_d(2:7));
expData.AS160_dr.sem_n(2:6) = mean(expData.AS160_dr.sem_n(2:6));
expData.AS160_dr.sem_d(2:5) = mean(expData.AS160_dr.sem_d(2:5));
expData.S6_dr.sem_n(2:7)=mean(expData.S6_dr.sem_n(2:7));
expData.ERK_dr.sem_n(2:end)=mean(expData.ERK_dr.sem_n(2:end));
expData.ERK_dr.sem_d(2:end)=mean(expData.ERK_dr.sem_d(2:end));
expData.Foxo_raw_dr.sem_n(1:end-1) = mean(expData.Foxo_raw_dr.sem_n(1:end-1));
expData.Foxo_raw_dr.sem_d(1:end-1) = mean(expData.Foxo_raw_dr.sem_d(1:end-1));
expData.GLUCOSE_dr.sem_n(2:end)=mean(expData.GLUCOSE_dr.sem_n(2:end));
expData.GLUCOSE_dr.sem_d(2:end)=mean(expData.GLUCOSE_dr.sem_d(2:end));

expData.IR_2p.sem_n(:)=mean(expData.IR_2p.sem_n);
expData.IR_2p.sem_d(:)=mean(expData.IR_2p.sem_d);
expData.IR_time.sem(:)=mean(expData.IR_time.sem);
expData.IRS1_2p.sem_n(:)=mean(expData.IRS1_2p.sem_n);
expData.IRS1_2p.sem_d(:)=mean(expData.IRS1_2p.sem_d);
expData.IRS1_time_3.sem(:)=mean(expData.IRS1_time_3.sem);
expData.IRS1_time_30.sem(:)=mean(expData.IRS1_time_30.sem);
expData.IRS1_ds.sem(:)=mean(expData.IRS1_ds.sem);
expData.IRS1307_time.sem_n(:)=mean(expData.IRS1307_time.sem_n);
expData.IRS1307_time.sem_d([2:4 6:end])=mean(expData.IRS1307_time.sem_d([2:4 6:end]));
expData.PKB308_time.sem_n(:)=mean(expData.PKB308_time.sem_n);
expData.PKB308_time.sem_d(1:end-1)=mean(expData.PKB308_time.sem_d(1:end-1));
expData.GLUCOSE_2p.sem_n(:)=mean(expData.GLUCOSE_2p.sem_n);
expData.GLUCOSE_2p.sem_d(:)=mean(expData.GLUCOSE_2p.sem_d);
expData.S6K_time.sem_n(:)=mean(expData.S6K_time.sem_n);
expData.S6K_time.sem_d(:)=mean(expData.S6K_time.sem_d);
expData.S6_time.sem_n(:)=mean(expData.S6_time.sem_n);
expData.S6_time.sem_d(:)=mean(expData.S6_time.sem_d);
expData.ERK_time.sem_n(:)=mean(expData.ERK_time.sem_n);
expData.ERK_time.sem_d(:)=mean(expData.ERK_time.sem_d);
expData.Elk1_time.sem_n(:)=mean(expData.Elk1_time.sem_n);
expData.Elk1_time.sem_d(:)=mean(expData.Elk1_time.sem_d);
expData.Foxo_time.sem_n(:)=mean(expData.Foxo_time.sem_n);
expData.Foxo_time.sem_d(:)=mean(expData.Foxo_time.sem_d);

% Weight down renormalized data
expData.PKB473_time_renorm.sem_n(:)=mean(expData.PKB473_time_renorm.sem_n(expData.PKB473_time_renorm.sem_n~=0));
expData.PKB473_time_renorm.sem_d(:)=mean(expData.PKB473_time_renorm.sem_d(expData.PKB473_time_renorm.sem_d~=0));
expData.AS160_time.sem_n(:)=mean(expData.AS160_time.sem_n);
expData.AS160_time.sem_d(:)=mean(expData.AS160_time.sem_d);

% Data to ignore
expData.IR_dr.sem_n(1)=inf;
expData.IR_dr.sem_d(1)=inf;
expData.IRS1_dr.sem_n([1,7])=inf;
expData.IRS1_dr.sem_d([1:3 7])=inf;
expData.IRS1307_dr.sem_n([1,7])=inf;
expData.IRS1307_dr.sem_d([1,7])=inf;
expData.PKB308_dr.sem_n([1,7])=inf;
expData.PKB308_dr.sem_d(1)=inf;
expData.PKB473_dr.sem_n(1)=inf;
expData.PKB473_dr.sem_d(1)=inf;
expData.AS160_dr.sem_n([1,7])=inf;
expData.AS160_dr.sem_d([1,6,7])=inf;
expData.S6_dr.sem_n(1)=inf;
expData.ERK_dr.sem_n(1)=inf;
expData.ERK_dr.sem_d(1)=inf;
expData.Foxo_raw_dr.sem_n(end)=inf;
expData.Foxo_raw_dr.sem_d(end)=inf;
expData.GLUCOSE_dr.sem_n(1)=inf;
expData.GLUCOSE_dr.sem_d(1)=inf;

expData.IRS1307_time.sem_d([1 5]) = inf;
expData.PKB308_time.sem_d(end)=inf;
expData.PKB473_time_renorm.sem_n(9)=inf;
expData.PKB473_time_renorm.sem_d(9)=inf;

%% setup degrees of freedom
% Dose response dgfs
dgfIR_dr=sum(~isinf([expData.IR_dr.sem_n; expData.IR_dr.sem_d]));
dgfIRS1_dr=sum(~isinf([expData.IRS1_dr.sem_n; expData.IRS1_dr.sem_d]));
dgfIRS1307_dr=sum(~isinf([expData.IRS1307_dr.sem_n; expData.IRS1307_dr.sem_d]));
dgfPKB308_dr=sum(~isinf([expData.PKB308_dr.sem_n; expData.PKB308_dr.sem_d]));
dgfPKB473_dr=sum(~isinf([expData.PKB473_dr.sem_n expData.PKB473_dr.sem_d]));
dgfAS160_dr=sum(~isinf([expData.AS160_dr.sem_n; expData.AS160_dr.sem_d]));
dgfS6_dr=sum(~isinf(expData.S6_dr.sem_n));
dgfERK_dr=sum(~isinf([expData.ERK_dr.sem_n; expData.ERK_dr.sem_d]));
dgfFoxo_raw_dr=sum(~isinf([expData.Foxo_raw_dr.sem_n; expData.Foxo_raw_dr.sem_d]));
dgfGLUCOSE_dr=sum(~isinf([expData.GLUCOSE_dr.sem_n; expData.GLUCOSE_dr.sem_d]));
dgfPKB473_time_renorm=sum(~isinf([expData.PKB473_time_renorm.sem_n; expData.PKB473_time_renorm.sem_d]));

dgfDR=dgfIR_dr+dgfIRS1_dr+dgfIRS1307_dr+dgfPKB308_dr+dgfPKB473_dr+dgfAS160_dr+dgfS6_dr+dgfERK_dr+dgfFoxo_raw_dr+dgfGLUCOSE_dr+dgfPKB473_time_renorm;

% Time series
dgfIR_2p=sum(~isinf([expData.IR_2p.sem_n; expData.IR_2p.sem_d]));
dgfIR_time=sum(~isinf(expData.IR_time.sem));
dgfIR_int=sum(~isinf(expData.IR_int.sem));
dgfIRS1_2p=sum(~isinf([expData.IRS1_2p.sem_n; expData.IRS1_2p.sem_d]));
dgfIRS1_time_3=sum(~isinf(expData.IRS1_time_3.sem));
dgfIRS1_time_30=sum(~isinf(expData.IRS1_time_30.sem));
dgfIRS1_ds=sum(~isinf(expData.IRS1_ds.sem));
dgfIRS1307_time=sum(~isinf([expData.IRS1307_time.sem_n; expData.IRS1307_time.sem_d]));
dgfPKB308_time=sum(~isinf([expData.PKB308_time.sem_n; expData.PKB308_time.sem_d]));
dgfPKB473_time_renorm=sum(~isinf([expData.PKB473_time_renorm.sem_n; expData.PKB473_time_renorm.sem_d])); 
dgfAS160_time=sum(~isinf([expData.AS160_time.sem_n; expData.AS160_time.sem_d])); 
dgfGLUCOSE_2p=sum(~isinf([expData.GLUCOSE_2p.sem_n; expData.GLUCOSE_2p.sem_d]));

if any(strcmp(variables,'measuredS6K'))
    dgfS6K_time=sum(~isinf([expData.S6K_time.sem_n; expData.S6K_time.sem_d]));
else
    dgfS6K_time=0;
end
if any(strcmp(variables,'measuredS6'))
    dgfS6_time=sum(~isinf([expData.S6_time.sem_n; expData.S6_time.sem_d]));
else
    dgfS6_time=0;
end
if any(strcmp(variables,'measuredERK'))
    dgfERK_time=sum(~isinf([expData.ERK_time.sem_n; expData.ERK_time.sem_d]));
else
    dgfERK_time=0;
end
if any(strcmp(variables,'measuredElk1'))
    dgfElk1_time=sum(~isinf([expData.Elk1_time.sem_n expData.Elk1_time.sem_d]));
else
    dgfElk1_time=0;
end
if any(strcmp(variables,'measuredFOXO'))
    dgfFoxo_time=sum(~isinf([expData.Foxo_time.sem_n; expData.Foxo_time.sem_d]));
else
    dgfFoxo_time=0;
end

dgfTS=dgfIR_2p+dgfIR_time+dgfIR_int+dgfIRS1_2p+dgfIRS1_time_3+dgfIRS1_time_30+dgfIRS1_ds+dgfIRS1307_time+dgfPKB308_time+dgfPKB473_time_renorm+dgfAS160_time+dgfGLUCOSE_2p+dgfS6K_time+dgfS6_time+dgfERK_time+dgfElk1_time+dgfFoxo_time;


dgf=dgfDR+dgfTS;
dgf(2)=0; %No explicit data comparisions for validation

%% Parameter bound setup
[pNames, ~] = IQMparameters(model);
nParams=find(strcmp(pNames,'pip'))-1;
expInd=~cellfun(@isempty,(regexp(pNames(1:nParams-2),'^k.+'))); 

lb=repmat(1e-6, 1, nParams-2);
ub=repmat(1e6, 1, nParams-2);

[lb,ub] = CustomBounds(pNames, lb, ub);

lb(expInd) = log(lb((expInd)));
ub(expInd) = log(ub((expInd)));

lb  = [lb 0]; % Adding diabetes
ub  = [ub 1];

%% Setup optimization
opts.ndiverse     = 50; 
opts.maxtime      = 750; 
opts.maxeval      = 1e4;
opts.log_var      = [];

opts.local.solver = 'dhc'; 
opts.local.finish = opts.local.solver;
opts.local.bestx = 0;
opts.local.tol = 1;
opts.local.iterprint = 1;

opts.dim_refset   = 'auto'; %

end

