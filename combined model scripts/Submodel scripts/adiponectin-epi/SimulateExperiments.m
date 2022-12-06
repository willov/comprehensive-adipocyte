function [ simulatedExperiments, simulations, scale, maxcAMP] = SimulateExperiments(param, times, ic, model, expData, experiments)

%% Initial setups

pNames=IQMparameters(model);
param(ismember(pNames,{'kLdrift', 'kLclear'}))=0; % Disables the drift in experimental data and clearance in vivo of glycerol/FA
param(contains(pNames, 'diab'))=1; %Disables diabetes. If diabetes is missing, two parameters are added

design=matlab.lang.makeValidName({'Ca+ATP';'Ca+noATP';'noCa+ATP';'noCa+noATP'}); %These are the experiments used for scaling. 
%   Ca,           ATP,         cAMP,   Epi, CL
stimulus=[...
    0.0015         3           0.1     0    0;   %'Ca+ATP'
    0.0015         0           0.1     0    0;   %'Ca+noATP'
    0              3           0.1     0    0;   %'noCa+ATP'
    0              0           0.1     0    0];  %'noCa+noATP'
if nargin<6 || isempty(experiments)
    experiments=table(stimulus(:,1:3), stimulus(:,4:end), 'VariableNames',{'Pipette','Agonist'},'RowNames',design);
end

missing=[];
if any(~ismember(design,experiments.Properties.RowNames))
    missing=~ismember(design,experiments.Properties.RowNames);
    experiments(design(missing),:)=table(stimulus(missing,1:3), stimulus(missing,4:end));
end
if nargin==1
    times=expData{'Time',:}; %The last row in expdata is the time
end

scaleExp=experiments.Properties.RowNames;

scaleTime=unique(expData{'Time',scaleExp});
simTimes=unique([0 times scaleTime]);

stateNames=IQMstates(model);
pipIdx = contains(stateNames,'pip');
ic(ismember(stateNames,'Adiponectin'))=0;% Adiponection release is set to zero.
simulations=[];
simulatedExperiments=table(nan(height(experiments),length(simTimes)),'VariableNames',{'Measures'},'RowNames',experiments.Properties.RowNames);
simulatedExperiments.PeakTime=nan(height(simulatedExperiments),1);
maxcAMP=-1;
cAMPInd=strcmp(stateNames,'cAMP');
statesInd=ismember(stateNames,{'Adiponectin','BETA3a','cAMP','Rel'});

%% Simulate experiments
simulatedExperiments.States=nan(height(experiments), length(times), sum(statesInd));

% input from other models not used in this model
diab = param(end-1:end);
param(end-1:end)=[];

phe0 = 0; % µM
iso0 = 0; % µM
gluc0 = 0;
ins0 = 0; % nM

pip = 1;

for i=1:height(experiments)
    ic(pipIdx) = experiments{i,'Pipette'};
    sim=model(simTimes, ic, [param diab pip phe0 iso0 experiments{i,'Agonist'} gluc0 ins0]);
    maxcAMP=max([maxcAMP; sim.statevalues(:,cAMPInd)]);
    simulations=[simulations sim]; % Collect the full simulation structure, but only for the time points requested.
    simulatedExperiments{i,'Measures'}=simulations(end).variablevalues(:,1)';
    simulatedExperiments{i,'States'}=reshape(sim.statevalues(:,statesInd), 1, length(times), sum(statesInd));
    simulatedExperiments.PeakTime(i)=sim.time(find(sim.variablevalues(:,1)==max(sim.variablevalues(:,1)),1));
end


%% Do scaling
meanValues=expData{'Mean',scaleExp}';
SEMValues=expData{'SEM',scaleExp}';
sims=[];
for exp = scaleExp'
    sim=simulatedExperiments{exp,'Measures'};
    sims=[sims; sim(:,ismember(simTimes,expData{'Time',exp}))'];
end

sims(isnan(meanValues))=[];
SEMValues(isnan(meanValues))=[];
meanValues(isnan(meanValues))=[];
scale=lscov(sims,meanValues, SEMValues.^-2);
if scale<0
    scale=1;
end

simulatedExperiments{:,'Measures'}=simulatedExperiments{:,'Measures'}.*scale;

%% Setup final output. 

simulatedExperiments(end+1,:)={simulations(end).time nan nan(1, length(times), sum(statesInd))};
simulatedExperiments.Properties.RowNames{end}='Time';
tInd=ismember(simTimes,times);
simulatedExperiments.Measures(:,~tInd)=[]; %Removes extra timepoints only used for scaling.
simulatedExperiments(design(missing),:)=[]; % we might have added simulations not asked for because we needed them for scaling
simulations(end-sum(missing)+1:end)=[]; % See above. The number of missing is the amount we need to remove. 

end



