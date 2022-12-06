function [intialconditions,SteadyState_cost] = SimulateSteadyState(param, model, SSsim_step, betaReduction)

pNames=IQMparameters(model);
param(ismember(pNames,{'kLdrift', 'kLclear'}))=0; % Disables the drift in experimental data and clearance in vivo of glycerol/FA
param(contains(pNames, 'diab'))=1; %Disables diabetes. If diabetes is missing, two parameters are added

%% Simulate steadystate
simTimeEnd=300000/60;
if nargin<3
    simTime=[0 simTimeEnd-100/60 simTimeEnd];
else
    simTime = unique([0:SSsim_step:simTimeEnd simTimeEnd-100/60]);
end

states=IQMstates(model); 
ic=IQMinitialconditions(model);
% ic(contains(states,'BETA2'))=ic(contains(states,'BETA2'))/100;

if nargin > 4 && ~isempty(betaReduction)
   ic(strcmp(states,'BETA3'))=ic(strcmp(states,'BETA3'))*(1-betaReduction);
end

% input from other models not used in this model
diab = param(end-1:end);
param(end-1:end)=[];
phe0 = 0; % µM
iso0 = 0; % µM
gluc0 = 0;
ins0 = 0; % nM

% input for this model
epi0 = 0; % µM
CL0 = 0;
pip = 0;

SSsim=model(simTime,ic,[param diab pip phe0 iso0 epi0 CL0 gluc0 ins0]); % Steady state simulation.

intialconditions=SSsim.statevalues(end,:);

%% Calculate steady state cost.
SteadyState_cost=0;
adiStates={'BETA3','BETA3a', 'BETA3de', 'Ca', 'ATP', 'cAMP', 'Res', 'Rel', 'PM', 'Endo', 'pipCa', 'pipATP', 'pipcAMP'};    
SSderv= abs((SSsim.statevalues(end,:)-SSsim.statevalues(end-1,:))./(SSsim.time(end) - SSsim.time(end-1)));
SSderv(~ismember(IQMstates(model),adiStates))=[];

largeInd=SSderv > 1e-8;
if any(largeInd)
    SteadyState_cost= 200+sum((SSderv(largeInd)*1e8).^2); 
end


end

