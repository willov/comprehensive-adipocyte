function [responses, steady,  responsesRaw, basal] = simulateInVitro(model, params, diab,  stimulus, inhib)
if nargin<5, inhib=''; end

pNames=IQMparameters(model);
params(ismember(pNames,{'kLdrift', 'kLclear'}))=0; % Disables the drift in experimental data and clearance in vivo of glycerol/FA

noDiabetes = [1 1]; % 1 = no effect of diabetes
[response, ~, steady]=simulate_dose_response(model,params, noDiabetes,stimulus);
responsesRaw.Normal=response;
vars=ismember(response.Properties.VariableNames, {'Glycerol','FA','HSL','PKB473', 'PKB308'});
basal=response(1,vars);
response{:,vars}=response{:,vars}./basal{1,:}*100;

responses.Normal=response;

if ~isequal(diab, [1 1])
    responses.Normal.Reesterification = 100 * (3*responses.Normal.Glycerol - responses.Normal.FA)./(3*responses.Normal.Glycerol);
    [response]=simulate_dose_response(model,params, diab,stimulus);
    responsesRaw.Diabetes=response;
    response{:,vars}=response{:,vars}./basal{1,:}*100;
    responses.Diabetes=response;
    responses.Diabetes.Reesterification = 100 * (3*responses.Diabetes.Glycerol - responses.Diabetes.FA)./(3*responses.Diabetes.Glycerol);
end

if any(ismember(inhib,'rapamycin'))
    idxInhib=ismember(pNames, {'kG5a1', 'kG5a2'});
    inhibParams=params;
    inhibParams(idxInhib)=params(idxInhib)/2;
    [responseInhib, ~, steady]=simulate_dose_response(model,inhibParams, noDiabetes,stimulus);
    responseInhib{:,vars}=responseInhib{:,vars}./basal{1,:}*100;
    responseInhib.Reesterification = 100 * (3*responseInhib.Glycerol - responseInhib.FA)./(3*responseInhib.Glycerol);

    for i=2:7
        inhibParams(idxInhib)=params(idxInhib)/(i*2);
        [response, ~, steady]=simulate_dose_response(model,inhibParams, noDiabetes,stimulus);
        response{:,vars}=response{:,vars}./basal{1,:}*100;
        response.Reesterification = 100 * (3*response.Glycerol - response.FA)./(3*response.Glycerol);
        responseInhib=ConcatenateTableColums(responseInhib, response(:,2:end));
    end
    responses.Rapamycin=responseInhib;
end
if any(ismember(inhib,'torin'))
    idxInhib1=ismember(pNames, {'kG5a1', 'kG5a2'});
    inhibParams=params;
    inhibParams(idxInhib1)=params(idxInhib1)/2;
    idxInhib2=ismember(pNames, {'kG5c'});
    inhibParams(idxInhib2)=params(idxInhib2)/(2*10);
    [responseInhib, ~, steady]=simulate_dose_response(model,inhibParams, noDiabetes,stimulus);
    responseInhib{:,vars}=responseInhib{:,vars}./basal{1,:}*100;
    responseInhib.Reesterification = 100 * (3*responseInhib.Glycerol - responseInhib.FA)./(3*responseInhib.Glycerol);

    for i=2:7
        inhibParams(idxInhib1)=params(idxInhib1)/(i*2);
        inhibParams(idxInhib2)=params(idxInhib2)/(i*20);

        [response, ~, steady]=simulate_dose_response(model,inhibParams, noDiabetes,stimulus);
        response{:,vars}=response{:,vars}./basal{1,:}*100;
        response.Reesterification = 100 * (3*response.Glycerol - response.FA)./(3*response.Glycerol);
        responseInhib=ConcatenateTableColums(responseInhib, response(:,2:end));
    end
    responses.Torin=responseInhib;
end
if any(ismember(inhib,'PD'))
    idxInhib=ismember(pNames, {'kG10basal', 'kG10a1', 'kG10a2'});
    inhibParams=params;
    inhibParams(idxInhib)=params(idxInhib)/2;
    [responseInhib, ~, steady]=simulate_dose_response(model,inhibParams, noDiabetes,stimulus);
    responseInhib{:,vars}=responseInhib{:,vars}./basal{1,:}*100;
    responseInhib.Reesterification = 100 * (3*responseInhib.Glycerol - responseInhib.FA)./(3*responseInhib.Glycerol);

    for i=2:7
        inhibParams(idxInhib)=params(idxInhib)/(i*2);
        [response, ~, steady]=simulate_dose_response(model,inhibParams, noDiabetes,stimulus);
        response{:,vars}=response{:,vars}./basal{1,:}*100;
        response.Reesterification = 100 * (3*response.Glycerol - response.FA)./(3*response.Glycerol);
        responseInhib=ConcatenateTableColums(responseInhib, response(:,2:end));
    end
    responses.PD=responseInhib;

end



end
function [response, basal, steady, sims] = simulate_dose_response(model, params, diab,stimulus)

variables=IQMvariables(model);
observableInd=contains(variables,'y_');
observables=strrep(variables(observableInd),'y_',''); % Find variables containing "y_", and then remove that string from those variables.
inputs=variables(~contains(variables,'y_'));
response=array2table(nan(height(stimulus),length(observables)),'VariableNames',observables);
IR=array2table(nan(height(stimulus),length(inputs)),'VariableNames',inputs);
washInd=ismember(IQMstates(model),{'Gly','FA'});

% input from other models not used in this model
pip0 = 0;
phe0 = 0; % µM
epi0 = 0; % µM
iso0 = 0; % µM
CL0 = 0;
gluc0 = 0;
ins0 = 0; % nM

sim0=model([0 30000-10 30000],[],[params diab pip0 phe0 iso0 epi0 CL0 gluc0 ins0]); % simulate rest before start of experiment

lipoStates = ismember(sim0.states, {'BETA2', 'BETA2a', 'BETA', 'BETAa', 'ALPHA', 'ALPHAa', ...
    'AC', 'ACa', 'PKB_L', 'PKBp', 'PDE3B', 'PDE3Ba', 'cAMP', ...
    'HSL', 'HSLp'});% ignores FA and Glycerol released since they are washed away
steady = abs(sim0.statevalues(end-1,lipoStates)-sim0.statevalues(end,lipoStates)); 
steady = sum(steady(steady>1e-6))*1e9;
ic=sim0.statevalues(end,:);
ic(washInd)=0; % wash away extracellular stuff

for i = fliplr(1:height(stimulus))

    sim0=model(0:15,ic,[params diab pip0 phe0 iso0 epi0 CL0 gluc0 ins0]); %simulate 15 minutes without stimulus

    simIncub=model(0:30,sim0.statevalues(end,:),[params diab pip0 phe0 iso0 epi0 CL0 gluc0 ins0]);% Others
    simIns=model(simIncub.time(end):simIncub.time(end)+15,simIncub.statevalues(end,:),[params diab pip0 phe0 iso0 epi0 CL0 gluc0 stimulus.Ins(i)]);%
    simIso=model(simIns.time(end):simIns.time(end)+10,simIns.statevalues(end,:),[params diab pip0 phe0 stimulus.Iso(i) epi0 CL0 gluc0 stimulus.Ins(i)]);

    response{i,:}=simIso.variablevalues(end,observableInd);

    IR{i,:}=simIso.variablevalues(end,end-length(inputs)+1:end);
    sims(i)=CatSims([simIncub simIns simIso]);
    if i==1
        basal=simIso.statevalues(end,:);
    end
end
response=[stimulus(:,1:2) IR response];
end


% end