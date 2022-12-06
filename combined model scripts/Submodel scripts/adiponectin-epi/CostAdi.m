function [costTot, costVal, sim] = CostAdi(param, expInd, model, expData,estimation, limit, validation)

if nargin<7
    validation=[];
end

 if ~isempty(expInd)
    param(expInd) = exp(param(expInd));
 end


pNames = IQMparameters(model);
param(strcmp(pNames, 'kLdrift'))=0; %disable the drift parameter that is only usefull in lipolysis in vivo simulations
param(contains(pNames, 'diab'))=1; %Disables diabetes
time=unique([0:1/60:max(expData{'Time',:}) expData{'Time',:}]);

SSsim_step=10000/60; % Resolution of the steady state simulation
% SSsim_step = [];

experiments=estimation.Properties.RowNames';

try
    [intialconditions,SteadyState_cost] = SimulateSteadyState(param, model, SSsim_step);
    SS_cost=SteadyState_cost; 
    sim=SimulateExperiments(param, time, intialconditions, model, expData, estimation); % Simulates the experiments. 
    
    peakCost=0;
    cost=0;
    for e = experiments
        tInd=ismember(sim{'Time','Measures'},expData{'Time',e}); % Finds which time points to use when comparing simulation to experimental data
        y=sim{e,'Measures'};
        cost=cost+sum((y(tInd)-expData{'Mean',e}).^2./expData{'SEM',e}.^2); 
        
        if strcmp(e,'EPI_ATP')
            if  sim{e,'PeakTime'}<1.1 || sim{e,'PeakTime'}>1.7
                peakCost=peakCost+limit*(1+abs(1.4-sim{e,'PeakTime'}));
            end
        elseif strcmp(e,'CL_Ca')
            if  sim{e,'PeakTime'}<1.6 || sim{e,'PeakTime'}>2.0
                peakCost=peakCost+limit*(1+abs(1.8-sim{e,'PeakTime'}));
            end
        end
    end
    
    costTot = cost + SS_cost + peakCost; 
    
    if isnan(costTot)
        costTot=1e25;
    end
    
    costVal=0;
    if ~isempty(validation)
        simVal=SimulateExperiments(param, time, intialconditions, model, expData, validation); % Simulates the experiments. 

        for e = validation.Properties.RowNames(1:end-1)
            tInd=ismember(simVal{'Time','Measures'},expData{'Time',e}); % Finds which time points to use when comparing simulation to experimental data
            y=simVal{e,'Measures'};
            costVal=costVal+sum((y(tInd)-expData{'Mean',e}).^2./expData{'SEM',e}.^2); 
        end    
    end

catch err
    %disp(getReport(err))
    costTot=1e25;
    costVal=1e25;
    sim=1e25;
end

end
