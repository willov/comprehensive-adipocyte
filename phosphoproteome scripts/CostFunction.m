function [cost, simulatedValues] = CostFunction(param)
global inputParams
global steadyStateParam
global data
global x
global modelName

if any(isnan(param))
    error('Parameter with NaN value')
end
param=exp(param);
try
    cost=0;
    model=str2func(modelName);
    steadyStateSim=model([0 x(end)*1000-10 x(end)*1000],[],[steadyStateParam 1 param]);
    steadyState=steadyStateSim.statevalues(end,:);
    derivative=abs((steadyStateSim.statevalues(end-1,:)-steadyStateSim.statevalues(end,:))/10);
    limit=1e-4;
    if any(derivative>limit)
        cost=cost+chi2inv(0.95,length(x))*((sum(derivative>limit)>1)+1/limit*sum(derivative));
    end
    sim=model(x,steadyState,[inputParams 1 param]);
    normFactor = sim.statevalues(1,2);
    if normFactor<1e-2
        cost=cost+1/abs(normFactor);
    end
    y = sim.statevalues(:,2)'./normFactor;
    cost=cost+ sum((y-data.meanValues).^2./data.SEMValues.^2);

    if cost==inf
        cost=1e99;
    end
    simulatedValues=sim.statevalues(:,2);
catch err
    cost=1e99;
end

end