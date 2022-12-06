function [ structure, parameters, interactions, matchedData] = GenerateLayer( interactions, matchedData )

global dgf
global x

x=[0 15/60 30/60 1 2 5 10 20 60];
structure=[];
parameters=[];
testedInteractions=[];
acceptedData=[];
interactions=shuffle(interactions);
interactions = sortrows(interactions,'Conf','descend');
uniqueTargets=unique(interactions.Target);
uniqueTargets=shuffle(uniqueTargets);
parfor i=1:length(uniqueTargets)
    TargetIndex=ismember(interactions.Target,uniqueTargets{i});
    proteinInteractions=interactions(TargetIndex,:);
    currentMatchedData=matchedData(ismember(matchedData.Genenamesprimary,[proteinInteractions.Source; proteinInteractions.Target]),:);
    [proteinStructure, proteinParameters, tested, cost,simulatedValues]=AddProtein(proteinInteractions, currentMatchedData, dgf,x);

    testedInteractions=[testedInteractions; tested];
    if ~isempty(proteinStructure)
        if cost<chi2inv(0.95,dgf)
            structure=[structure; proteinStructure];
            parameters=[parameters; proteinParameters];
            tmp=matchedData(ismember(matchedData.Genenamesprimary,uniqueTargets{i}),:);
            tmp.simulatedValues=simulatedValues';
            acceptedData=[acceptedData; tmp];
        else
            fprintf('Cost for none-empty struct: %.2f\n',cost)
        end
    end
end

interactions=[testedInteractions; interactions];
[~,ind]=unique(interactions(:,1:2),'rows');
interactions=interactions(ind,:);

matchedData=[acceptedData; matchedData];
[~,ind]=unique(matchedData(:,1:7),'rows');
matchedData=matchedData(ind,:);
end

function [structure, parameters, interactions,cost,simulatedValues]=AddProtein(interactions, matchedData, dgfs, times)
valid=false;
global x
global dgf
dgf=dgfs;
x=times;
structure=[];
parameters=[];
cost=1e99;
simulatedValues=[];
newInteractions=interactions(interactions.Tested==inf,:);
oldInteractions=interactions(interactions.Tested~=inf,:);

if ~any(interactions.Tested<chi2inv(0.95,dgf))&& ~isempty(newInteractions)
    for i =1:height(newInteractions)

        Source=matchedData(ismember(matchedData.Genenamesprimary,newInteractions.Source{i}),:);
        Target=matchedData(ismember(matchedData.Genenamesprimary,newInteractions.Target{i}),:);
        if ~(isempty(Source) || isempty(Target)) && newInteractions.Tested(i)==inf
            disp([newInteractions.Target{i} ' by ' newInteractions.Source{i} ' - Pairwise'])
            [pairStructure, pairParameters, pairCost]=Pairwise(Source, Target);
            newInteractions.Tested(i)=pairCost;
            newInteractions.Parameters{i,1}=pairParameters.Value';
            newInteractions.Type(i)=pairStructure.Type;
            if pairCost<chi2inv(0.95,dgf)
                valid=true;
                structure=pairStructure;
                parameters=pairParameters;
                break
            end
        elseif (isempty(Source) || isempty(Target))
            newInteractions.Tested(i)=1e99;
            break;
        else

        end
    end
    interactions=[oldInteractions; newInteractions];

    if ~valid && height(interactions)>1 && ~isempty(Target)
        disp('Testing multiple input.')
        [interactions,sorting] = sortrows(interactions,'Tested','ascend');
        Sources=matchedData(ismember(matchedData.Genenamesprimary,interactions.Source),:);
        Sources=sortrows(Sources,'Genenamesprimary','ascend');

        Sources=Sources(sorting,:);
        [multiStructure, multiParameters, multiCost]=MultiInput(Sources, Target, interactions);
        if multiCost<chi2inv(0.95,dgf)
            valid=true;
            ind=ismember(interactions.Source,strrep(multiStructure.Source,'_p',''));
            interactions.Parameters(ind)={multiParameters.Value};
            interactions.Tested(ind)=multiCost;
            structure=multiStructure;
            parameters=multiParameters;
        end

    end
    if ~valid && height(newInteractions)>0 && ~ isempty(Target)
        for i =1:height(newInteractions)
            Source=matchedData(ismember(matchedData.Genenamesprimary,newInteractions.Source{i}),:);
            Target=matchedData(ismember(matchedData.Genenamesprimary,newInteractions.Target{i}),:);
            if ~(isempty(Source) || isempty(Target)) && newInteractions.Tested(i)~=1e99
                disp([newInteractions.Target{i} ' by ' newInteractions.Source{i} ' - Plan B'])
                [poolStructure, poolParameters, poolCost]=PlanB(Source, Target);
                if poolCost<chi2inv(0.95,dgf)
                    structure=poolStructure;
                    parameters=poolParameters;

                    newInteractions.Tested(i)=poolCost;
                    newInteractions.Parameters{i,1}=parameters.Value';
                    break
                end
            elseif (isempty(Source) || isempty(Target))
                newInteractions.Tested(i)=1e99;
            else

            end
        end
        interactions=[oldInteractions; newInteractions];
    end

    if ~isempty(parameters)
        [cost,simulatedValues]=CostFunction(parameters.Value');
    elseif exist('pairCost','var')
        cost=pairCost;
    else
        cost=1e99;
    end
    warning('off','MATLAB:DELETE:FileNotFound')
    if ~isempty(Target)
        delete(['Multi_' Target.Genenamesprimary{1} '.*'])
    end
    warning('on','MATLAB:DELETE:FileNotFound')

elseif ~isempty(newInteractions)
    newInteractions.Tested=repmat(1e99,height(newInteractions),1);
    interactions=[oldInteractions; newInteractions];
end
end

function [structure, parameters, cost]=Pairwise(Source, target, useSmooth)
% This function takes data for a Source and a target given in table format and calculates the
% optimal structure for a Pairwise interaction. It outputs this structure,
% as well as optimal parameters.
% The input "useSmooth" is used if data should be smoothed. This is done
% using matlab function csaps with default settings.

global data
global x
global inputParams
global steadyStateParam
global dgf
global modelName
global lbValue
global ubValue

lbValue=log(1e-6);
ubValue=log(1e6);

if nargin <3
    useSmooth=false;
end

format compact
format long

options = optimset('Display','off','algorith','sqp');

try
    yo=Source.simulatedValues;
    if useSmooth
        y=csaps(x,yo,[],x);
    else
        y=yo;
    end
    steadyStateParam=[x repmat(y(1),1,length(x))];
    inputParams=[x y];
    data=target;
    targetName=target.Genenamesprimary{:};
    SourceName=Source.Genenamesprimary{:};

    %% Phosphorylation
    [cost, optParam]=EvaluateModel('phos_model');

    if cost>chi2inv(0.95, dgf)
        [costDephos, optParamDephos]=EvaluateModel('dephos_model');
    else
        costDephos=inf;
    end

    if cost>chi2inv(0.95, dgf) && costDephos>chi2inv(0.95, dgf)
        [costMM, optParamMM]=EvaluateModel('mm_model', [optParam 2]);
    else
        costMM=inf;
    end

    [~,bestFit]=min([costDephos, costMM, cost]);

    %% select which model to return
    if bestFit==1 %% dephos
        modelName='dephos_model';
        cost=costDephos;
        optParam=optParamDephos;

        structure=cell2table({[SourceName '_p'],targetName,'p', 'D'},'VariableNames',{'Source','Target','Site', 'Type'}); % lowercase p is in this instance the "site", which could be a specific site if no mean curves are used.
        parameters=table({...
            ['k_' targetName '_' targetName '_p_f'];...
            ['k_' targetName '_' targetName '_p_b']}, optParam','VariableNames',{'Name','Value'});
    elseif bestFit==2%% mm
        modelName='mm_model';
        cost=costMM;
        optParam=optParamMM;

        structure=cell2table({[SourceName '_p'],targetName,'p', 'PM'},'VariableNames',{'Source','Target','Site', 'Type'}); % lowercase p is in this instance the "site", which could be a specific site if no mean curves are used.
        parameters=table({...
            ['k_' targetName '_' targetName '_p_f'];...
            ['k_' targetName '_' targetName '_p_b'];
            ['km_' targetName '_' targetName '_p']}, optParam','VariableNames',{'Name','Value'});
    else %% phos
        modelName='phos_model';
        structure=cell2table({[SourceName '_p'],targetName,'p', 'P'},'VariableNames',{'Source','Target','Site', 'Type'}); % lowercase p is in this instance the "site", which could be a specific site if no mean curves are used.
        parameters=table({...
            ['k_' targetName '_' targetName '_p_f'];...
            ['k_' targetName '_' targetName '_p_b']}, optParam','VariableNames',{'Name','Value'});
    end

catch err
    disp(getReport(err))
    structure=[];
    parameters=[];
    cost=inf;
end

    function [cost, optParam]=EvaluateModel(model, startGuess)

        modelName=model;
        if nargin<2
            [~,startGuess]=IQMparameters(modelName);
            startGuess=startGuess(find(strcmp(IQMparameters(modelName),'useInterp'))+1:end)';
        end

        lb = lbValue*ones(1,length(startGuess));
        ub = ubValue*ones(1,length(startGuess));
        try
            [~, startCost]=fmincon(@CostFunction,startGuess,[],[], [],[],lb,ub,[],  options);
        catch err
            save('Results/FAILBACKUP-fmincon-EVALUATEMODEL')
            rethrow(err)
        end

        if startCost>=1e90
            particleoptions = optimoptions('particleswarm','SwarmSize',100,'HybridFcn',@fmincon,'Display','iter');
            startGuessLocal=particleswarm(@CostFunction,length(startGuess),lb,ub,particleoptions);
            startCost=CostFunction(startGuessLocal);
        end
        if startCost~=inf
            problem.f = 'CostFunction';
            problem.x_L = lb;
            problem.x_U = ub;
            problem.x_0=startGuess;
            opts.maxtime      = 750; % In cess this option will be overwritten
            opts.maxeval      = 1e4;
            opts.iterprint    = 0;
            opts.local.iterprint = 0;
            opts.local.solver = 'dhc'; %'dhc'; %'fmincon'; %'nl2sol'; %'mix'; %
            opts.local.finish = opts.local.solver;
            opts.local.bestx = 0;
            opts.local.tol = 1;
            opts.dim_refset   = 'auto'; %
            opts.local.use_gradient_for_finish = 0; %DW: provide gradient to fmincon
            opts.local.check_gradient_for_finish = 0; %DW: gradient checker
            optim_algorithm = 'ess'; % 'multistart'; %  'cess'; %

            try
                results = MEIGO(problem,opts,optim_algorithm);
            catch err
                save('Results/FAILBACKUP.mat')
                rethrow(err)
            end

            X = results.xbest;
            try
                optParam=fmincon(@CostFunction,X(1,:),[],[], [],[],lb,ub,[],  options);
            catch err
                disp(getReport(err))
                optParam=X(1,:);
            end
            cost=CostFunction(optParam);
        else
            optParam=startGuess;
        end
    end
end

function [structure, parameters, cost]=MultiInput(Sources, target, interactionsInput)

global data
global steadyStateParam
global inputParams
global modelName
global dgf
global x
data=target;

lbValue=log(1e-6);
ubValue=log(1e6);

format compact
format long

options = optimset('Display','off','Algorithm','sqp');

cost=inf;
n=height(Sources);
optParam=interactionsInput.Parameters{1,:};
optParam=optParam(1:2);

interpStructure=table(repmat({[]},n,1),strcat(Sources.Genenamesprimary,'_p'),repmat({[]},n,1),repmat({'Interp9'},n,1),...
    'VariableNames',{'Source','Target','Site','Type'});

type=repmat({'PO'},n,1);
dephosInd=contains(interactionsInput.Type,'D'); % find interactions that are dephosphorylations
type(dephosInd)=repmat({'DO'},sum(dephosInd),1);

interactions=table(strcat(Sources.Genenamesprimary,'_p'),repmat(target.Genenamesprimary,n,1),repmat({'p'},n,1),type,...
    'VariableNames',{'Source','Target','Site','Type'});

firstOneway=find(strcmp(interactions.Type,'PO') | strcmp(interactions.Type,'DO'),1);
interactions.Type(firstOneway)=strrep(interactions.Type(firstOneway),'O','');
for i = 2:min(height(Sources), 25)
    disp(['Trying with the ' num2str(i) ' with lowest cost, of ' num2str(height(Sources)) ' possible'])
    input=Sources([1 i],:);
    structure=interactions([1 i],:);

    if strcmp(structure.Type{1},'P') && strcmp(structure.Type{2},'PO')
        modelName='Multi_Phos';
    elseif strcmp(structure.Type{1},'D') && strcmp(structure.Type{2},'PO')
        modelName='Multi_Combined';
        structure=[structure(2,:); structure(1,:)];
        structure.Type{1}='PO';
        structure.Type{2}='DO';
        optParam=[interactionsInput.Parameters{2}(1) interactionsInput.Parameters{1}(2)];
        optParam=optParam(1:2);
        input = input([2,1],:);
    elseif (strcmp(structure.Type{1},'P') && strcmp(structure.Type{2},'DO')) || strcmp(structure.Type{1},'D') && strcmp(structure.Type{2},'PO')
        modelName='Multi_Combined';
        structure.Type{1}='PO';
        optParam=[interactionsInput.Parameters{1}(1) interactionsInput.Parameters{2}(2)];
    elseif strcmp(structure.Type{1},'D') && strcmp(structure.Type{2},'DO')
        modelName='Multi_Dephos';
    else
        clear mex
        modelName=@tmpmodel;
        GenerateModel([interpStructure([1 i],:); structure], 'tmpmodel',pwd,100)
        IQMmakeMEXmodel(IQMmodel([modelName '.txt']))
    end
    [paramNames,~]=IQMparameters(modelName);
    paramNames=paramNames(find(strcmp(paramNames,'useInterp'))+1:end)';
    if ~strcmp(modelName, 'Multi_Combined')
        startGuess=[optParam(1:2) lbValue];
    else
        startGuess=optParam(1:2);
    end
    interpParams=[repmat(x,height(input),1) repmat(input.simulatedValues(:,1),1,size(input.simulatedValues,2))];
    steadyStateParam=reshape(interpParams',numel(interpParams),1)';

    interpParams=[repmat(x,height(input),1) input.simulatedValues];
    inputParams=reshape(interpParams',numel(interpParams),1)';

    %% Optimize Model parameters
    lb = repmat(lbValue,1,length(startGuess));
    ub = repmat(ubValue,1,length(startGuess));

    problem.f = 'CostFunction';
    problem.x_L = lb;
    problem.x_U = ub;
    problem.x_0=startGuess;
    opts.maxtime      = 750; % In cess this option will be overwritten
    opts.maxeval      = 1e4;
    opts.iterprint    = 0;
    opts.local.iterprint = 0;
    opts.local.solver = 'dhc'; %'dhc'; %'fmincon'; %'nl2sol'; %'mix'; %
    opts.local.finish = opts.local.solver;
    opts.local.bestx = 0;
    opts.local.tol = 1;
    opts.dim_refset   = 'auto'; %
    opts.local.use_gradient_for_finish = 0; %DW: provide gradient to fmincon
    opts.local.check_gradient_for_finish = 0; %DW: gradient checker
    optim_algorithm = 'ess'; % 'multistart'; %  'cess'; %

    try
        results = MEIGO(problem,opts,optim_algorithm);
    catch err
        save('Results/FAILBACKUP-meigo-combined')
        rethrow(err)
    end

    X = results.xbest;
    try
        optParam=fmincon(@CostFunction,X(1,:),[],[], [],[],lb,ub,[],  options);
    catch err
        disp(getReport(err));
    end
    cost=CostFunction(optParam);
    if cost<chi2inv(0.95,dgf)
        break;
    else
        delete([modelName '.*'])
    end
end


try
    if cost>=chi2inv(0.95, dgf)
        parameters=[];
    else
        parameters=table(paramNames', optParam','VariableNames',{'Name','Value'});
    end
catch err
    disp(getReport(err))
end
end

function [structure, parameters, cost]=PlanB(Source, target, useSmooth)
% [structure, parameters, cost]Pairwise(Source, target, useSmooth)
% tThis function takes data for a Source and a target given in table format and calculates the
% optimal structure for a Pairwise interaction. It outputs this structure,
% as well as optimal parameters.
% The input "useSmooth" is used if data should be smoothed. This is done
% using matlab function csaps with default settings.

global data
global x
global inputParams
global steadyStateParam
global dgf
global modelName

lbValue=log(1e-6);
ubValue=log(1e6);

if nargin <3
    useSmooth=false;
end

format compact
format long

options = optimset('Display','off','Algorithm','sqp');

yo=Source.simulatedValues;
if useSmooth
    y=csaps(x,yo,[],x);
else
    y=yo;
end
steadyStateParam=[x repmat(y(1),1,length(x))];
inputParams=[x y];
data=target;
targetName=target.Genenamesprimary{:};
SourceName=Source.Genenamesprimary{:};

modelName='pool_model';
[~,startGuess]=IQMparameters(modelName);
startGuess=startGuess(find(strcmp(IQMparameters(modelName),'useInterp'))+1:end)';
lb = repmat(lbValue,1,length(startGuess));
ub = repmat(ubValue,1,length(startGuess));

problem.f = 'CostFunction';
problem.x_L = lb;
problem.x_U = ub;
problem.x_0=startGuess;
opts.maxtime      = 750; % In cess this option will be overwritten
opts.maxeval      = 1e4;
opts.iterprint    = 0;
opts.local.iterprint = 0;
opts.local.solver = 'dhc'; %'dhc'; %'fmincon'; %'nl2sol'; %'mix'; %
opts.local.finish = opts.local.solver;
opts.local.bestx = 0;
opts.local.tol = 1;
opts.dim_refset   = 'auto'; %
opts.local.use_gradient_for_finish = 0; %DW: provide gradient to fmincon
opts.local.check_gradient_for_finish = 0; %DW: gradient checker
optim_algorithm = 'ess'; % 'multistart'; %  'cess'; %

try
    results = MEIGO(problem,opts,optim_algorithm);
catch err
    save('Results/FAILBACKUP-meigo-PLANB')
    rethrow(err)
end

X = results.xbest;
optParamPool=fmincon(@CostFunction,X(1,:),[],[], [],[],lb,ub,[],  options);

costPool=CostFunction(optParamPool);
cost=costPool;

%% select which model to return
structure=cell2table({...
    [SourceName '_p'],targetName,'p', 'P';... % Creates the first phosphorylation interaction
    '',[targetName '_p'],[targetName '_U'], 'O';...% Creates the transfer from phosphorylated to unavailable pool (U)
    '',[targetName '_U'], targetName, 'O'},'VariableNames',{'Source','Target','Site', 'Type'}); % creates interaction from U to unphosphorylated.

parameters=table({...
    ['k_' targetName '_' targetName '_p_f'];...
    ['k_' targetName '_' targetName '_p_d'];...
    ['k_' targetName '_p_' targetName '_U'];...
    ['k_' targetName '_U_' targetName]}, optParamPool','VariableNames',{'Name','Value'});
end