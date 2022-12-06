function [] =  PlotAdi(optParam, modelName, expInd, choice, path, useCL_ATP)

if nargin < 1, load('Results/adr_endo, opt(28.3531).mat', 'optParam'); end
if nargin < 2, modelName='adr_endo'; end
if nargin < 4, choice=[]; end
if nargin < 3,  expInd = []; end
if nargin<5, path = './'; end
if nargin<6, useCL_ATP=1; end


colorNorm = [2 64 167]/256;

%% Setup things
fileDir=strrep(path, 'Results','Results-PPL');
bestparam=optParam;

if useCL_ATP
    toEstimateOn={'EPI_ATP', 'CL_ATP','CL_Ca','highCa_ATP', 'Ca_ATP','Ca_noATP','noCa_ATP','noCa_noATP'};
    toValidateOn={};
else
    toEstimateOn={'EPI_ATP','CL_Ca', 'Ca_ATP','Ca_noATP','noCa_ATP','noCa_noATP','highCa_ATP'};
    toValidateOn={'CL_ATP'};
end

[model, expData]=Init_adi(modelName, toEstimateOn, toValidateOn);
toValidateOn={'CL_ATP'};
if isempty(expInd)
    [pNames, ~] = IQMparameters(model);
    nParams=find(strcmp(pNames,'pip'))-1;
    expInd=~cellfun(@isempty,(regexp(pNames(1:nParams-2),'^k.+'))); % -2 for ignoring the diabetes parameters
end

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

experiments=table(stimulus(:,1:3), stimulus(:,4:end), 'VariableNames',{'Pipette','Agonist'},'RowNames',design);

%% Select what to plot
if isempty(choice)
    choice = input('What to Plot?  \n 1. Estimation (Fig 3)  \n 2. Validation (Fig 4)   \n 3. Prediction (Fig 5)   \n 4. Insights (Fig 6)  \n Choice: ');
end
%% Load parameters
if choice ~=5   
    files=dir(sprintf('%s/**/*%s*.mat', fileDir , modelName));
    nFiles=size(files,1);
    params=nan(nFiles,length(bestparam));
    if nFiles > 0
        for p=1:nFiles
            load([files(p).folder '/' files(p).name],'optParam');
            params(p,:)=optParam;
        end
    end
    params(~any(~isnan(params),2),:)=[]; % Removes dummy parameter sets

    disp('Collapsing to unique parameter sets')
    params=unique(params,'rows');
end

%% Get uncertainty and plot
if choice == 1 % Estimation data
    [ boundry, bestSim] = MinMax(model, params, bestparam, expData, experiments(toEstimateOn,:));
    Plot_Exocytosis_Uncertainty (bestSim, boundry, expData)
elseif choice == 2 % Validation
    [ boundry, bestSim] = MinMax(model, params, bestparam, expData, experiments(toValidateOn,:));
    Plot_Exocytosis_Uncertainty (bestSim, boundry, expData)
end

%% Plotting functions
    function [] = Plot_Exocytosis_Uncertainty(sim, boundry, expData)
    design = sim.Properties.RowNames';
    design(strcmp(design,'Time'))=[];

    if length(design)<4
        figure(42);
        hold on
        shadeColor=colorNorm;
        Plotcolor=colorNorm;
        alpha=0.5;
        simTime=sim{'Time','Measures'};
        subplot(2,1,1)
        hold on
        title('Simulation')
        xx=[simTime fliplr(simTime) ];
        yy=[boundry{'CL_ATP','Min'} fliplr(boundry{'CL_ATP','Max'})];
        h=fill(xx,yy,shadeColor);
        set(h, 'EdgeColor','None', 'facealpha',alpha); 
        plot(simTime,sim{'CL_ATP','Measures'},'color', Plotcolor , 'LineWidth' , 1.5 ) %Plot the best parameter set
        xlabel('Time (min)')
        ylabel('\DeltaC/\Deltat (fF/s)')
        axis([-0.2 10 -2.5 20])
        box off

        subplot(2,1,2)
        errorbar(expData{'Time','CL_ATP'},expData{'Mean','CL_ATP'},expData{'SEM','CL_ATP'},'ko', 'MarkerFaceColor','auto')
        title('Experimental data','Interpreter', 'none')
        xlabel('Time (min)')
        ylabel('\DeltaC/\Deltat (fF/s)')
        axis([-0.2 10 -2.5 20])
        box off
        sgtitle('External stimulus, CL no Ca')
        set(figure(42), 'outerposition',[964 50 630 758], 'PaperType','a4')

        figure(53)
        subplot(1,2,1)
        y=sim{'CL_ATP','Measures'};
       
        plot(simTime,y,'--', 'linewidth', 2.5, 'color', colorNorm)
        hold on
        errorbar(expData{'Time','CL_ATP'},expData{'Mean','CL_ATP'},expData{'SEM','CL_ATP'},'^', 'color',colorNorm, 'MarkerFaceColor','white', 'linewidth', 2.5, 'capsize',12)
        xlabel('Time (min)')
        ylabel('\DeltaC/\Deltat (fF/s)')
        axis([-0.2 10 -2.5 20])
        box off
        set(figure(53), 'outerposition',[310 50 1914 566], 'PaperType','a4')
        set(gca, 'FontSize', 15)
    else
        figure(41)
        for pos=1:length(design)
            subplot(4,2,pos)
            Plot_Subplot(boundry, sim, design{pos}, expData)
        end
        figure(41), SetSameAxes(41,[0,2],[],[],0)
        set(gcf, 'OuterPosition',[590,5,804,1336])
    end
    end

    function []=Plot_Subplot(boundry, bestSim,e, data, shadeColor)
    hold on
    if nargin<5 || isempty(shadeColor)
        shadeColor=colorNorm;
        alpha=0.5;
    else
        alpha=0.5;
    end
    bestColor=colorNorm;
    simTime=bestSim{'Time','Measures'};
    if any(strcmp(boundry.Properties.RowNames,e))
        xx=[simTime fliplr(simTime) ];
        yy=[boundry{e,'Min'} fliplr(boundry{e,'Max'})];
        h=fill(xx,yy,shadeColor);
        set(h, 'EdgeColor','None', 'facealpha',alpha);
    end
    if any(strcmp(bestSim.Properties.RowNames,e))
        plot(simTime,bestSim{e,'Measures'},'color', bestColor , 'LineWidth' , 2.5 ) %Plot the best parameter set
    end
    if ~isempty(data) && any(strcmp(data.Properties.VariableNames,e))
        errorbar(data{'Time',e},data{'Mean',e},data{'SEM',e},'o', 'color', colorNorm, 'MarkerFaceColor','auto','linewidth', 2, 'capsize', 12)
    end
    set(gca, 'FontSize', 15)
    title(e,'Interpreter', 'none')
    xlabel('Time (min)')
    ylabel('\DeltaC/\Deltat (fF/s)')
    axis tight
    box off
    end

%% Estimating uncertainty
    function [ boundry, sim] = MinMax(model, params, bestparam, expData, experiments, tend, betaReduction)
    if nargin<6
        tend=12;
    end
    if nargin<7
        betaReduction=[];
    end
    params=[params; bestparam];
    if size(params,1)>1
        fprintf('\nSimulating the uncertainty for the experiments of the adiponectin exocytosis submodel:\n')
    else
        fprintf('\nSimulating the experiments of the adiponectin exocytosis submodel:\n')
    end
    for ii = 1 : size(params,1)
        try
            param = params(ii,:);
            param(expInd)=exp(param(expInd));
            intialconditions = SimulateSteadyState(param, model, 10000/60, betaReduction);
            [sim,~,~, maxcAMP]=SimulateExperiments(param, 0:1/60:tend, intialconditions, model, expData, experiments); % Simulates the experiments.

            sim_valuesStates=sim.States;
            idxTime = ismember(sim{'Time','Measures'}, 15);

            if any(idxTime)
                idxAdi = 4;
                idxCL = strcmp(sim.Properties.RowNames,'CL_ATP');
                sim_valuesStates(1:end-1,:,idxAdi)=sim_valuesStates(1:end-1,:,idxAdi)./sim_valuesStates(idxCL,idxTime,idxAdi);
            end
            if ii == 1
                minValues=sim.Measures;
                maxValues=minValues;

                minValuesStates=sim_valuesStates;
                maxValuesStates=minValuesStates;
            end
            if maxcAMP<inf
                sim_values=sim.Measures;
                minValues(minValues>sim_values)=sim_values(minValues>sim_values);
                maxValues(maxValues<sim_values)=sim_values(maxValues<sim_values);


                minValuesStates(minValuesStates>sim_valuesStates)=sim_valuesStates(minValuesStates>sim_valuesStates);
                maxValuesStates(maxValuesStates<sim_valuesStates)=sim_valuesStates(maxValuesStates<sim_valuesStates);
            end

        catch err
            disp('error in sim');
        end
        if size(params,1)>1
            if ii==1
                fprintf('%i of %i \n|',ii,size(params,1))
            elseif mod(ii,50)==0
                fprintf(' %i of %i \n',ii,size(params,1))
            else
                fprintf('|')
            end
        end
    end
    fprintf('\n\n')
    boundry=table(maxValues, minValues,maxValuesStates, minValuesStates, 'VariableNames',{'Max','Min','MaxStates','MinStates'},'RowNames',sim.Properties.RowNames);
    end
end