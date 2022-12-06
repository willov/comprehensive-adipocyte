function [diabetesEffect, totalInhibitions, inhibitionValuesDirection, diabeticSimData]=PlotFinalModel(parameters, time, expData, modelName, coreModel, list, simulateDiabetes, plotNormal, plotDiabetes, plotInhib, simulateInhib)

if nargin<7
    simulateDiabetes=false;
end
if nargin<8
    plotNormal=false;
end
if nargin<9
    plotDiabetes=false;
end
if nargin<10
    plotInhib=false;
end
if nargin<11
    simulateInhib=false;
end

if plotDiabetes
    simulateDiabetes=true;
end

if ~exist('../Phosphoproteome-plots/Normal','dir')  && plotNormal
    mkdir('../Phosphoproteome-plots/Normal')
end
if ~exist('../Phosphoproteome-plots/Diabetes','dir') && plotDiabetes
    mkdir('../Phosphoproteome-plots/Diabetes')
end
if ~exist('../Phosphoproteome-plots/Inhibitions/MK','dir') && plotInhib
    mkdir('../Phosphoproteome-plots/Inhibitions/MK')
end
if ~exist('../Phosphoproteome-plots/Inhibitions/LY','dir') && plotInhib
    mkdir('../Phosphoproteome-plots/Inhibitions/LY')
end

close all

timeHR = unique([0:0.1:time(end) time]);
model=str2func(modelName);

% Update names to handle problematic characters
expData.Genenamesprimary=strrep(expData.Genenamesprimary,'-','');
list{:,1:2}=strrep(list{:,1:2},'-','');

range=1:5;
[~, normalSimData, Bmarkinterp, diabeticSimData, inhibSimData]=GetInputs(expData, time, coreModel, plotInhib || simulateInhib, range);
diabeticSimData.Effect = diabeticSimData.simulatedValues(:,time==20)./normalSimData.simulatedValues(:,time==20);
parameters.Value=exp(parameters.Value);
normalInput=innerjoin(table(regexprep(Bmarkinterp.Target,'_[up]+','')),normalSimData,'LeftKey','Var1','RightKey','HumphreyName');

interpParams=[repmat(time,height(normalInput),1) repmat(normalInput.simulatedValues(:,1),1,size(normalInput.simulatedValues,2))];
interpParams=reshape(interpParams',numel(interpParams),1);
normalSteadyStateParam=[interpParams; 1; parameters.Value];

interpParams=[repmat(time,height(normalInput),1) normalInput.simulatedValues];
interpParams=reshape(interpParams',numel(interpParams),1);
normalParam=[interpParams; 1; parameters.Value];

disp('Simulating the model under normal conditions')
normalSteadyStateSim=model([0 time(end)*1000],[],normalSteadyStateParam);
normalSS=normalSteadyStateSim.statevalues(end,:);
normalSim=model(timeHR,normalSS,normalParam);
normIC = normalSim.statevalues(1,:);
normalSim.statevalues=normalSim.statevalues./normIC;

if simulateDiabetes %also performs the diabetic simulation if no plotting are true, to be able to generate the list of interactions with diabetes effect without plotting figures
    diabeticInput=innerjoin(table(regexprep(Bmarkinterp.Target,'_[up]+','')),diabeticSimData,'LeftKey','Var1','RightKey','HumphreyName');

    interpParams=[repmat(time,height(diabeticInput),1) repmat(diabeticInput.simulatedValues(:,1),1,size(diabeticInput.simulatedValues,2))];
    interpParams=reshape(interpParams',numel(interpParams),1);
    diabeticSteadyStateParam=[interpParams; 1; parameters.Value];

    interpParams=[repmat(time,height(diabeticInput),1) diabeticInput.simulatedValues];
    interpParams=reshape(interpParams',numel(interpParams),1);
    diabeticParam=[interpParams; 1; parameters.Value];

    disp('Simulating the model under diabetic conditions')
    diabeticSteadyStateSim=model([0 time(end)*1000],[],diabeticSteadyStateParam);
    diabeticSS=diabeticSteadyStateSim.statevalues(end,:);
    diabeticSim=model(timeHR,diabeticSS,diabeticParam);
    diabeticSim.statevalues=diabeticSim.statevalues./normIC;
    [phosphorylations,~]=regexp(IQMstates(modelName),'(\w+)_[p]','tokens','match');
    ind=~cellfun(@isempty,phosphorylations);
    tind=ismember(diabeticSim.time,20);
    names=cell2table([phosphorylations{ind}]','VariableNames',{'Protein'});
    change=log2(diabeticSim.statevalues(tind,ind)./normalSim.statevalues(tind,ind));

    diabetesEffect=[names table(change','VariableNames',{'Effect'})];
else
    diabetesEffect=[];
    diabTime=[];
    diabVal=[];
end

if plotInhib || simulateInhib% simulate inhibitions
    inhibSims=table();
    inhibitions = ['Normal';inhibSimData.Properties.VariableNames(1:end-1)'];
    sims=cell(height(inhibitions)-1,1);
    disp('Simulating inhibitions')

    for i=1:width(inhibSimData)-1
        inhibInput=innerjoin(table(regexprep(Bmarkinterp.Target,'_[up]+','')),inhibSimData(:,[i end]),'LeftKey','Var1','RightKey','HumphreyName');

        values=inhibInput{:,2};
        interpParams=[repmat(time,height(inhibInput),1) repmat(values(:,1),1,size(values,2))];
        interpParams=reshape(interpParams',numel(interpParams),1);
        inhibSteadyStateParam=[interpParams; 1; parameters.Value];

        interpParams=[repmat(time,height(inhibInput),1) values];
        interpParams=reshape(interpParams',numel(interpParams),1);
        inhibParam=[interpParams; 1; parameters.Value];

        inhibSteadyStateSim=model([0 time(end)*1000],[],inhibSteadyStateParam);
        inhibSS=inhibSteadyStateSim.statevalues(end,:);
        inhibSim = model(timeHR,inhibSS,inhibParam);
        inhibSim.statevalues=inhibSim.statevalues./normIC;
        sims{i}= inhibSim;
        disp(['Done with: ' num2str(i) ' of ' num2str(width(inhibSimData)-1)])
    end
    inhibSims.Sim=vertcat({normalSim}, sims);
    inhibSims.Properties.RowNames=inhibitions;
end

%% Plot final model
[phosphorylations,~]=regexp(normalSim.states,'(\w+)_p','tokens','match');
count=0;
for i=1:length(phosphorylations)
    if ~isempty(regexprep(char(phosphorylations{i}{:}), '_.+','')) %Ignores the sites without an identified gene
        protein=char(phosphorylations{i}{1});
        if regexp(protein,'^_[0-9A-Z]{2,}')
            protein=protein(2:end);
        end
        interactions = list(ismember(list.Target,protein),:);
        interactions(interactions.Tested~=min(interactions.Tested),:)=[]; %Remove all tested, not best, interactions
        layer = unique(interactions.ExpansionIteration);
        maxLayer=max(list.ExpansionIteration(list.Tested<chi2inv(0.95,length(time)-1))); % only works if loading the final layer
        if maxLayer==5
            maxLayer = 77; % To manage the colors if only plotting the first
        end
        data=expData(ismember(expData.Genenamesprimary,protein),:);
        data.simulatedValues=data.simulatedValues/data.simulatedValues(1);

        cost=sum((data.meanValues-normalSim.statevalues(ismember(normalSim.time,time),i)').^2./data.SEMValues.^2);
        count=count+1;

        normTime=normalSim.time;
        normVal=normalSim.statevalues(:,i);
        if plotDiabetes
            diabTime=diabeticSim.time;
            diabVal=diabeticSim.statevalues(:,i);
        end

        inhibSim=[];
        if plotInhib
            for j=2:height(inhibSims)
                sim=inhibSims.Sim{j};
                inhibSim=[inhibSim; sim.statevalues(:,i)'];
            end
            inhibSim=[inhibSim; sim.time];
        end
        source = regexprep(sprintf('%s+', interactions.Source{:}),'+$','');
        if sum(~cellfun(@isempty, (regexp(normalSim.states,[protein '_U']))))>0
            titleStr=sprintf('%s -> %s - cost: %.2f, layer: %i (pool)',source, protein, cost, layer);
        else
            titleStr=sprintf('%s -> %s - cost: %.2f, layer: %i',source, protein, cost, layer);
        end

        %% Plot
        if plotNormal
            PlotSims(protein, layer, maxLayer, time, data, normTime,normVal,[],[], [],'', titleStr); % plot normal
        end
        if plotDiabetes
            PlotSims(protein, layer, maxLayer, time, data, normTime,normVal,diabTime,diabVal, [],'', titleStr); % plot diabetes
        end
        if plotInhib
            type.type='MK';
            type.values=round((1-1./(2.^range))*100);
            PlotSims(protein, layer, maxLayer, time, data, normTime,normVal,[],[], inhibSim([1:5 end],:),type, titleStr); % plot MK inhibitions
            type.type='LY';
            PlotSims(protein, layer, maxLayer, time, data, normTime,normVal,[],[], inhibSim(6:end,:),type,titleStr); % plot plot LY inhibitions
        end
    end
    if mod(i,100)==0
        close(gcf)
    end
end

if simulateInhib
    [totalInhibitions, inhibitionValuesDirection] = SimulateInhib(expData, modelName, inhibSims);
else
    totalInhibitions=[];
    inhibitionValuesDirection=[];
end
end

function [] = PlotSims(protein, layer, maxLayer, time, data, normTime,normVal,diabTime,diabVal, inhibSim, inhibType, titleStr)
clf

colorNorm = [2 64 167]/256;
colorDiab = [0.6350 0.0780 0.1840];
colorInhib = [32 237 178]/256;

cStart = [8 69 148]/256;
cEnd = [230 238 247]/256;
layerGradient = [linspace(cStart(1),cEnd(1),maxLayer)', linspace(cStart(2),cEnd(2),maxLayer)', linspace(cStart(3),cEnd(3),maxLayer)'];

cStart = layerGradient(layer,:);
cEnd = colorInhib;
inhibGradient = [linspace(cStart(1),cEnd(1),6)', linspace(cStart(2),cEnd(2),6)', linspace(cStart(3),cEnd(3),6)'];

plot(normTime,normVal,'linewidth',3,'Color',layerGradient(layer,:))
hold on

if ~isempty(diabVal)
    plot(diabTime,diabVal,'linewidth',3,'Color',colorDiab,'linewidth',3)
    hold on
end

if ~isempty(inhibSim)
    for j=1:5
        plot(inhibSim(end,:),inhibSim(j,:),'color',inhibGradient(j+1,:),'linewidth',3)
        hold on
    end

    if strcmp(inhibType.type,'MK')
        value=data.meanValues(1,ismember(time,20))*(data.MKMean);
        upper=data.meanValues(1,ismember(time,20))*(data.MKMean+data.MKSEM);
        lower=data.meanValues(1,ismember(time,20))*(data.MKMean-data.MKSEM);
    elseif strcmp(inhibType.type,'LY')
        value=data.meanValues(1,ismember(time,20))*(data.LYMean);
        upper=data.meanValues(1,ismember(time,20))*(data.LYMean+data.LYSEM);
        lower=data.meanValues(1,ismember(time,20))*(data.LYMean-data.LYSEM);
    else
        disp('Something went wrong.. Incorrect inhibition type!')
    end
    errorbar(20,value,upper-lower ,'o','linewidth',3,'Color',colorInhib, 'CapSize',12, 'MarkerFaceColor','auto')
end

errorbar(time,data.meanValues,data.SEMValues,'o','linewidth',3,'Color',colorNorm,'CapSize',12, 'MarkerFaceColor','auto') % Plot data

axis tight
box off
xlabel('Time (min)')
ylabel('Fold over basal')
set(gca,'fontsize',26)
set(gcf, 'Units','pixels', 'outerposition',[0 0 2560 1440], 'PaperType','a4')

%% setup legend
entries={'Normal'}';

if ~isempty(diabVal)
    entries=[entries; {'Diabetic'}];
end
if ~isempty(inhibSim)
    inhibEntries=cellstr(strcat(inhibType.type,num2str(inhibType.values')));
    entries=[entries; inhibEntries; {[inhibType.type ' data']}];
end
entries=[entries; {'Normal data'}];

legend(entries,'location','best')
title(titleStr, 'Interpreter','none')
%% name figure
if ~isempty(diabVal)
    name=sprintf('../Phosphoproteome-plots/Diabetes/%i_%s_diabetes',layer, protein);
elseif ~isempty(inhibType) && strcmp(inhibType.type,'MK')
    name=sprintf('../Phosphoproteome-plots/Inhibitions/MK/%i_%s_MK', layer, protein);
elseif ~isempty(inhibType) && strcmp(inhibType.type,'LY')
    name=sprintf('../Phosphoproteome-plots/Inhibitions/LY/%i_%s_LY', layer, protein);
else
    name=sprintf('../Phosphoproteome-plots/Normal/%i_%s',layer, protein);
end
print('-djpeg', '-r0', [name '.png'])
exportgraphics(gcf, [name '.pdf'], 'ContentType','vector')
end

function [totalInhibitions, inhibitionValuesDirection] = SimulateInhib(expData, modelName, sims)
expData(:,strcmp(expData.Properties.VariableNames, 'Cluster'))=[];

%% Simulate inhibitions
[phosphorylations,~]=regexp(IQMstates(modelName),'(\w+)_[p]','tokens','match');
ind=~cellfun(@isempty,phosphorylations);
inhibitionValues=cell2table([phosphorylations{ind}]','VariableNames',{'Protein'});

nInhibs=height(sims)-1;

for i=1:nInhibs+1
    sim=sims.Sim{i};
    inhibitionValues=[inhibitionValues table(sim.statevalues(sim.time==20,ind)','VariableNames',sims.Properties.RowNames(i))];
end

inhibitionValues{:,3:end}=inhibitionValues{:,3:end}./repmat(inhibitionValues{:,2},1,nInhibs); % Normalized to fold over insulin without inhibition
inhibitionValues.pSim = sign(log2(inhibitionValues.Normal));
inhibitionValues.Normal=[];
inhibCol={'Genenamesprimary', 'meanValues', 'LYDirection', 'LYMean','LYSEM', 'MKDirection', 'MKMean', 'MKSEM'};
inhibitionValues=innerjoin(inhibitionValues,expData(:,inhibCol),'LeftKey','Protein','RightKey','Genenamesprimary');

inhibitions=[repmat(inhibitionValues.LYMean,1,nInhibs/2) repmat(inhibitionValues.MKMean,1,nInhibs/2)];
SEMs=[repmat(inhibitionValues.LYSEM,1,nInhibs/2) repmat(inhibitionValues.MKSEM,1,nInhibs/2)];

%% Get inhibitions with same directions
inhibitionValuesDirection=inhibitionValues(:,1:nInhibs+1);
simulatedDirections=inhibitionValuesDirection{:,2:nInhibs+1};
simulatedDirections(simulatedDirections<1)=-1;
simulatedDirections(simulatedDirections>1)=1;

directions=zeros(size(inhibitions));
directions(inhibitions+SEMs<1)=-1;
directions(inhibitions-SEMs>1)=1;

inhibitionValuesDirection{:,2:nInhibs+1}=abs(simulatedDirections+directions)==2;% & directions==-1;

%% Convert to log2 scale
inhibitionValues{:,2:nInhibs+1} = log2(inhibitionValues{:,2:nInhibs+1});

%% Collect into final table
totalInhibitions = CompileList();

%% set unclear directions to nan
agreementDirections = inhibitionValuesDirection{:,2:nInhibs+1};
agreementDirections(directions==0)=nan;
inhibitionValuesDirection{:,2:nInhibs+1}=agreementDirections;% set unclear directions to nan for use in InhibitionsForAllLayers.

%% Subfunctions
    function [totalInhibitions]=CompileList()
        totalInhibitions=inhibitionValues(1,2:nInhibs/2+1);
        totalInhibitions.Properties.VariableNames=strrep(totalInhibitions.Properties.VariableNames,'MK','I');

        % Correct direction
        totalInhibitions{1,:}=sum(inhibitionValuesDirection{:,2:nInhibs/2+1})./sum(directions(:,1:nInhibs/2)~=0);
        totalInhibitions{2,:}=sum(inhibitionValuesDirection{:,nInhibs/2+2:nInhibs+1})./sum(directions(:,nInhibs/2+1:end)~=0);
        totalInhibitions{3,:}=sum(inhibitionValuesDirection{:,nInhibs/2+2:nInhibs+1} & inhibitionValuesDirection{:,2:nInhibs/2+1})./sum(directions(:,1:nInhibs/2)~=0 & directions(:,nInhibs/2+1:end)~=0);

        % Unclear
        totalInhibitions{4,:}=sum(directions(:,1:nInhibs/2)==0)./height(directions);
        totalInhibitions{5,:}=sum(directions(:,nInhibs/2+1:nInhibs)==0)./height(directions);
        totalInhibitions{6,:}=sum(directions(:,1:nInhibs/2)==0 & directions(:,nInhibs/2+1:nInhibs)==0)./height(directions);

        totalInhibitions.Properties.RowNames={'direction_MK';'direction_LY';'direction_MK+LY';'Unclear_MK';'Unclear_LY'; 'Unclear_MK+LY'};
    end

end
