
function [list, parameters, structure, expData]=ExpandModel(plotNormal, plotDiabetes, plotInhib)

if nargin<1
    plotNormal=false;
end
if nargin<2
    plotDiabetes=false;
end
if nargin<3
    plotInhib=false;
end

coreModel='combined';
run('../combined model scripts/Setup_combined')
Setup_phosphoproteome()

try
    global x
    currentPath=pwd;

    %% Load stuff
    matchedData=[];

    matchedData=GetPhosphoproteomeData(matchedData);
    matchedData(ismember(matchedData.Genenamesprimary, ''),:)=[];
    repeats(:,:,1)=matchedData.TC_Exp1;
    repeats(:,:,2)=matchedData.TC_Exp2;
    repeats(:,:,3)=matchedData.TC_Exp3;
    matchedData(any(sum(~isnan(repeats),3)<2,2),:)=[]; % Remove data having less than two repeats
    matchedData(strcmp(matchedData.Cluster,''),:)=[]; % Remove unnecessary column

    %% Do inital stuff
    global dgf
    global nBetterMulti
    nBetterMulti=0;
    matchedData.SEMValues(:,1)=inf;

    list=[];

    x=[0 15/60 30/60 1 2 5 10 20 60];
    dgf=length(x)-1;
    [expData, ~, Bmarkinterp]=GetInputs(matchedData, x, coreModel);

    ind=[];
    [~,ind]=unique(expData.Genenamesprimary);
    ind=~ismember(1:height(expData),ind)' & ismember(expData.Genenamesprimary,strrep(Bmarkinterp.Target,'_p',''));
    expData(ind,:)=[];
    [~, idx] = unique(expData.Genenamesprimary); %Remove duplicate sites
    expData=expData(idx,:);

    structure=Bmarkinterp;
    parameters=[];
    layer=0;
    initialProteins=strrep(structure.Target,'_p','');

    %% Setup OmniPath interactions
    load('Interactions/interactionsOmni.mat','interactionsOmni')
    coreIdx = ismember(expData.Genenamesprimary, initialProteins);

    sites=table();
    sites.Target=regexprep(expData.Genenamesprimary(~coreIdx),'_.+','');
    sites.new = expData.Genenamesprimary(~coreIdx);

    interactionsOmni=innerjoin(interactionsOmni, sites);
    interactionsOmni.Target=interactionsOmni.new;
    interactionsOmni.new=[];

    sites.Properties.VariableNames{1}='Source';
    core=table(structure.Target, structure.Target, 'VariableNames',{'Source','new'});
    core{:,:}=strrep(core{:,:},'_p','');
    core.Source=regexprep(core.Source,'_.+','');
    sites=[core;sites];
    interactionsOmni=innerjoin(interactionsOmni, sites);
    interactionsOmni.Source=interactionsOmni.new;
    interactionsOmni.new=[];

    % Divide responding and nonresponding data and interactions
    load('responders_sequence.mat','responders_sequence');
    expDataNonResponders=expData(~ismember(expData.SequenceWindow,responders_sequence.SequenceWindow) & ~ismember(expData.Source,'Start'),:);
    expData = expData(ismember(expData.SequenceWindow,responders_sequence.SequenceWindow) | ismember(expData.Source,'Start'),:);

    interactionsResponders=interactionsOmni(all(ismember(interactionsOmni{:,1:2}, expData.Genenamesprimary),2),:); % Interactions where both source and  target are in the responder data
    interactionsResponders.DataSelection(:) = {'Responders'};

    interactionsNonResponders=interactionsOmni(~all(ismember(interactionsOmni{:,1:2}, expData.Genenamesprimary),2),:); % Interactions where both source and  target are in the responder data
    interactionsNonResponders.DataSelection(:) = {'Nonresponders'};

    %% Expand using the data of sites responding to insulin, and with decreasing reliability of interactions
    interactionsNew=[];
    confLevels = sort(unique(interactionsResponders.Conf),'descend')';
    for rel = confLevels
        interactionsNew = [interactionsNew; interactionsResponders(interactionsResponders.Conf==rel,:)];
        list=EvaluteListInteractions(list, interactionsNew);
    end
    save('Results/PriorDriven-responders.mat') % Can be used to restart the algorithm from this point later on if necessary

    % load('Results/PriorDriven-responders.mat','list','expData','structure','parameters') % Load to restart the algorithm from this point later on if necessary

    %% Expand using all available data
    list(list.Tested==inf,:)=[];
    expData=[expData; expDataNonResponders];

    interactionsNew=[];
    confLevels = sort(unique(interactionsNonResponders.Conf),'descend')';
    for rel = confLevels
        interactionsNew = [interactionsNew; interactionsNonResponders(interactionsNonResponders.Conf==rel,:)];
        list=EvaluteListInteractions(list, interactionsNew);
    end
    save('Results/PriorDriven-all.mat') % Can be used to restart the algorithm from this point later on if necessary

    % load('Results/PriorDriven-all.mat','list','expData','structure','parameters') % Load to restart the algorithm from this point later on if necessary

    %% Save prior driven results and run datadriven extensions
    rel = -1;
    list=ExpandDatadriven(list, expData);
    save('Results/DataDriven.mat') % Can be used to restart the algorithm from this point later on if necessary

    % load('Results/DataDriven.mat','list','expData','structure','parameters') % Load to restart the algorithm from this point later on if necessary

    %% Cleanup after model expansion
    layers = unique(list.ExpansionIteration);
    layers(layers==inf)=[];

    for l = 1:length(layers)
        list.ExpansionIteration(list.ExpansionIteration==layers(l))=l;
    end

    maxLayer = max(list.ExpansionIteration(~isinf(list.ExpansionIteration)));
    list.ExpansionIteration(list.Tested>chi2inv(0.95,dgf) & list.Tested<inf)=maxLayer+10; %Add all tested, not valid interactions to a distant layer.

    %% Generate, simulate and plot final model
    % load('Results/PriorDriven-responders.mat','list','expData','structure','parameters')

    modelName='FinalExpandedModel';

    if ~isempty(structure)
        GenerateModel(structure,modelName,pwd,100);
        IQMmakeMEXmodel(IQMmodel([modelName '.txt']))
    end

    simulateDiabetes=1; % Diabetes must be simulated for saving the results later
    [diabetesEffect, ~, ~,  diabeticSimData]=PlotFinalModel(parameters,x, expData, modelName, coreModel, list, simulateDiabetes, plotNormal,plotDiabetes,plotInhib, plotAllInhib);

    %% Save list of interactions
    testedInteractions=list(list.Tested<chi2inv(0.95,dgf),[1:2 5 3 3 4 4 7 7 8 8]);
    if ~isempty(diabetesEffect)
        testedInteractions=outerjoin(testedInteractions,diabetesEffect,'Leftkey','Target','rightkey','Protein','Type','left');
        testedInteractions.Protein=[];
    end

    load('Translation.mat', 'Translation')
    Translation.BmarkName = regexprep(Translation.BmarkName, {'measured', 'y_'},{'',''});
    testedInteractions.Source=regexprep(testedInteractions.Source, Translation.HumphreyName, Translation.BmarkName);
    testedInteractions.Source=regexprep(testedInteractions.Source, {'Pde3b'}, {'PDE3B'});

    load('Interactions/coreInteractions.mat', 'coreInteractions');
    diabeticSimData.simulatedValues=[];
    diabeticSimData(cellfun(@isempty, diabeticSimData.HumphreyName),:)=[];
    diabeticSimData.HumphreyName=regexprep(diabeticSimData.HumphreyName, Translation.HumphreyName, Translation.BmarkName);
    diabeticSimData.HumphreyName=regexprep(diabeticSimData.HumphreyName, {'Pde3b'}, {'PDE3B'});
    diabeticSimData=unique(diabeticSimData,'rows');
    coreInteractions=outerjoin(coreInteractions, diabeticSimData,'LeftKeys','Target', 'RightKeys','HumphreyName');
    coreInteractions.HumphreyName=[];
    writetable([coreInteractions; testedInteractions],'Results/testedInteractions.xlsx','WriteMode','overwrite')

    cd(currentPath);
    clear mex

catch err
    disp(getReport(err))
    save('Results/Crash Backup')
end

%% Subfunctions
    function [list]=ExpandDatadriven(list, expData)
        if ~ismember('Priority',list.Properties.VariableNames)
            list.Priority=zeros(height(list),1);
        end
        while true

            addedInd=list.Tested<chi2inv(0.95,dgf);
            notAdded=expData(~ismember(expData.Genenamesprimary,unique([list.Source(addedInd); list.Target(addedInd)])),:);
            added=expData(~ismember(expData.Genenamesprimary,notAdded.Genenamesprimary),:);
            tmpList=[];
            for i=1:height(notAdded)
                y=repmat(notAdded.meanValues(i,:),height(added),1);
                SEM=repmat(notAdded.SEMValues(i,:),height(added),1);
                costs=sum(((y-added.simulatedValues).^2)./SEM.^2,2);
                [costs,costInd]=sort(costs,'ascend');
                Sources=added.Genenamesprimary(costInd);
                validInd=costs<chi2inv(0.95,dgf)+5;
                Sources(~validInd | costs==inf | isnan(costs))=[];
                tmpList=[tmpList; table(Sources, repmat(notAdded.Genenamesprimary(i),length(Sources),1), costs(validInd),'VariableNames',{'Source','Target','Priority'})];
                disp(['Generating interactions for ' num2str(i) ' of ' num2str(height(notAdded))]);
            end

            tmpList.Conf(:)=-1;
            tmpList.DataSelection(:) = unique(list.DataSelection(list.ExpansionIteration==max(list.ExpansionIteration(~isinf(list.ExpansionIteration)))));
            tmpList.Tested=inf(height(tmpList),1);
            tmpList.Parameters=repmat({''},height(tmpList),1);
            tmpList.ExpansionIteration=inf(height(tmpList),1);
            tmpList.Type(:)={''};
            tmpList=tmpList(:,[1 2 4:end 3]);
            list=[list; tmpList];

            [~,ind]=unique(list(:,1:2),'rows');
            list=list(ind,:);

            interactions=GetPriorList(structure, list, expData);
            if isempty(interactions) || sum(interactions.Tested==inf)==0
                break;
            end
            [layerStructure, layerParameters, interactions, expData]=GenerateLayer(interactions, expData);
            layer=layer+1;
            structure=[structure; layerStructure];
            parameters=[parameters; layerParameters];
            interactions.ExpansionIteration(interactions.ExpansionIteration==inf & interactions.Tested<chi2inv(0.95,dgf))=layer;
            interactions(ismember(interactions(interactions.Tested>chi2inv(0.95,dgf),1:2), list(list.Tested<chi2inv(0.95,dgf),1:2)),:)=[]; %if a datadriven interaction is also prior, remove it if it does not pass a chi2-test.
            list(ismember(list(:,1:2),interactions(interactions.Tested~=inf,1:2),'rows'),:)=interactions(interactions.Tested~=inf,:);

            added=interactions(interactions.Tested<chi2inv(0.95,dgf),:);
            if ~isempty(layerStructure)
                save(['Results/layer' num2str(layer) '.mat'])
            else
                save(['Results/empty' num2str(layer) '.mat'])
            end
        end
    end

    function [list]=EvaluteListInteractions(list, newInteractions)
        tmpList=newInteractions;
        tmpList.Tested=inf(height(tmpList),1);
        tmpList.Parameters=repmat({''},height(tmpList),1);
        tmpList.ExpansionIteration=inf(height(tmpList),1);
        tmpList.Type=repmat({''},height(tmpList),1);
        if ~isempty(list)
            tmpList(ismember(tmpList.Target,list.Target(list.Tested<chi2inv(0.95,dgf))),:)=[];
        end
        list=[list; tmpList];
        [~,ind]=unique(list(:,1:2),'rows');
        list=list(ind,:);

        while true
            interactions=GetPriorList(structure, list, expData);
            if isempty(interactions) || sum(interactions.Tested==inf)==0
                break;
            end
            [ layerStructure, layerParameters, interactions, expData]=GenerateLayer(interactions, expData);
            layer=layer+1;
            structure=[structure; layerStructure];
            parameters=[parameters; layerParameters];
            interactions.ExpansionIteration(interactions.ExpansionIteration==inf & interactions.Tested<chi2inv(0.95,dgf))=layer;
            list(ismember(list(:,1:2),interactions(interactions.Tested~=inf,1:2),'rows'),:)=interactions(interactions.Tested~=inf,:);
            added=interactions(interactions.Tested<chi2inv(0.95,dgf),:);
            if ~isempty(layerStructure)
                save(['Results/layer' num2str(layer) '.mat'])
            else
                save(['Results/empty' num2str(layer) '.mat'])
            end
        end
    end

end

