function [expData, inputData, Bmarkinterp, diabeticInput, inhibitionData]=GetInputs(expData, x, coreModel, simulateInhibitions, range)
load('Translation.mat','Translation')
load('core.mat','core')

if nargin<4
    simulateInhibitions=false;
    range=0:5;
elseif nargin <5
    range=0:5;
end

[normalSim, diabeticSim]=CoreSimulation(coreModel);

% Extract and format the correct variables
inputData=table(normalSim.variables','VariableNames',{'Target'});
inputData.simulatedValues=normalSim.variablevalues(ismember(normalSim.time,x),:)';
inputData=outerjoin(inputData, Translation, 'LeftKey','Target','RightKey','BmarkName','Type','Left','MergeKeys',true);
inputData(:,1)=[];

diabeticInput=table(diabeticSim.variables','VariableNames',{'Target'});
diabeticInput.simulatedValues=diabeticSim.variablevalues(ismember(diabeticSim.time,x),:)';
diabeticInput=outerjoin(diabeticInput, Translation, 'LeftKey','Target','RightKey','BmarkName','Type','Left','MergeKeys',true);
diabeticInput(:,1)=[];

if ~ismember(expData.Properties.VariableNames, 'simulatedValues') % Add the simulated values to the table of experimental data.
    expData(ismember(expData.Index, core.Index),:)=[]; % Remove all sites corresponding to the core data from the experimental data.
    expData=[expData; core];
    expData.Genenamesprimary = strrep(strcat(expData.Genenamesprimary, '_', expData.PhosphorylatedAminoAcid, num2str(expData.IPI368PositioninProtein)),' ', '');
    expData=outerjoin(expData, inputData, 'LeftKey','Genenamesprimary','RightKey','HumphreyName','Type','Left');
    expData(:,end)=[];
end

expData(ismember(expData.Genenamesprimary, {'Mapk1_Y185', 'Mapk3_T203', 'Rps6_S236'}),:)=[]; %Remove entries in the data originating from the core model, which have the same simulations in other entry.

Bmarkinterp=table(repmat({''},height(Translation),1),strcat(Translation.HumphreyName,'_p'),repmat({''},height(Translation),1),repmat({'Interp9'},height(Translation),1),...
    'VariableNames',{'Source','Target','Site','Type'});
Bmarkinterp=sortrows(Bmarkinterp,'Target');

if simulateInhibitions
    inhibitionData=table(normalSim.variables','VariableNames',{'Target'});

    for n=range
        [inhibitionSimulation, ~]=CoreSimulation(coreModel,[1/(2^n) 1]);
        inhibitionData=[inhibitionData table(inhibitionSimulation.variablevalues(ismember(normalSim.time,x),:)','VariableNames',{['MK' num2str(round((1-1/(2^n))*100))]}) ];
    end
    for n=range
        [inhibitionSimulation, ~]=CoreSimulation(coreModel,[1/(2^n) 1/(2^n)]);
        inhibitionData=[inhibitionData table(inhibitionSimulation.variablevalues(ismember(normalSim.time,x),:)','VariableNames',{['LY' num2str(round((1-1/(2^n))*100))]}) ];
    end
    inhibitionData=outerjoin(inhibitionData, Translation, 'LeftKey','Target','RightKey','BmarkName','Type','Left','MergeKeys',true);
    inhibitionData(:,1)=[];
else
    inhibitionData=[];
end

end