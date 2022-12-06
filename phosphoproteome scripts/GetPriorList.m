function [Prior]=GetPriorList(structure, SuperPrior, matchedData)

%% This function takes a structure, an interaction list and a list of data as input.
%  As output it sends a list of interactions that should be tested.

if ~any([isempty(structure), isempty(SuperPrior), isempty(matchedData)])
    list=structure; %just so we dont overwrite what we may want to keep.

    list.Source=regexprep(list.Source,'_[pPuU]+','');
    list.Target=regexprep(list.Target,'_[pPuU]+','');
    interpIndex=find(~cellfun('isempty',strfind(list.Type,'Interp')));
    list.Source(interpIndex)=list.Target(interpIndex);
    list.Target(interpIndex)=repmat({''},length(interpIndex),1);

    interactionsWithData=SuperPrior(ismember(SuperPrior.Source,matchedData.Genenamesprimary) & ismember(SuperPrior.Target,matchedData.Genenamesprimary),:); % Only use interactions that have data
    interactionsWithData(strcmp(interactionsWithData.Source,interactionsWithData.Target),:)=[]; %Remove all autophosphorylations
    Prior=interactionsWithData(ismember(interactionsWithData.Source,[list.Target, list.Source]),:); % add all possible new interactions based on present Source/Targets.
    Prior(ismember(Prior.Target, [list.Target, list.Source]),:)=[]; % Removes all interactions, where the Target is already added.

    Prior=sortrows(Prior,'Tested','ascend');
    [~,ind]=unique(Prior(:,1:2),'rows');
    Prior=Prior(ind,:);
else
    Prior=[];
end
end