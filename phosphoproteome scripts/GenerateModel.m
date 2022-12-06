function [  ] = GenerateModel( interactions, modelName, folderLocation, maxAmount)
%GENERATEMODEL This model automatically generates IQM tools compatible
%models. Inputs are: interactions - a table of interactions, modelName -
%the desired name of the model, folderLocation - desired location.
%giving name/location=[] uses default values (model / current directory)
%   Detailed explanation goes here

%% Handle inputs
if nargin<1
    disp 'Need at least 1 input.'
    return
elseif nargin<2
    disp 'Assuming model name is \'model\' '
    modelName='model';
    disp 'Assuming folder should be pwd/ModelFiles'
    folderLocation=pwd;
    maxAmount=100;
elseif nargin<3
    disp 'Assuming folder should be pwd/ModelFiles'
    folderLocation=pwd;
    maxAmount=100;
elseif nargin<4
    maxAmount=100;
end

if isempty(modelName)
    modelName='model';
end

if isempty(folderLocation)
    folderLocation=pwd;
else
    warning off
    mkdir(folderLocation)
    warning on
end

interactions{:,1:2}=strrep(interactions{:,1:2},'-','');


%% Generate states
states=cell2table(repmat({''},1,2),'VariableNames',{'State','Reaction'});
states(ismember(states.State,''),:)=[];

reactions=cell2table(repmat({''},1,2),'VariableNames',{'Reaction','Equation'});
reactions(ismember(reactions.Reaction,''),:)=[];
variables=cell2table(repmat({''},1,2),'VariableNames',{'Variable','Equation'});
variables(ismember(variables.Variable,''),:)=[];
parameters={}; %=[];
interpParam=[];

% loop over all interactions
for i=1:height(interactions)
    state=cell2table(cell(2,2),'VariableNames',{'State','Reaction'});
    reaction=cell2table(cell(1,2),'VariableNames',{'Reaction','Equation'});
    variable={};
    params={};

    try
        type=interactions.Type{i};
        Source=interactions.Source{i};
        if ~isempty(Source) && ~isempty(regexp(Source(1),'[0-9]'))
            Source=['_' Source];
        end

        from=interactions.Target{i};
        if ~isempty(from) && ~isempty(regexp(from(1),'[0-9]'))
            from=['_' from];
        end

        if sum(strcmp(interactions.Site{i},interactions.Target))>0
            to=interactions.Site{i};
        else
            to=[interactions.Target{i} '_' interactions.Site{i}];
        end
        if ~isempty(to) && ~isempty(regexp(to(1),'[0-9]'))
            to=['_' to];
        end
    catch err
        disp(getReport(err));

    end
    reaction.Reaction=['v_' from '__' to];

    % Check what type of interaction, and perform appropriate actions.
    if strfind(type,'Interp')
        disp('Found interpolation.')
        from='';
        to='';
        reaction=cell2table(repmat({''},1,2),'VariableNames',{'Reaction','Equation'});
        reaction(ismember(reaction.Reaction,''),:)=[];
        variable=cell2table(cell(1,2),'VariableNames',{'Variable','Equation'});
        n=str2double(type(7:end));

        x=strrep(strcat({[interactions.Target{i} 'x' ]},num2str((1:n)')),' ','');
        y=strrep(strcat({[interactions.Target{i} 'y' ]},num2str((1:n)')),' ','');
        variable.Variable=interactions.Target(i);
        variable.Equation={['(1-useInterp)+useInterp*interpcsIQM([' strjoin(x,', ') '], [' strjoin(y,', ') '], time)']};
        interpParam=[interpParam; x;y];
    elseif any(regexp(type,'PO[M]*'))
        disp('Found irreversible phosphorylation')
        params={['k_' from '_' to]};
        duplicateIndex=find(ismember(params, parameters));

        for j=1:length(duplicateIndex)
            duplicateCount=sum(~cellfun(@isempty,regexp(parameters,strcat(params{j}, '_\d'))));
            params(duplicateIndex(j))=strcat(params(duplicateIndex(j)), '_', num2str(duplicateCount+1));
        end
        reaction.Equation=['-' from '*' Source '*' params{1}];
    elseif any(regexp(type,'P[M]*')) || any(regexp(type,'PR[M]*'))
        disp('Found reversible phosphorylation')
        params={['k_' from '_' to '_f']; ['k_' from '_' to '_b']};
        duplicateIndex=ismember(params, parameters);
        if sum(duplicateIndex>0)
            params(duplicateIndex)=strcat(params(duplicateIndex), '_1');
        end
        reaction.Equation=['-'  from '*' Source '*' params{1} '  +  ' to '*' params{2}];
    elseif any(regexp(type,'IO[M]*'))
        disp('Found irreversible inhibition')
        params={['k_' from '_' to]};
        duplicateIndex=ismember(params, parameters);
        for j=1:length(duplicateIndex)
            duplicateCount=sum(~cellfun(@isempty,regexp(parameters,strcat(params{j}, '_\d'))));
            params(duplicateIndex(j))=strcat(params(duplicateIndex(j)), '_', num2str(duplicateCount+1));
        end
        reaction.Equation=['-'  from '/' Source '*' params{1}];
    elseif any(regexp(type,'I[M]*')) ||any(regexp(type,'IR[M]*'))
        disp('Found reversible phosphorylation')
        params={['k_' from '-' to '_f']; ['k_' from '_' to '_b']};
        duplicateIndex=ismember(params, parameters);
        for j=1:length(duplicateIndex)
            duplicateCount=sum(~cellfun(@isempty,regexp(parameters,strcat(params{j}, '_\d'))));
            params(duplicateIndex(j))=strcat(params(duplicateIndex(j)), '_', num2str(duplicateCount+1));
        end
        reaction.Equation=['-'  from '/' Source '*' params{1} '  +  ' to '*' params{2}];
    elseif any(regexp(type,'R[M]*'))
        disp('Found reversible interaction')
        params={['k_' from '_' to '_f']; ['k_' from '_' to '_b']};
        duplicateIndex=ismember(params, parameters);
        for j=1:length(duplicateIndex)
            duplicateCount=sum(~cellfun(@isempty,regexp(parameters,strcat(params{j}, '_\d'))));
            params(duplicateIndex(j))=strcat(params(duplicateIndex(j)), '_', num2str(duplicateCount+1));
        end
        reaction.Equation=['-'  from '*' params{1} '  +  ' to '*' params{2}];
    elseif any(regexp(type,'DO[M]*'))
        disp('Found irreversible dephosphorylation')
        params={['k_' from '_' to]};
        duplicateIndex=find(ismember(params, parameters));

        for j=1:length(duplicateIndex)
            duplicateCount=sum(~cellfun(@isempty,regexp(parameters,strcat(params{j}, '_\d'))));
            params(duplicateIndex(j))=strcat(params(duplicateIndex(j)), '_', num2str(duplicateCount+1));
        end
        reaction.Equation=[to '*' Source '*' params{1}];
    elseif any(regexp(type,'D[M]*'))
        disp('Found reversible dephosphorylation')
        params={['k_' from '_' to '_f']; ['k_' from '_' to '_b']};
        duplicateIndex=ismember(params, parameters);
        for j=1:length(duplicateIndex)
            duplicateCount=sum(~cellfun(@isempty,regexp(parameters,strcat(params{j}, '_\d'))));
            params(duplicateIndex(j))=strcat(params(duplicateIndex(j)), '_', num2str(duplicateCount+1));
        end
        reaction.Equation=['-'  from  '*' params{1} '  +  ' to '*' Source '*'  params{2}];
    elseif any(regexp(type,'O[M]*'))
        disp('Found irreversible interaction')
        params={['k_' from '_' to]};
        duplicateIndex=ismember(params, parameters);
        for j=1:length(duplicateIndex)
            duplicateCount=sum(~cellfun(@isempty,regexp(parameters,strcat(params{j}, '_\d'))));
            params(duplicateIndex(j))=strcat(params(duplicateIndex(j)), '_', num2str(duplicateCount+1));
        end
        reaction.Equation=['-' from '*' params{1}];
    else
        disp('Unknown interaction type. Assuming reversible MA.');
    end


    if any(regexp(type,'P\w*[M]'))
        params=[params; {['km_' from '_' to]}];

        duplicateIndex=find(ismember(params, parameters));
        for j=1:length(duplicateIndex)
            duplicateCount=sum(~cellfun(@isempty,regexp(parameters,strcat(params{j}, '_\d'))));
            params(duplicateIndex(j))=strcat(params(duplicateIndex(j)), '_', num2str(duplicateCount+1));
        end
        reaction.Equation=regexprep(reaction.Equation,Source,[ Source '/(' params{end} ' + ' Source ')']);
    end


    state.State={from; to};
    state.Reaction={reaction.Reaction;strcat('-', reaction.Reaction)};
    state(ismember(state.State,''),:)=[];

    if ~isempty(reaction.Reaction)
        if isempty(reactions) || ~ismember(reaction.Reaction,reactions.Reaction)
            reactions=[reactions; reaction];
        else
            [~,ind]=ismember(reaction.Reaction,reactions.Reaction);
            reactions.Equation{ind}=[reactions.Equation{ind} ' + ' reaction.Equation];

        end
    end

    % Handle new states
    newIndex=~ismember(state.State,states.State);
    states=[states; state(newIndex,:)];
    variables=[variables; variable];
    parameters=[parameters; params];

    % Handle existing states
    [~,ind]=ismember(state.State(~newIndex),states.State);
    for j=1:length(ind)
        if ~isequal(states.Reaction{ind(j)},reaction.Reaction) && ~isequal( states.Reaction{ind(j)},['-' reaction.Reaction])
            reactionString=state.Reaction{ismember(state.State,states.State{ind(j)})};
            if reactionString(1)~='-'
                reactionString=[' + ' reactionString];
            else
                reactionString=[' ' reactionString(1) ' ' reactionString(2:end)];
            end
            states.Reaction{ind(j)}=[states.Reaction{ind(j)} reactionString];
        end
    end
end % end interactions loop

%% Print to File
fid=fopen([modelName '.txt'],'wt');

% Print header
fprintf(fid,'********** MODEL NAME\n');
fprintf(fid,[modelName '\n\n']);
fprintf(fid,'********** MODEL NOTES\n');
fprintf(fid,'\n');

% Print states
fprintf(fid,'********** MODEL STATES\n');
for i=1:height(states)
    fprintf(fid,['d/dt(' states.State{i} ')=' states.Reaction{i} '\n']);
end
fprintf(fid,'\n');

for i=1:height(states)
    usedStates=strrep(states.Reaction{i},'-','');
    usedStates=strrep(usedStates,'v_','');
    usedStates=strrep(strrep(usedStates, '+' ,' '),'__',' ');

    usedStates=unique(strsplit(usedStates));

    fprintf(fid,[states.State{i} '(0)=' num2str(maxAmount/length(usedStates)) '\n']);
end
fprintf(fid,'\n');

% Print parameters
fprintf(fid,'********** MODEL PARAMETERS\n');

if sum(~cellfun(@isempty,regexpi(interactions.Type,'interp')))>0
    for i=1:length(interpParam)
        fprintf(fid,[interpParam{i} '=1\n']);
    end
    fprintf(fid,'\n');
    fprintf(fid, 'useInterp=1\n');
    fprintf(fid,'\n');
end
for i=1:length(parameters)
    fprintf(fid,[parameters{i} '=1\n']);
end
fprintf(fid,'\n');


% Print variables
fprintf(fid,'********** MODEL VARIABLES\n');
for i=1:height(variables)
    fprintf(fid,[variables.Variable{i} '=' variables.Equation{i} '\n']);
end
fprintf(fid,'\n');

% Print reactions
fprintf(fid,'********** MODEL REACTIONS\n');
for i=1:height(reactions)
    fprintf(fid,[reactions.Reaction{i} '=' reactions.Equation{i} '\n']);
end
fprintf(fid,'\n');

% Print footer
fprintf(fid,'********** MODEL FUNCTIONS\n\n');
fprintf(fid,'********** MODEL EVENTS\n\n');
fprintf(fid,'********** MODEL MATLAB FUNCTIONS\n\n');

fclose(fid);
if ~strcmp(folderLocation,pwd)
    movefile([modelName '.txt'],folderLocation,'f')
end
end

