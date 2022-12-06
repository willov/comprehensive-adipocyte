function [] = Setup_combined(modelName)
if nargin<1, modelName=[]; end

if ~exist('IQMsimulate','file')
    fprintf('\n\nFor these scripts to work, the IQM tools toolbox have to be compiled.\n')
    disp('To do this, a valid C-compiler is necessary')
    disp('If the compilation of the toolbox does not, it is likely that a valid C-compiler is missing.')
    input('Press any key to continue and compile the toolbox.\nPress ctrl + C to abort');
    run('../Tools/IQM Tools/installIQMtoolsInitial.m')
end

addpath('../Tools')
addpath('../Data/')
addpath(genpath('../Models'))
addpath(genpath('.'))

%% Setup the model

%% Compile the models
if ~isempty(modelName)
    cwd = pwd;
    try
        cd('../Models');
        IQMmakeMEXmodel(IQMmodel([modelName,'.txt']));
        cd(cwd)
    catch
        cd(cwd)
    end
end

end