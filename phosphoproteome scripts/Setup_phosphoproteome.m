function [] = Setup_phosphoproteome()

if ~exist('IQMsimulate','file')
    fprintf('\n\nFor these scripts to work, the IQM tools toolbox have to be compiled.\n')
    disp('To do this, a valid C-compiler is necessary')
    disp('If the compilation of the toolbox does not, it is likely that a valid C-compiler is missing.')
    input('Press any key to continue and compile the toolbox.\nPress ctrl + C to abort');
    run('../Tools/IQM Tools/installIQMtoolsInitial.m')
end

addpath('../Tools')
addpath(genpath('../Data/'))
addpath(genpath('../Models'))
addpath(genpath('../combined model scripts'))
addpath(genpath('.'))

if ~exist("../Phosphoproteome-plots",'dir')
    mkdir("../Phosphoproteome-plots")
end

%% Compile the models
cwd = pwd;
try
    cd('../Models');
    files=dir('./*.txt');
    for i = 1:length(files)
        IQMmakeMEXmodel(IQMmodel(files(i).name))
    end

    cd('Pairwise')
    files=dir('./*.txt');
    for i = 1:length(files)
        IQMmakeMEXmodel(IQMmodel(files(i).name))
    end

    cd(cwd)
catch
    cd(cwd)
end
end