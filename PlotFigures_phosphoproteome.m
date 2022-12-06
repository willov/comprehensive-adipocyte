function [] = PlotFigures_phosphoproteome()

% This function will plot all figures corresponding to the automatic model
% expansion

%% Setup initial things
close all
warning off

choice = input("\n\nDo you want to simulate the supplemental figures in addition to the main figures? \nThis will take a long time. \nPress ctrl+C to cancel. \nSimulate supplemental figures (y/N): ", 's');
plotSupplemental =  lower(choice) == 'y';

format shortg
cwd = pwd;
cd('phosphoproteome scripts')
try
    %% Setup the model
    coreModelName = 'combined';
    Setup_phosphoproteome()

    %% Generate the layerwise inibition prediction plot.
    InhibitionsForAllLayers(coreModelName)

    %% Plot supplemental timeseries
    disp('Plotting diabetes time series for the first 5 layers.')
    load('Results/layer6.mat','x', 'list','expData','structure','parameters') % Note: the second added layer is empty, therefore the fifth layer correspond to the saved layer6.
    list.ExpansionIteration(list.ExpansionIteration>2) = list.ExpansionIteration(list.ExpansionIteration>2)-1; % Adjust the layer names since layer 2 is empty
    modelName='layer5';
    disp('Starting to build the model file.')
    GenerateModel(structure,modelName,pwd,100);
    disp('Starting to compile the model. This may take a while')
    IQMmakeMEXmodel(IQMmodel([modelName '.txt']))
    PlotFinalModel(parameters,x, expData, modelName, coreModelName, list, 1, 0,1,0);
    disp('Finished plotting the figures');
    disp('The main figures are available as PDFs in the root folder.')

    figFiles = dir('../Phosphoproteome-plots/Diabetes/*.pdf');
    for i =1:length(figFiles)
        copyfile([figFiles(i).folder '/' figFiles(i).name], ['../Fig. 9 ' figFiles(i).name])
    end
    rmdir('../Phosphoproteome-plots/Diabetes/','s')
    close all

    if plotSupplemental
        disp('Plotting all supplemental figures. This might take a while. ')
        load('Results/DataDriven.mat','x', 'list','expData','structure','parameters') % Load to restart the algorithm from this point later on if necessary
        modelName='FinalExpandedModel';
        GenerateModel(structure,modelName,pwd,100);
        IQMmakeMEXmodel(IQMmodel([modelName '.txt']))
        PlotFinalModel(parameters,x, expData, modelName, coreModelName, list, 1, 1,1,1);
        disp('Finished plotting the supplemental figures figures');
        disp('The supplemental figures are available as PDFs in the Phosphoproteome-plots folder.')
        disp('Note that some layers are empty, and thus no plots for that layer is saved (e.g., in layer 2).')
    end

    %% Final cleanup
    close all
    warning on
    cd(cwd)

catch err
    disp('Something went wrong, returning to the original directory')
    cd(cwd)
    rethrow(err)
end
end