function [] = PlotFigures_combined()

% This function will plot all figures corresponding to the combined core
% model.
%% Setup initial things
close all
warning off

choice = input("\n\nDo you want to simulate dose responses in high resolution (as in the paper), \nor faster simulations with a lower resolution? \nPress ctrl+C to cancel. \nHigh resolution? (y/N): ", 's');
if lower(choice) == 'y'
    disp('Using the original resolution (slow)')
    res = 0.01; %use 0.01 for the same high resolution from the original paper
else
    disp('Using a lower resolution (faster)')
    res = 1; % 1 for a lower resolution that decreases the time to simulate
end

choice = input("\n\nDo you want to simulate the supplemental figures in addition to the main figures? \npress ctrl+C to cancel. \nSimulate supplemental figures (y/N): ", 's');
plotSupplemental =  lower(choice) == 'y';

format shortg
fcost = fopen('costs.txt','w');

%% Setup the model
modelName = 'combined';
run('combined model scripts/Setup_combined(modelName)')

%% Setup values corresponding to the estimation and testing datasets
useCL_ATP=0;
useHSL = 0;
baseFolder = './combined model scripts/Results/Without HSL and CL_ATP';

[model,data, ~, ~, ~, expInd, stimLipo, stimAdi, ~, dgf] = Init_combined(modelName, useHSL, useCL_ATP);

%% Calculate the cost for the estimation dataset
load('./combined model scripts/Results/Without HSL and CL_ATP/combined, opt(616.501639) 220204-184028.mat', 'optParam')

limitEstimation = chi2inv(0.95, dgf(1));
costEstimation = CostCombined(optParam,model, expInd, data, stimLipo, stimAdi, limitEstimation, useHSL);

fprintf('Cost-estimation: %.2f, limit: %.2f (dgf=%i), pass: %d\n', costEstimation, limitEstimation, dgf(1), costEstimation<limitEstimation)
fprintf(fcost, 'Cost-estimation: %.2f, limit: %.2f (dgf=%i), pass: %d\n', costEstimation, limitEstimation, dgf(1), costEstimation<limitEstimation);

%% Calculate the cost for the test data set
[model,data, ~, ~, ~, expInd, stimLipo, stimAdi,stimAdiPred, dgf] = Init_combined(modelName, useHSL, useCL_ATP);
toPredict.submodel='lipoadi';
toPredict.polarity='min';

load('./combined model scripts/Results/Without HSL and CL_ATP/combined, opt_test(20.6510) 220401-124306.mat', 'optParam')
[costTest, costEstimationTest] = CostCombinedPred(optParam,model, expInd, data, stimLipo, stimAdi, limitEstimation, useHSL, stimAdiPred, toPredict);

limitTest = chi2inv(0.95, dgf(2));

assert(costEstimationTest <= limitEstimation, 'The cost calculated with the two functions are not the same, something has gone wrong.')
fprintf('Cost-test: %.2f, limit: %.2f (dgf=%i), pass: %d\n', costTest, limitTest, dgf(2), costTest<limitTest)
fprintf(fcost, 'Cost-test: %.2f, limit: %.2f (dgf=%i), pass: %d\n', costTest, limitTest, dgf(2), costTest<limitTest);

%% Plot the agreements to the estimation and test data (supplementary)
if plotSupplemental
    fprintf("\n\nPlotting the model with uncertainty for the supplemental figues\n")
    PlotAgreementGlu(optParam, model, expInd, data.Rajan, baseFolder);
    PlotAgreementLipo(optParam, modelName, res, useHSL, 0, baseFolder)
    PlotAdi(optParam, modelName, expInd, 1,baseFolder,useCL_ATP);

    exportgraphics(figure(61), 'S1. Estimation-glucose.pdf', 'ContentType','vector')
    exportgraphics(figure(11), 'S2. Estimation-lipolysis.pdf', 'ContentType','vector')
    exportgraphics(figure(41), 'S3. Estimation-adiponectin.pdf', 'ContentType','vector')

    PlotInhibitorsGlu(optParam, model, expInd, data.Rajan);
    PlotAgreementLipo(optParam, modelName, res, useHSL, 1, baseFolder)
    PlotAdi(optParam, modelName, expInd, 2, baseFolder, useCL_ATP);

    exportgraphics(figure(53), 'S4. Test-lipoadi.pdf', 'ContentType','vector')
    exportgraphics(figure(62), 'S5. Test-glucose-rajan.pdf', 'ContentType','vector')
    exportgraphics(figure(63), 'S6. Test-glucose-brännmark.pdf', 'ContentType','vector')

    close all
end
%% Setup values corresponding to the full datasets
useCL_ATP=1;
useHSL = 1;
baseFolder = './combined model scripts/Results/With HSL and CL_ATP' ;

[model,data, ~, ~, ~, expInd, stimLipo, stimAdi,~, dgf] = Init_combined(modelName, useHSL, useCL_ATP);

limit = chi2inv(0.95, dgf(1));

%% Calculate the cost for the full dataset
load('./combined model scripts/Results/With HSL and CL_ATP/combined, opt(654.812129) 220214-162036.mat', 'optParam')

costCombined = CostCombined(optParam,model, expInd, data, stimLipo, stimAdi, limit, useHSL);
fprintf('Cost-total: %.2f, limit: %.2f (dgf=%i), pass: %d\n', costCombined, limit, dgf(1), costCombined<limit)
fprintf(fcost, 'Cost-total: %.2f, limit: %.2f (dgf=%i), pass: %d\n', costCombined, limit, dgf(1), costCombined<limit);

%% Plot the agreements to the extended data set
fprintf("\n\nPlotting the model with uncertainty for the main figues\n")

PlotIsoEffectGlu(optParam, model, expInd, data.Rajan)

PlotAgreementGlu(optParam, model, expInd, data.Rajan, baseFolder);
PlotAgreementLipo(optParam, modelName, res, useHSL, 0, baseFolder)
PlotAdi(optParam, modelName, expInd, 1,baseFolder, useCL_ATP);

exportgraphics(figure(11), 'Fig. 3 Estimation-lipolysis.pdf', 'ContentType','vector')
exportgraphics(figure(41), 'Fig. 4 Estimation-adiponectin.pdf', 'ContentType','vector')
exportgraphics(figure(61), 'Fig. 5 Estimation-glucose.pdf', 'ContentType','vector')

%% Final cleanup
close all
warning on
disp('Finished plotting the figures. The figures are available as PDFs in the root folder.')
disp('The costs have been printed during the execution, and are also printed to the "costs.txt" file')
fclose(fcost);
end