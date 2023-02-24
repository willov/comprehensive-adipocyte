function [] = PlotIsoEffectGlu(optParam, model, expInd, expData)

colorNorm = [2 64 167]/256;
colorDiab = [0.6350 0.0780 0.1840];

time=unique([expData.AllTimes 0:0.01:60]);

variables = IQMvariables(model);
idxGLUCOSE   = strcmp(variables,'measuredGLUCOSE');

fprintf('\nSimulating the effect of iso on glucose uptake.\n')
isoDoses = [0.0001, 0.001, 0.01, 0.1, 1];

optParam(expInd)=exp(optParam(expInd));

%% Simulate with no iso
[sim, tmpCost] = SimulateGlucose(optParam, model, time);
idxT30=find(sim.n_60.time==30);

GLUCOSE = [sim.dr_gluc_n0.variablevalues(end,idxGLUCOSE); sim.dr_gluc_n6.variablevalues(idxT30,idxGLUCOSE)];
scaleGLUCOSE = lscov(GLUCOSE,expData.GLUCOSE_2p.response_n);
GLUCOSE_n(:,1) = scaleGLUCOSE*[sim.dr_gluc_n0.variablevalues(end,idxGLUCOSE); sim.dr_gluc_n5.variablevalues(idxT30,idxGLUCOSE)];
GLUCOSE_d(:,1) = scaleGLUCOSE*[sim.dr_gluc_d0.variablevalues(end,idxGLUCOSE); sim.dr_gluc_d5.variablevalues(idxT30,idxGLUCOSE)];

%% Simulate with iso doses
n_doses = length(isoDoses);
for i = 1:n_doses


    [sim, tmpCost] = SimulateGlucose(optParam, model, time, isoDoses(i));
    if tmpCost >= 1e20
        disp('Simulation has crashed.');
        return
    end
    %% Cost calculations
    % Glucose uptake two point data
    % 0 and 100 nM insulin 15 min
    % then 0.05 nM glucose 30 min
    GLUCOSE_n(:,i+1) = scaleGLUCOSE*[sim.dr_gluc_n0.variablevalues(end,idxGLUCOSE); sim.dr_gluc_n5.variablevalues(idxT30,idxGLUCOSE)];
    GLUCOSE_d(:,i+1) = scaleGLUCOSE*[sim.dr_gluc_d0.variablevalues(end,idxGLUCOSE); sim.dr_gluc_d5.variablevalues(idxT30,idxGLUCOSE)];

    if i==1
        fprintf('%i of %i \n|',i,n_doses)
    elseif mod(i,50)==0
        fprintf(' %i of %i \n',i,n_doses)
    else
        fprintf('|')
    end
end
fprintf('\n')

%% Normalize to fold over normal 0 ins + 0 iso
GLUCOSE_d = GLUCOSE_d/GLUCOSE_n(1,1);
GLUCOSE_n = GLUCOSE_n/GLUCOSE_n(1,1);

%% Plot figures
figsize = [450 75 900 600];

isoDoses = [0 isoDoses];
n_doses = length(isoDoses);

set(figure(599),'Position', figsize);
subplot(2,1,1)
bar(1:n_doses,  GLUCOSE_n(1,:), 'facecolor', colorNorm)
hold on
bar((1:n_doses)+n_doses+2,GLUCOSE_n(2,:), 'facecolor', colorNorm)
xticks([1:n_doses (1:n_doses)+n_doses+2])
xticklabels([string(isoDoses) string(isoDoses)])
title('Glucose uptake, normal')
ylabel('Fold over normal basal')
yticks([0 1 2])
text(3,-0.3,'0 nM ins','FontSize',18,'FontWeight','Bold')
text(11,-0.3,'10 nM ins.','FontSize',18,'FontWeight','Bold')
text(0,-0.1,'iso (nM):', 'FontSize',12)
text(8,-0.1,'iso (nM):', 'FontSize',12)
box off
set(gca,'FontSize',11)
y_limits = get(gca,'ylim');

subplot(2,1,2)
bar(1:n_doses,  GLUCOSE_d(1,:), 'facecolor', colorDiab)
hold on
bar((1:n_doses)+n_doses+2,GLUCOSE_d(2,:), 'facecolor', colorDiab)
xticks([1:n_doses (1:n_doses)+n_doses+2])
xticklabels([string(isoDoses) string(isoDoses)])
yticks([0 1 2])
title('Glucose uptake, type 2 diabetes')
ylabel('Fold over normal basal')
text(3,-0.3,'0 nM ins','FontSize',18,'FontWeight','Bold')
text(11,-0.3,'10 nM ins.','FontSize',18,'FontWeight','Bold')
text(0,-0.1,'iso (nM):', 'FontSize',12)
text(8,-0.1,'iso (nM):', 'FontSize',12)
box off
set(gca,'FontSize',11, 'YLim', y_limits)
set(figure(599), 'outerposition',[249 1 1727 1408], 'PaperType','a4') % Vertical mode

end

