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


[glucoseDR, ins_conc] = SimulateGlucoseDR(optParam, model, 0);

GLUCOSE_dr_n(:,1) = glucoseDR(:,1);
GLUCOSE_dr_d(:,1) = glucoseDR(:,2);

glucoseDR_iso = SimulateGlucoseDR(optParam, model, 0.01);

GLUCOSE_dr_n(:,2) = glucoseDR_iso(:,1);
GLUCOSE_dr_d(:,2) = glucoseDR_iso(:,2);

% normalize to fold over basal. 
GLUCOSE_dr_d = GLUCOSE_dr_d./GLUCOSE_dr_n(end,end)*100;
GLUCOSE_dr_n = GLUCOSE_dr_n./GLUCOSE_dr_n(end,end)*100;

figure(598)
hold on
ins_conc = ins_conc*1e-9; % Scaling from nM to M

% conc_glucose = expData.GLUCOSE_dr.conc;
plot(ins_conc,GLUCOSE_dr_n(:,1),'-', 'color', colorNorm, 'LineWidth',2)
plot(ins_conc,GLUCOSE_dr_d(:,1),'-', 'color', colorDiab, 'LineWidth',2)
% plot([conc_glucose(end)*7.5, conc_glucose(end)*12.5],[GLUCOSE_dr_n(1,1), GLUCOSE_dr_n(1,1)],'-', 'color', colorNorm, 'LineWidth',2)
% plot([conc_glucose(end)*7.5, conc_glucose(end)*12.5],[GLUCOSE_dr_d(1,1), GLUCOSE_dr_d(1,1)],'-', 'color', colorDiab, 'LineWidth',2)
plot(ins_conc,GLUCOSE_dr_n(:,2),'--', 'color', colorNorm, 'LineWidth',2)
plot(ins_conc,GLUCOSE_dr_d(:,2),'--', 'color', colorDiab, 'LineWidth',2)
legend({'normal, 0 iso', 't2d, 0 iso', 'normal, 10 nM iso', 't2d, 10 nmM iso'})
% errorbar(expData.GLUCOSE_dr.conc,expData.GLUCOSE_dr.response_n,expData.GLUCOSE_dr.sem_n,'o', 'color', colorNorm,'MarkerFaceColor','auto','MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
% errorbar(expData.GLUCOSE_dr.conc,expData.GLUCOSE_dr.response_d,expData.GLUCOSE_dr.sem_d,'o', 'color', colorDiab,'MarkerFaceColor','auto','MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
axis('tight')
set(gca, 'XScale', 'log','YTick',[0 50 100],'XTick',[1e-11 1e-9 1e-7]);
box off
xlabel('[Insulin] (M)')
ylabel('% of max')
title({'Glucose uptake'}) %, 0.05 nM ins.
set(gca,'FontSize',11)



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
conc = 1e-9*[0.001 0.01 0.1 0.3 1 10 100];
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


function [glucoseDR, ins_conc] = SimulateGlucoseDR(param, model, iso)
    simOptions              = [];
    simOptions.method       = 'stiff';
    simOptions.maxnumsteps  = 1e6;
    simOptions.reltol       = 1e-12;
    simOptions.abstol       = 1e-9;
    simOptions.minstep      = 1e-10;

    paramNames=IQMparameters(model);
    param(ismember(paramNames,{'kLdrift', 'kLclear'}))=0; % Disables the drift in experimental data and clearance in vivo of glycerol/FA
    index_diabetes = contains(paramNames,'diab')==1;
    index_gluc = ismember(paramNames,'gluc')==1;

    variables = IQMvariables(model);
    idxIR        = find(strcmp(variables,'measuredIR'));        
    idxIRS1      = find(strcmp(variables,'measuredIRS1'));     
    idxIRS1307   = find(strcmp(variables,'measuredIRS1307'));
    idxPKB308    = find(strcmp(variables,'measuredPKB308'));
    idxPKB473    = find(strcmp(variables,'measuredPKB473'));
    idxAS160     = find(strcmp(variables,'measuredAS160'));
    idxS6        = find(strcmp(variables,'measuredS6'));
    idxGLUCOSE   = find(strcmp(variables,'measuredGLUCOSE'));
    idxERK       = find(strcmp(variables,'measuredERK'));
    idxElk1      = find(strcmp(variables,'measuredElk1'));
    idxFOXO      = find(strcmp(variables,'measuredFOXO'));

    tempParam = [param 0 0 0 0 0 0 0]; % Inserting values used in other models

    sim=struct();

    %% Simulate steady state
    param_n = tempParam;                                % Parameter values
        param_n(index_diabetes)     = 1;                    % No diabetes
        sim.n_0 = model(0:1:500, [], param_n, simOptions);  % baseline simulation controls

        param_d = tempParam;                                % Parameter values for diabetics
        
        initcond_0_d = sim.n_0.statevalues(end,:);          % Initial conditions
        index_IR = ismember(sim.n_0.states,'IR')==1;        % index for state IR
        initcond_0_d(index_IR) = 55;                        % Diabetics IR amount 55 % of controls
        index_GLUT4 = ismember(sim.n_0.states,'GLUT4')==1;  % index for state GLUT4
        initcond_0_d(index_GLUT4) = 50;                     % Diabetics IR amount 55 % of controls
        
        index_FOXO = ismember(sim.n_0.states,'FOXO')==1;
        initcond_0_d(index_FOXO) = 55;                      
        index_AS160 = ismember(sim.n_0.states,'AS160')==1;
        initcond_0_d(index_AS160) = 45; 
        index_S6 = ismember(sim.n_0.states,'S6')==1;
        initcond_0_d(index_S6) = 48;

        sim.d_0 = model([0 490 500], initcond_0_d, param_d, simOptions);  % baseline simulation diabetics
       
        initcond_n = sim.n_0.statevalues(end,:);    % using the steady state values from the baseline-simulation
        initcond_d = sim.d_0.statevalues(end,:);  

    %% Set iso, glucose and ins doses
    param_n(index_gluc)=0.05; % Stimulating with 0.05 mM glucose
    param_d(index_gluc)=0.05; % Stimulating with 0.05 mM glucose

    param_n(end-3) = iso;
    param_d(end-3) = iso;

    ins_conc = 10.^linspace(log10(0.001), log10(100), 100);
    glucoseDR = nan(length(ins_conc),2);

    %% Simulate dose response experiments, with increasing doses of insulin.
    for idx = 1:length(ins_conc)
        % Simulate normal conditions
        param_n(end) = ins_conc(idx);
        sim_normal = model([0, 30], initcond_n, param_n, simOptions);   % 30 min simulation with 0.01 nM insulin
        glucoseDR(idx,1) = sim_normal.variablevalues(end,idxGLUCOSE);

        % Simulate diabetic conditions
        param_d(end) = ins_conc(idx);
        sim_diabetes = model([0, 30], initcond_d, param_d, simOptions);   % 30 min simulation with 0.01 nM insulin
        glucoseDR(idx,2) = sim_diabetes.variablevalues(end,idxGLUCOSE);
    end

end