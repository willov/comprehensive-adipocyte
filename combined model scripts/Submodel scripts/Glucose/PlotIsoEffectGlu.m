function [] = PlotIsoEffectGlu(optParam, model, expInd)

colorNorm = [2 64 167]/256;
colorDiab = [0.6350 0.0780 0.1840];

fprintf('\nSimulating the effect of iso on glucose uptake.\n')
optParam(expInd)=exp(optParam(expInd));

%% Simulate with no iso
[glucoseDR, ins_conc] = SimulateGlucoseDR(optParam, model, 0);
GLUCOSE_dr_n(:,1) = glucoseDR(:,1);
GLUCOSE_dr_d(:,1) = glucoseDR(:,2);

%% Simulate with 10 nM (0.01 µM) iso
glucoseDR_iso = SimulateGlucoseDR(optParam, model, 0.01);
GLUCOSE_dr_n(:,2) = glucoseDR_iso(:,1);
GLUCOSE_dr_d(:,2) = glucoseDR_iso(:,2);

% normalize to fold over basal. 
GLUCOSE_dr_d = GLUCOSE_dr_d./GLUCOSE_dr_n(end,end)*100;
GLUCOSE_dr_n = GLUCOSE_dr_n./GLUCOSE_dr_n(end,end)*100;

figure(64)
hold on
ins_conc = ins_conc*1e-9; % Scaling from nM to M

plot(ins_conc,GLUCOSE_dr_n(:,1),'-', 'color', colorNorm, 'LineWidth',2)
plot(ins_conc,GLUCOSE_dr_n(:,2),'--', 'color', colorNorm, 'LineWidth',2)
plot(ins_conc,GLUCOSE_dr_d(:,1),'-', 'color', colorDiab, 'LineWidth',2)
plot(ins_conc,GLUCOSE_dr_d(:,2),'--', 'color', colorDiab, 'LineWidth',2)
legend({'normal, 0 nM iso', 'normal, 10 nM iso', 't2d, 0 nM iso',  't2d, 10 nM iso'})
axis('tight')
set(gca, 'XScale', 'log','YTick',[0 50 100],'XTick',[1e-11 1e-9 1e-7]);
box off
xlabel('[Insulin] (M)')
ylabel('% of max')
title({'Glucose uptake'})
set(gca,'FontSize',11)

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
    idxGLUCOSE   = find(strcmp(variables,'measuredGLUCOSE'));

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