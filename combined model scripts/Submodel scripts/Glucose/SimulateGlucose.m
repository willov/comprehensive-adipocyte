function [sim, cost] = SimulateGlucose(param, model, time, isoIn)
    if nargin<3
        time=0:0.01:60;
    end
    if nargin<4, isoIn=[]; end

    cost = 0;
    sim=struct();
    
    simOptions              = [];
    simOptions.method       = 'stiff';
    simOptions.maxnumsteps  = 1e6;
    simOptions.reltol       = 1e-12;
    simOptions.abstol       = 1e-9;
    simOptions.minstep      = 1e-10;
    
    paramNames=IQMparameters(model);
    param(ismember(paramNames,{'kLdrift', 'kLclear'}))=0; % Disables the drift in experimental data and clearance in vivo of glycerol/FA

    %% Indices
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
    
    idxT4=find(time==4);
    idxT30=find(time==30);
    idxT60=find(time==60);
    %% Baseline simulation
    
    % input from other models not used in this model
    pip0 = 0;
    phe0 = 0; % µM
    epi0 = 0; % µM
    iso0 = 0; % µM
    CL0 = 0;
    gluc0 = 0;
    ins0 = 0; % nM

    if length(param)==find(index_diabetes,1,'first') % If diab_reest is missing, add it. 
        param = [param 1]; 
    end

    tempParam = [param pip0 phe0 iso0 epi0 CL0 gluc0 ins0]; % Inserting values used in other models
    try
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
    catch err
        cost = 1e24;
        return
    end

    initcond_n = sim.n_0.statevalues(end,:);    % using the steady state values from the baseline-simulation
    initcond_d = sim.d_0.statevalues(end,:);  

    %% Simulate 30 minutes with no insulin input if iso is given
    if ~isempty(isoIn)
        param_n(end-3) = isoIn;
        param_d(end-3) = isoIn;
    end
    %% Insulin stimulation simulations
    
    ins_conc = [0.001 0.01 0.03 0.1 0.3 1 10 100];
    try
        % Controls



        sim.dr_n001 = model(time(1:idxT30), initcond_n, [param_n(1:end-1) ins_conc(1)], simOptions);   % 30 min simulation with 0.01 nM insulin
        sim.dr_n1a = model(time(1:idxT30), initcond_n, [param_n(1:end-1) ins_conc(2)], simOptions);   % 15 min simulation with 0.01 nM insulin
        sim.dr_n1b = model(time(1:idxT30), initcond_n, [param_n(1:end-1) ins_conc(3)], simOptions);   % 15 min simulation with 0.03 nM insulin
        sim.dr_n2 = model(time(1:idxT30), initcond_n, [param_n(1:end-1) ins_conc(4)], simOptions);    % 15 min simulation with 0.1 nM insulin
        sim.dr_n3 = model(time(1:idxT30), initcond_n, [param_n(1:end-1) ins_conc(5)], simOptions);    % 15 min simulation with 0.3 nM insulin
        sim.dr_n4 = model(time(1:idxT30), initcond_n, [param_n(1:end-1) ins_conc(6)], simOptions);    % 15 min simulation with 1 nM insulin
        sim.n_60 = model(time(1:idxT60), initcond_n, [param_n(1:end-1) ins_conc(7)], simOptions);     % 60 min simulation with 10 nM insulin
        sim.n_180 = model(time(1:idxT30), initcond_n, [param_n(1:end-1) ins_conc(8)], simOptions);   % 180 min simulation with 100 nM insulin %NOTE: used only for dosereponse at specific times, so 30 is fine
        sim.ds1_n = model(time(1:idxT4), initcond_n, [param_n(1:end-1) 1.2], simOptions);       % 4 min simulation with 1.2 nM insulin for double step data
        initcond2_n = sim.ds1_n.statevalues(end,:);                                 % using the statevalues from the ds1 simulation as initial conditions
        sim.ds2_n = model(time(1:idxT4), initcond2_n, [param_n(1:end-1) 10], simOptions);        % 4 min simulation with 10 nM insulin for double step data
        
        % Diabetics
        sim.dr_d001 = model(time(1:idxT30), initcond_d, [param_n(1:end-1) ins_conc(1)], simOptions);   % 30 min simulation with 0.01 nM insulin
        sim.dr_d1a = model(time(1:idxT30), initcond_d, [param_d(1:end-1) ins_conc(2)], simOptions);   % 15 min simulation with 0.01 nM insulin   
        sim.dr_d1b = model(time(1:idxT30), initcond_d, [param_d(1:end-1) ins_conc(3)], simOptions);   % 15 min simulation with 0.03 nM insulin
        sim.dr_d2 = model(time(1:idxT30), initcond_d, [param_d(1:end-1) ins_conc(4)], simOptions);    % 15 min simulation with 0.1 nM insulin
        sim.dr_d3 = model(time(1:idxT30), initcond_d, [param_d(1:end-1) ins_conc(5)], simOptions);    % 15 min simulation with 0.3 nM insulin
        sim.dr_d4 = model(time(1:idxT30), initcond_d, [param_d(1:end-1) ins_conc(6)], simOptions);    % 15 min simulation with 1 nM insulin
        sim.d_60  = model(time(1:idxT60), initcond_d, [param_d(1:end-1) ins_conc(7)], simOptions);     % 60 min simulation with 10 nM insulin
        sim.d_180 = model(time(1:idxT30), initcond_d, [param_d(1:end-1) ins_conc(8)], simOptions);   % 180 min simulation with 100 nM insulin %NOTE: used only for dosereponse at specific times, so 30 is fine
        sim.ds1_d = model(time(1:idxT4), initcond_d, [param_d(1:end-1) 1.2], simOptions);       % 4 min simulation with 1.2 nM insulin for double step data
        initcond2_d = sim.ds1_d.statevalues(end,:);                                 % using the statevalues from the ds1 simulation as initial conditions
        sim.ds2_d = model(time(1:idxT4), initcond2_d, [param_d(1:end-1) 10], simOptions);        % 4 min simulation with 10 nM insulin for double step data

    catch err
        cost = 1e24;
        return
    end
    
    idxT10=find(sim.n_60.time==10);
    idxT15=find(sim.n_60.time==15);
    idxT30=find(sim.n_60.time==30);

    %% Glucose stimulation simulations
    
    try
        % Controls
        param_n(index_gluc)=0.05; % Stimulating with 0.05 mM glucose
        sim.dr_gluc_n0 = model(time(1:idxT30), initcond_n, [param_n(1:end-1) 0], simOptions);     % initial conditions from baseline simulation, 30 min simulation without insulin
        initcond_n1a = sim.dr_n1a.statevalues(idxT15,:);                                   % initial conditions from 0.01 nM insulin simulation (dr_n1a)
        sim.dr_gluc_n1a = model(time(1:idxT30), initcond_n1a, [param_n(1:end-1) ins_conc(2)], simOptions);  % 30 min simulation with 0.01 nM insulin
        initcond_n1b = sim.dr_n1b.statevalues(idxT15,:);                                   % initial conditions from 0.03 nM insulin simulation (dr_n1b)
        sim.dr_gluc_n1b = model(time(1:idxT30), initcond_n1b, [param_n(1:end-1) ins_conc(3)], simOptions);  % 30 min simulation with 0.03 nM insulin 
        initcond_n2 = sim.dr_n2.statevalues(idxT15,:);                                     % initial conditions from 0.1 nM insulin simulation (dr_n2)
        sim.dr_gluc_n2 = model(time(1:idxT30), initcond_n2, [param_n(1:end-1) ins_conc(4)], simOptions);    % 30 min simulation with 0.1 nM insulin 
        initcond_n3 = sim.dr_n3.statevalues(idxT15,:);                                     % initial conditions from 0.3 nM insulin simulation (dr_n3)
        sim.dr_gluc_n3 = model(time(1:idxT30), initcond_n3, [param_n(1:end-1) ins_conc(5)], simOptions);    % 30 min simulation with 0.3 nM insulin 
        initcond_n4 = sim.dr_n4.statevalues(idxT15,:);                                     % initial conditions from .01 nM insulin simulation (dr_n4)
        sim.dr_gluc_n4 = model(time(1:idxT30), initcond_n4, [param_n(1:end-1) ins_conc(6)], simOptions);    % 30 min simulation with 1 nM insulin 
        initcond_n5 = sim.n_60.statevalues(idxT15,:);                                      % initial conditions from 10 nM insulin simulation (60_n)
        sim.dr_gluc_n5 = model(time(1:idxT30), initcond_n5, [param_n(1:end-1) ins_conc(7)], simOptions);    % 30 min simulation with 10 nM insulin 
        initcond_n6 = sim.n_180.statevalues(idxT15,:);                                     % initial conditions from 100 nM insulin simulation (180_n)
        sim.dr_gluc_n6 = model(time(1:idxT30), initcond_n6, [param_n(1:end-1) ins_conc(8)], simOptions);    % 30 min simulation with 100 nM insulin 
        
        % Diabetics
        param_d(index_gluc)=0.05; % Stimulating with 0.05 mM glucose
        sim.dr_gluc_d0 = model(time(1:idxT30), initcond_d, [param_d(1:end-1) 0], simOptions);     % initial conditions from baseline simulation, 30 min simulation without insulin
        initcond_d1a = sim.dr_d1a.statevalues(idxT15,:);                                   % initial conditions from 0.01 nM insulin simulation (dr_d1a)
        sim.dr_gluc_d1a = model(time(1:idxT30), initcond_d1a, [param_d(1:end-1) ins_conc(2)], simOptions);  % 30 min simulation with 0.01 nM insulin 
        initcond_d1b = sim.dr_d1b.statevalues(idxT15,:);                                   % initial conditions from 0.03 nM insulin simulation (dr_d1b)
        sim.dr_gluc_d1b = model(time(1:idxT30), initcond_d1b, [param_d(1:end-1) ins_conc(3)], simOptions);  % 30 min simulation with 0.03 nM insulin 
        initcond_d2 = sim.dr_d2.statevalues(idxT15,:);                                     % initial conditions from 0.1 nM insulin simulation (dr_d2)
        sim.dr_gluc_d2 = model(time(1:idxT30), initcond_d2, [param_d(1:end-1) ins_conc(4)], simOptions);    % 30 min simulation with 0.1 nM insulin 
        initcond_d3 = sim.dr_d3.statevalues(idxT15,:);                                     % initial conditions from 0.3 nM insulin simulation (dr_d3)
        sim.dr_gluc_d3 = model(time(1:idxT30), initcond_d3, [param_d(1:end-1) ins_conc(5)], simOptions);    % 30 min simulation with 0.3 nM insulin 
        initcond_d4 = sim.dr_d4.statevalues(idxT15,:);                                     % initial conditions from .01 nM insulin simulation (dr_d4)
        sim.dr_gluc_d4 = model(time(1:idxT30), initcond_d4, [param_d(1:end-1) ins_conc(6)], simOptions);    % 30 min simulation with 1 nM insulin 
        initcond_d5 = sim.d_60.statevalues(idxT15,:);                                      % initial conditions from 10 nM insulin simulation (60_d)
        sim.dr_gluc_d5 = model(time(1:idxT30), initcond_d5, [param_d(1:end-1) ins_conc(7)], simOptions);    % 30 min simulation with 10 nM insulin 
        initcond_d6 = sim.d_180.statevalues(idxT15,:);                                     % initial conditions from 100 nM insulin simulation (180_d)
        sim.dr_gluc_d6 = model(time(1:idxT30), initcond_d6, [param_d(1:end-1) ins_conc(8)], simOptions);    % 30 min simulation with 100 nM insulin 
    catch err
        cost = 1e24;
        return
    end
    
    % Collecting simulation data at 10 min for variable values to use for dose response cost calculations
    dose_n = [sim.n_0.variablevalues(end,:);sim.dr_n1a.variablevalues(idxT10,:); sim.dr_n2.variablevalues(idxT10,:); sim.dr_n3.variablevalues(idxT10,:); sim.dr_n4.variablevalues(idxT10,:); sim.n_60.variablevalues(idxT10,:); sim.n_180.variablevalues(idxT10,:)];
    dose_d = [sim.d_0.variablevalues(end,:);sim.dr_d1a.variablevalues(idxT10,:); sim.dr_d2.variablevalues(idxT10,:); sim.dr_d3.variablevalues(idxT10,:); sim.dr_d4.variablevalues(idxT10,:); sim.d_60.variablevalues(idxT10,:); sim.d_180.variablevalues(idxT10,:)];
    
    % Saving dose response data for IRp + IRip (variablevalue 1)
    %                               IRS1p (variablevalue 2)
    %                               IRS1307 (variablevalue 3)
    %                               PKB308 (variablevalue 4)
    %                               PKB473 (variablevalue 5)
    %                               AS160 (variablevalue 6)
    %                               S6 (variablevalue 9)
    % The doseresponse curves are normalized to 100 at maximum value
    % Baselines removed since data was treated that way
    % All dose response curves in this loop have insulin concentrations:
    % 1E-11,1E-10,3E-10,1E-9,1E-8,1E-7
    
    %Controls (standard)
    sim.dr_n.IR = dose_n(:,idxIR) - sim.n_0.variablevalues(end,idxIR); %baseline removal
    sim.dr_n.IR = sim.dr_n.IR*100/max(sim.dr_n.IR);
    sim.dr_n.IRS1 = dose_n(:,idxIRS1) - sim.n_0.variablevalues(end,idxIRS1); %baseline removal
    sim.dr_n.IRS1 = sim.dr_n.IRS1*100/max(sim.dr_n.IRS1);
    sim.dr_n.IRS1307 = dose_n(:,idxIRS1307) - sim.n_0.variablevalues(end,idxIRS1307); %baseline removal
    sim.dr_n.IRS1307 = sim.dr_n.IRS1307*100/max(sim.dr_n.IRS1307);
    sim.dr_n.PKB308 = dose_n(:,idxPKB308) - sim.n_0.variablevalues(end,idxPKB308); %baseline removal
    sim.dr_n.PKB308 = sim.dr_n.PKB308*100/max(sim.dr_n.PKB308);
    sim.dr_n.PKB473 = dose_n(:,idxPKB473) - sim.n_0.variablevalues(end,idxPKB473); %baseline removal
    sim.dr_n.PKB473 = sim.dr_n.PKB473*100/max(sim.dr_n.PKB473);
    sim.dr_n.AS160 = dose_n(:,idxAS160) - sim.n_0.variablevalues(end,idxAS160); %baseline removal
    sim.dr_n.AS160 = sim.dr_n.AS160*100/max(sim.dr_n.AS160);

    %Diabetes (standard)
    sim.dr_d.IR = dose_d(:,idxIR) - sim.d_0.variablevalues(end,idxIR); %baseline removal
    sim.dr_d.IR = sim.dr_d.IR*100/max(sim.dr_d.IR);
    sim.dr_d.IRS1 = dose_d(:,idxIRS1) - sim.d_0.variablevalues(end,idxIRS1); %baseline removal
    sim.dr_d.IRS1 = sim.dr_d.IRS1*100/max(sim.dr_d.IRS1);
    sim.dr_d.IRS1307 = dose_d(:,idxIRS1307) - sim.d_0.variablevalues(end,idxIRS1307); %baseline removal
    sim.dr_d.IRS1307 = sim.dr_d.IRS1307*100/max(sim.dr_d.IRS1307);
    sim.dr_d.PKB308 = dose_d(:,idxPKB308) - sim.d_0.variablevalues(end,idxPKB308); %baseline removal
    sim.dr_d.PKB308 = sim.dr_d.PKB308*100/max(sim.dr_d.PKB308);
    sim.dr_d.PKB473 = dose_d(:,idxPKB473) - sim.d_0.variablevalues(end,idxPKB473); %baseline removal
    sim.dr_d.PKB473 = sim.dr_d.PKB473*100/max(sim.dr_d.PKB473);
    sim.dr_d.AS160 = dose_d(:,idxAS160) - sim.d_0.variablevalues(end,idxAS160); %baseline removal
    sim.dr_d.AS160 = sim.dr_d.AS160*100/max(sim.dr_d.AS160);
    
    % Optional dose responses
    if any(idxS6)
        sim.dr_n.S6 = dose_n(:,idxS6) - sim.n_0.variablevalues(end,idxS6); %baseline removal
        sim.dr_n.S6 = sim.dr_n.S6*100/max(sim.dr_n.S6);
        sim.dr_d.S6 = dose_d(:,idxS6) - sim.d_0.variablevalues(end,idxS6); %baseline removal
        sim.dr_d.S6 = sim.dr_d.S6*100/max(sim.dr_d.S6);
    end

    if any(idxERK)
        sim.dr_n.ERK = dose_n(:,idxERK) - sim.n_0.variablevalues(end,idxERK); %baseline removal
        sim.dr_n.ERK = sim.dr_n.ERK*100/max(sim.dr_n.ERK);
        sim.dr_d.ERK = dose_d(:,idxERK) - sim.d_0.variablevalues(end,idxERK); %baseline removal
        sim.dr_d.ERK = sim.dr_d.ERK*100/max(sim.dr_d.ERK);
    end

    if any(idxElk1)
        dose_n30 = [sim.n_0.variablevalues(end,idxElk1);sim.dr_n1a.variablevalues(idxT30,idxElk1); sim.dr_n2.variablevalues(idxT30,idxElk1); sim.dr_n3.variablevalues(idxT30,idxElk1); sim.dr_n4.variablevalues(idxT30,idxElk1); sim.n_60.variablevalues(idxT30,idxElk1); sim.n_180.variablevalues(idxT30,idxElk1)];
        dose_d30 = [sim.d_0.variablevalues(end,idxElk1);sim.dr_d1a.variablevalues(idxT30,idxElk1); sim.dr_d2.variablevalues(idxT30,idxElk1); sim.dr_d3.variablevalues(idxT30,idxElk1); sim.dr_d4.variablevalues(idxT30,idxElk1); sim.d_60.variablevalues(idxT30,idxElk1); sim.d_180.variablevalues(idxT30,idxElk1)];
    
        sim.dr_n.Elk1 = dose_n30*100/max(dose_n30);
        sim.dr_d.Elk1 = dose_d30*100/max(dose_d30);
    end
    
    if any(idxFOXO)
        dose_n30 = [sim.n_0.variablevalues(end,idxFOXO);sim.dr_n001.variablevalues(idxT30,idxFOXO); sim.dr_n1a.variablevalues(idxT30,idxFOXO); sim.dr_n2.variablevalues(idxT30,idxFOXO); sim.dr_n3.variablevalues(idxT30,idxFOXO); sim.dr_n4.variablevalues(idxT30,idxFOXO); sim.n_60.variablevalues(idxT30,idxFOXO); sim.n_180.variablevalues(idxT30,idxFOXO)];
        dose_d30 = [sim.d_0.variablevalues(end,idxFOXO);sim.dr_d001.variablevalues(idxT30,idxFOXO); sim.dr_d1a.variablevalues(idxT30,idxFOXO); sim.dr_d2.variablevalues(idxT30,idxFOXO); sim.dr_d3.variablevalues(idxT30,idxFOXO); sim.dr_d4.variablevalues(idxT30,idxFOXO); sim.d_60.variablevalues(idxT30,idxFOXO); sim.d_180.variablevalues(idxT30,idxFOXO)];
           
        sim.dr_n.FOXO = dose_n30*100/max(dose_n30);
        sim.dr_d.FOXO = dose_d30*100/max(dose_d30);
        
        sim.dr_Foxo = [sim.n_0.variablevalues(end,idxFOXO);sim.dr_n001.variablevalues(idxT30,idxFOXO);sim.dr_n1a.variablevalues(idxT30,idxFOXO); sim.dr_n2.variablevalues(idxT30,idxFOXO); sim.dr_n3.variablevalues(idxT30,idxFOXO); sim.dr_n4.variablevalues(idxT30,idxFOXO); sim.n_60.variablevalues(idxT30,idxFOXO); sim.n_180.variablevalues(idxT30,idxFOXO)]./sim.n_180.variablevalues(idxT30,idxFOXO).*100;
        sim.dr_Foxo_d= [sim.d_0.variablevalues(end,idxFOXO);sim.dr_d001.variablevalues(idxT30,idxFOXO);sim.dr_d1a.variablevalues(idxT30,idxFOXO); sim.dr_d2.variablevalues(idxT30,idxFOXO); sim.dr_d3.variablevalues(idxT30,idxFOXO); sim.dr_d4.variablevalues(idxT30,idxFOXO); sim.d_60.variablevalues(idxT30,idxFOXO); sim.d_180.variablevalues(idxT30,idxFOXO)]./sim.d_180.variablevalues(idxT30,idxFOXO).*100;
    end


    %% Dose response with insulin and    GLUCOSE_cell (variablevalue 11)
    
    % The doseresponse curve is normalized to 100 at its maximum value
    % Insulin concentrations: 1E-11,3E-11,1E-10,3E-10,1E-9,1E-8,1E-7
    % Collecting simulation data for glucose after insulin and glucose stimulation
    GLUCOSE_dr_n1 = [sim.dr_gluc_n0.variablevalues(end,idxGLUCOSE) sim.dr_gluc_n1a.variablevalues(idxT30,idxGLUCOSE) sim.dr_gluc_n1b.variablevalues(idxT30,idxGLUCOSE) sim.dr_gluc_n2.variablevalues(idxT30,idxGLUCOSE) sim.dr_gluc_n3.variablevalues(idxT30,idxGLUCOSE) sim.dr_gluc_n4.variablevalues(idxT30,idxGLUCOSE) sim.dr_gluc_n5.variablevalues(idxT30,idxGLUCOSE) sim.dr_gluc_n6.variablevalues(idxT30,idxGLUCOSE)];
    GLUCOSE_dr_d1 = [sim.dr_gluc_d0.variablevalues(end,idxGLUCOSE) sim.dr_gluc_d1a.variablevalues(idxT30,idxGLUCOSE) sim.dr_gluc_d1b.variablevalues(idxT30,idxGLUCOSE) sim.dr_gluc_d2.variablevalues(idxT30,idxGLUCOSE) sim.dr_gluc_d3.variablevalues(idxT30,idxGLUCOSE) sim.dr_gluc_d4.variablevalues(idxT30,idxGLUCOSE) sim.dr_gluc_d5.variablevalues(idxT30,idxGLUCOSE) sim.dr_gluc_d6.variablevalues(idxT30,idxGLUCOSE)];
    
    GLUCOSE_dr_n = GLUCOSE_dr_n1 - sim.dr_gluc_n0.variablevalues(end,idxGLUCOSE);
    GLUCOSE_dr_d = GLUCOSE_dr_d1 - sim.dr_gluc_d0.variablevalues(end,idxGLUCOSE);
    GLUCOSE_dr_n = GLUCOSE_dr_n*100/max(GLUCOSE_dr_n);
    GLUCOSE_dr_d = GLUCOSE_dr_d*100/max(GLUCOSE_dr_d);
    
    sim.GLUCOSE_dr_n = GLUCOSE_dr_n;
    sim.GLUCOSE_dr_d = GLUCOSE_dr_d;

    %% EC50s
    EC50_n = zeros(7,1);
    EC50_d = zeros(7,1);
    ins_conc2 = [0.001 0.01 0.1 0.3 1 10 100];
    ins_conc3 = [0.001 0.003 0.01 0.1 0.3 1 10 100];
    
    try
        EC50_n(1) = 10.^(findX(log10(ins_conc2'),sim.dr_n.IR,50));
        EC50_d(1) = 10.^(findX(log10(ins_conc2'),sim.dr_d.IR,50));
        EC50_n(2) = 10.^(findX(log10(ins_conc2'),sim.dr_n.IRS1,50));
        EC50_d(2) = 10.^(findX(log10(ins_conc2'),sim.dr_d.IRS1,50));
        EC50_n(3) = 10.^(findX(log10(ins_conc2'),sim.dr_n.IRS1307,50));
        EC50_d(3) = 10.^(findX(log10(ins_conc2'),sim.dr_d.IRS1307,50));       
        EC50_n(4) = 10.^(findX(log10(ins_conc2'),sim.dr_n.PKB308,50));
        EC50_d(4) = 10.^(findX(log10(ins_conc2'),sim.dr_d.PKB308,50));        
        EC50_n(5) = 10.^(findX(log10(ins_conc2'),sim.dr_n.PKB473,50));
        EC50_d(5) = 10.^(findX(log10(ins_conc2'),sim.dr_d.PKB473,50));
        EC50_n(6) = 10.^(findX(log10(ins_conc2'),sim.dr_n.AS160,50));
        EC50_d(6) = 10.^(findX(log10(ins_conc2'),sim.dr_d.AS160,50));
        EC50_n(7) = 10.^(findX(log10(ins_conc3'),GLUCOSE_dr_n,50));
        EC50_d(7) = 10.^(findX(log10(ins_conc3'),GLUCOSE_dr_d,50));
    catch err
        cost = 1e24;
        return
    end
    sim.EC50_n=EC50_n;
    sim.EC50_d=EC50_d;
end