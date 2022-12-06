function [] = PlotAgreementGlu(optParam, model, expInd, expData, baseFolder)

if nargin<5
    baseFolder='';
end

colorNorm = [2 64 167]/256;
colorDiab = [0.6350 0.0780 0.1840];

modelName = char(model);

time=unique([expData.AllTimes 0:0.01:60]);

variables = IQMvariables(model);
idxIR=strcmp(variables,'measuredIR');
idxIRS1=strcmp(variables,'measuredIRS1');
idxIRS1307   = strcmp(variables,'measuredIRS1307');
idxPKB308    = strcmp(variables,'measuredPKB308');
idxPKB473    = strcmp(variables,'measuredPKB473');
idxAS160     = strcmp(variables,'measuredAS160');
idxmTORC1    = strcmp(variables,'measuredmTORC1');
idxS6K       = strcmp(variables,'measuredS6K');
idxS6        = strcmp(variables,'measuredS6');
idxGLUCOSE   = strcmp(variables,'measuredGLUCOSE');
idxIRi       = strcmp(variables,'fractionIRi');
idxmTORC2    = strcmp(variables,'measuredmTORC2');
idxERK       = strcmp(variables,'measuredERK');
idxElk1      = strcmp(variables,'measuredElk1');
idxFOXO      = strcmp(variables,'measuredFOXO');

% Combination: Multiply expData.time with expData.2p divided by expData.time at 10 min => expData.time(9)=1.
comb_expData.PKB473_time.response_n = expData.PKB473_time_renorm.response_n*expData.PKB473_2p.response_n(2); % normalising the 60 min data to the 10 min 2p data point
comb_expData.PKB473_time.sem_n = expData.PKB473_time_renorm.sem_n*expData.PKB473_2p.response_n(2);
comb_expData.PKB473_time.response_d = expData.PKB473_time_renorm.response_d*expData.PKB473_2p.response_d(2); % normalising the 60 min data to the 10 min 2p data point
comb_expData.PKB473_time.sem_d = expData.PKB473_time_renorm.sem_d*expData.PKB473_2p.response_d(2);

%% Loading all uncertainty parameter sets
baseFolder=strrep(baseFolder, 'Results','Results-PPL');
files=dir(sprintf('%s/**/*%s*.mat', baseFolder, modelName));
bestParam=optParam;
optParams=optParam;
for i = fliplr(1:length(files))
    load([files(i).folder '/' files(i).name],'optParam');
    optParams(i,:)=optParam;
end
optParams=[bestParam; unique(optParams,'rows')];

%%
fprintf('\nSimulating the uncertainty for the experiments of the glucose uptake submodel:\n')
for i = 1:size(optParams,1) 
    param=optParams(i,:); 
    param(expInd)=exp(param(expInd));

    [sim, tmpCost] = SimulateGlucose(param, model, time);
    if tmpCost > 0
        disp('Simulation has crashed.');
        return
    end
    idxT3=find(sim.n_60.time==3);
    idxT10=find(sim.n_60.time==10);
    idxT15=find(sim.n_60.time==15);
    idxT30=find(sim.n_60.time==30);

    %% Cost calculations
    % Extract dose response simulations

    if i ==1
        for f = fields(sim.dr_n)'
            dr_n.(f{:}) = repmat(sim.dr_n.(f{:}),1,3); 
            dr_d.(f{:}) = repmat(sim.dr_d.(f{:}),1,3);
        end
        GLUCOSE_dr_n = repmat(sim.GLUCOSE_dr_n',1,3);
        GLUCOSE_dr_d = repmat(sim.GLUCOSE_dr_d',1,3);        
        EC50_n = repmat(sim.EC50_n,1,3);
        EC50_d = repmat(sim.EC50_d,1,3);
    else
        for f = fields(sim.dr_n)'
            dr_n.(f{:})(:,2) = min(dr_n.(f{:})(:,2), sim.dr_n.(f{:})); 
            dr_d.(f{:})(:,2) = min(dr_d.(f{:})(:,2), sim.dr_d.(f{:}));
            dr_n.(f{:})(:,3) = max(dr_n.(f{:})(:,3), sim.dr_n.(f{:}));
            dr_d.(f{:})(:,3) = max(dr_d.(f{:})(:,3), sim.dr_d.(f{:}));
        end
        GLUCOSE_dr_n(:,2) = min(GLUCOSE_dr_n(:,2), sim.GLUCOSE_dr_n');
        GLUCOSE_dr_d(:,2) = min(GLUCOSE_dr_d(:,2), sim.GLUCOSE_dr_d');        
        GLUCOSE_dr_n(:,3) = max(GLUCOSE_dr_n(:,3), sim.GLUCOSE_dr_n');
        GLUCOSE_dr_d(:,3) = max(GLUCOSE_dr_d(:,3), sim.GLUCOSE_dr_d');    
        EC50_n(:,2) = min(EC50_n(:,2), sim.EC50_n);
        EC50_d(:,2) = min(EC50_d(:,2), sim.EC50_d);
        EC50_n(:,3) = max(EC50_n(:,3), sim.EC50_n);
        EC50_d(:,3) = max(EC50_d(:,3), sim.EC50_d);
    end
    %% Caluclate scaling factors
    % IRp + IRip

    % IR baseline/insulin stimulation 10 nM 10 min normal/diabetes
    IR_2p = [sim.n_0.variablevalues(end,idxIR);sim.n_60.variablevalues(idxT10,idxIR)];
    scaleIR_2p = lscov(IR_2p,expData.IR_2p.response_n);
    if i ==1
        IR_2p_n = repmat(scaleIR_2p*sim.n_60.variablevalues(:,idxIR),1,3);
        IR_2p_d = repmat(scaleIR_2p*sim.d_60.variablevalues(:,idxIR),1,3);
    else
        IR_2p_n(:,2)=min(IR_2p_n(:,2), scaleIR_2p*sim.n_60.variablevalues(:,idxIR));
        IR_2p_d(:,2)=min(IR_2p_d(:,2), scaleIR_2p*sim.d_60.variablevalues(:,idxIR));
        IR_2p_n(:,3)=max(IR_2p_n(:,3), scaleIR_2p*sim.n_60.variablevalues(:,idxIR));
        IR_2p_d(:,3)=max(IR_2p_d(:,3), scaleIR_2p*sim.d_60.variablevalues(:,idxIR));
    end

    % IRp 100nm insulin, 30 min
    % baseline subtracted
    IR_timep = sim.n_180.variablevalues(1:idxT30,idxIR)-sim.n_0.variablevalues(end,idxIR);
    tIdx=ismember(sim.n_180.time, expData.IR_time.time);
    scaleIR = lscov(IR_timep(tIdx),expData.IR_time.response);
    if i ==1
        IR_timep_n = repmat(scaleIR*(sim.n_180.variablevalues(:,idxIR)-sim.n_0.variablevalues(end,idxIR)),1,3);
        IR_timep_d = repmat(scaleIR*(sim.d_180.variablevalues(:,idxIR)-sim.d_0.variablevalues(end,idxIR)),1,3);
    else
        IR_timep_n(:,2)=min(IR_timep_n(:,2), scaleIR*(sim.n_180.variablevalues(:,idxIR)-sim.n_0.variablevalues(end,idxIR)));
        IR_timep_d(:,2)=min(IR_timep_d(:,2), scaleIR*(sim.d_180.variablevalues(:,idxIR)-sim.d_0.variablevalues(end,idxIR)));
        IR_timep_n(:,3)=max(IR_timep_n(:,3), scaleIR*(sim.n_180.variablevalues(:,idxIR)-sim.n_0.variablevalues(end,idxIR)));
        IR_timep_d(:,3)=max(IR_timep_d(:,3), scaleIR*(sim.d_180.variablevalues(:,idxIR)-sim.d_0.variablevalues(end,idxIR)));
    end

    % Internalized IR nM insulin, 10 min
    % percent of total IR
    if i ==1
        IRint_n = repmat(sim.n_60.variablevalues(:,idxIRi)*100,1,3);
        IRint_d = repmat(sim.d_60.variablevalues(:,idxIRi)*100,1,3);

    else
        IRint_n(:,2)=min(IRint_n(:,2), sim.n_60.variablevalues(:,idxIRi)*100);
        IRint_d(:,2)=min(IRint_d(:,2), sim.d_60.variablevalues(:,idxIRi)*100);
        IRint_n(:,3)=max(IRint_n(:,3), sim.n_60.variablevalues(:,idxIRi)*100);
        IRint_d(:,3)=max(IRint_d(:,3), sim.d_60.variablevalues(:,idxIRi)*100);
    end
    % IRS1
    % IRS1 baseline/insulin stimulation 10 nM 10 min normal/diabetes
    IRS1_2p = [sim.n_0.variablevalues(end,idxIRS1);sim.n_60.variablevalues(idxT10,idxIRS1)];
    scaleIRS1_2p = lscov(IRS1_2p,expData.IRS1_2p.response_n);
    if i ==1
        IRS1_2p_n = repmat(scaleIRS1_2p*sim.n_60.variablevalues(:,idxIRS1),1,3);
        IRS1_2p_d = repmat(scaleIRS1_2p*sim.d_60.variablevalues(:,idxIRS1),1,3);
    else
        IRS1_2p_n(:,2)=min(IRS1_2p_n(:,2), scaleIRS1_2p*sim.n_60.variablevalues(:,idxIRS1));
        IRS1_2p_d(:,2)=min(IRS1_2p_d(:,2), scaleIRS1_2p*sim.d_60.variablevalues(:,idxIRS1));
        IRS1_2p_n(:,3)=max(IRS1_2p_n(:,3), scaleIRS1_2p*sim.n_60.variablevalues(:,idxIRS1));
        IRS1_2p_d(:,3)=max(IRS1_2p_d(:,3), scaleIRS1_2p*sim.d_60.variablevalues(:,idxIRS1));
    end
    % IRS1 10nm insulin 3 min
    % baseline not subtracted
    % data normalised to max 100, simdata scaled with lscov
    tIdx=ismember(sim.n_60.time, expData.IRS1_time_3.time);
    IRS1_time3p = sim.n_60.variablevalues(1:idxT3,idxIRS1);
    scaleIRS1_10 = lscov(IRS1_time3p(tIdx),expData.IRS1_time_3.response);
    if i ==1
        IRS1_time3p_n = repmat(scaleIRS1_10*sim.n_60.variablevalues(:,idxIRS1),1,3);
        IRS1_time3p_d = repmat(scaleIRS1_10*sim.d_60.variablevalues(:,idxIRS1),1,3);
    else
        IRS1_time3p_n(:,2)=min(IRS1_time3p_n(:,2), scaleIRS1_10*sim.n_60.variablevalues(:,idxIRS1));
        IRS1_time3p_d(:,2)=min(IRS1_time3p_d(:,2), scaleIRS1_10*sim.d_60.variablevalues(:,idxIRS1));
        IRS1_time3p_n(:,3)=max(IRS1_time3p_n(:,3), scaleIRS1_10*sim.n_60.variablevalues(:,idxIRS1));
        IRS1_time3p_d(:,3)=max(IRS1_time3p_d(:,3), scaleIRS1_10*sim.d_60.variablevalues(:,idxIRS1));
    end
    % IRS1p 100nm insulin 30 min
    % baseline subtracted
    % data normalised to max 100, simdata scaled with lscov
    tIdx=ismember(sim.n_180.time, expData.IRS1_time_30.time);
    IRS1_time30p = sim.n_180.variablevalues(1:idxT30,idxIRS1)-sim.n_0.variablevalues(end,idxIRS1);
    scaleIRS1_30 = lscov(IRS1_time30p(tIdx),expData.IRS1_time_30.response);
    if i ==1
        IRS1_time30p_n = repmat(scaleIRS1_30*(sim.n_180.variablevalues(:,idxIRS1)-sim.n_0.variablevalues(end,idxIRS1)),1,3);
        IRS1_time30p_d = repmat(scaleIRS1_30*(sim.d_180.variablevalues(:,idxIRS1)-sim.d_0.variablevalues(end,idxIRS1)),1,3);
    else
        IRS1_time30p_n(:,2)=min(IRS1_time30p_n(:,2), scaleIRS1_30*(sim.n_180.variablevalues(:,idxIRS1)-sim.n_0.variablevalues(end,idxIRS1)));
        IRS1_time30p_d(:,2)=min(IRS1_time30p_d(:,2), scaleIRS1_30*(sim.d_180.variablevalues(:,idxIRS1)-sim.d_0.variablevalues(end,idxIRS1)));
        IRS1_time30p_n(:,3)=max(IRS1_time30p_n(:,3), scaleIRS1_30*(sim.n_180.variablevalues(:,idxIRS1)-sim.n_0.variablevalues(end,idxIRS1)));
        IRS1_time30p_d(:,3)=max(IRS1_time30p_d(:,3), scaleIRS1_30*(sim.d_180.variablevalues(:,idxIRS1)-sim.d_0.variablevalues(end,idxIRS1)));
    end
    % Cost for IRS1p double step (1.2 + 10 nm)
    % baseline not subtracted
    % simdata scaled with lscov
    tIdx_n=ismember(sim.ds1_n.time, expData.IRS1_ds.time);
    tIdx_d=ismember(sim.ds2_n.time(2:end)+4, expData.IRS1_ds.time);
    IRS1_ds = [sim.ds1_n.variablevalues(tIdx_n,idxIRS1); sim.ds2_n.variablevalues(tIdx_d,idxIRS1)];
    scaleIRS1ds = lscov(IRS1_ds,expData.IRS1_ds.response);
    if i ==1
        IRS1_dsp_n = repmat(scaleIRS1ds*[sim.ds1_n.variablevalues(:,idxIRS1); sim.ds2_n.variablevalues(:,idxIRS1)],1,3);
        IRS1_dsp_d = repmat(scaleIRS1ds*[sim.ds1_d.variablevalues(:,idxIRS1); sim.ds2_d.variablevalues(:,idxIRS1)],1,3);
    else
        IRS1_dsp_n(:,2)=min(IRS1_dsp_n(:,2), scaleIRS1ds*[sim.ds1_n.variablevalues(:,idxIRS1); sim.ds2_n.variablevalues(:,idxIRS1)]);
        IRS1_dsp_d(:,2)=min(IRS1_dsp_d(:,2), scaleIRS1ds*[sim.ds1_d.variablevalues(:,idxIRS1); sim.ds2_d.variablevalues(:,idxIRS1)]);
        IRS1_dsp_n(:,3)=max(IRS1_dsp_n(:,3), scaleIRS1ds*[sim.ds1_n.variablevalues(:,idxIRS1); sim.ds2_n.variablevalues(:,idxIRS1)]);
        IRS1_dsp_d(:,3)=max(IRS1_dsp_d(:,3), scaleIRS1ds*[sim.ds1_d.variablevalues(:,idxIRS1); sim.ds2_d.variablevalues(:,idxIRS1)]);
    end

    % IRS1307
    % Controls (n) and diabetes (d)
    % IRS1307 10nm insulin 60min
    % baseline not subtracted
    % simulation scaled with lscov
    tIdx=ismember(sim.n_60.time, expData.IRS1307_time.time);
    IRS1307_time = sim.n_60.variablevalues(tIdx,idxIRS1307);
    scaleIRS1307 = lscov(IRS1307_time,expData.IRS1307_time.response_n); %scale parameter for IRS1307
    if i ==1
        IRS1307_time_n = repmat(scaleIRS1307 * sim.n_60.variablevalues(:,idxIRS1307),1,3);
        IRS1307_time_d = repmat(scaleIRS1307 * sim.d_60.variablevalues(:,idxIRS1307),1,3);
    else
        IRS1307_time_n(:,2)=min(IRS1307_time_n(:,2), scaleIRS1307 * sim.n_60.variablevalues(:,idxIRS1307));
        IRS1307_time_d(:,2)=min(IRS1307_time_d(:,2), scaleIRS1307 * sim.d_60.variablevalues(:,idxIRS1307));
        IRS1307_time_n(:,3)=max(IRS1307_time_n(:,3), scaleIRS1307 * sim.n_60.variablevalues(:,idxIRS1307));
        IRS1307_time_d(:,3)=max(IRS1307_time_d(:,3), scaleIRS1307 * sim.d_60.variablevalues(:,idxIRS1307));
    end

    % PKB308
    tIdx=ismember(sim.n_60.time,expData.PKB308_time.time);
    PKB308_time = sim.n_60.variablevalues(tIdx,idxPKB308);
    scalePKB308 = lscov(PKB308_time,expData.PKB308_time.response_n);
    if i ==1
        PKB308_time_n = repmat(scalePKB308*sim.n_60.variablevalues(:,idxPKB308),1,3);
        PKB308_time_d = repmat(scalePKB308*sim.d_60.variablevalues(:,idxPKB308),1,3);
    else
        PKB308_time_n(:,2)=min(PKB308_time_n(:,2), scalePKB308*sim.n_60.variablevalues(:,idxPKB308));
        PKB308_time_d(:,2)=min(PKB308_time_d(:,2), scalePKB308*sim.d_60.variablevalues(:,idxPKB308));
        PKB308_time_n(:,3)=max(PKB308_time_n(:,3), scalePKB308*sim.n_60.variablevalues(:,idxPKB308));
        PKB308_time_d(:,3)=max(PKB308_time_d(:,3), scalePKB308*sim.d_60.variablevalues(:,idxPKB308));
    end
    % PKB473
    % 10nm insulin 60min normal (n) and diabetes (d)
    % baseline not subtracted, ExpData normalised to max 100. Normalised to 10 point in 2p ExpData to get relationship between n and d
    % simulation scaled with lscov

    tIdx = ismember(sim.n_60.time, expData.PKB473_time_renorm.time);
    PKB473_time = sim.n_60.variablevalues(tIdx,idxPKB473);
    scalePKB473 = lscov(PKB473_time,comb_expData.PKB473_time.response_n);
    if i ==1
        PKB473_time_n = repmat(scalePKB473*sim.n_60.variablevalues(:,idxPKB473),1,3);
        PKB473_time_d = repmat(scalePKB473*sim.d_60.variablevalues(:,idxPKB473),1,3);
    else
        PKB473_time_n(:,2)=min(PKB473_time_n(:,2), scalePKB473*sim.n_60.variablevalues(:,idxPKB473));
        PKB473_time_d(:,2)=min(PKB473_time_d(:,2), scalePKB473*sim.d_60.variablevalues(:,idxPKB473));
        PKB473_time_n(:,3)=max(PKB473_time_n(:,3), scalePKB473*sim.n_60.variablevalues(:,idxPKB473));
        PKB473_time_d(:,3)=max(PKB473_time_d(:,3), scalePKB473*sim.d_60.variablevalues(:,idxPKB473));
    end
    % AS160
    tIdx=ismember(sim.n_60.time, expData.AS160_time.time);
    AS160_time = sim.n_60.variablevalues(tIdx,idxAS160);
    scaleAS160 = lscov(AS160_time,expData.AS160_time.response_n);
    if i ==1
        AS160_time_n = repmat(scaleAS160*sim.n_60.variablevalues(:,idxAS160),1,3);
        AS160_time_d = repmat(scaleAS160*sim.d_60.variablevalues(:,idxAS160),1,3);
    else
        AS160_time_n(:,2)=min(AS160_time_n(:,2), scaleAS160*sim.n_60.variablevalues(:,idxAS160));
        AS160_time_d(:,2)=min(AS160_time_d(:,2), scaleAS160*sim.d_60.variablevalues(:,idxAS160));
        AS160_time_n(:,3)=max(AS160_time_n(:,3), scaleAS160*sim.n_60.variablevalues(:,idxAS160));
        AS160_time_d(:,3)=max(AS160_time_d(:,3), scaleAS160*sim.d_60.variablevalues(:,idxAS160));
    end
    % Glucose uptake two point data
    % 0 and 100 nM insulin 15 min
    % then 0.05 nM glucose 30 min
    GLUCOSE = [sim.dr_gluc_n0.variablevalues(end,idxGLUCOSE); sim.dr_gluc_n6.variablevalues(idxT30,idxGLUCOSE)];
    scaleGLUCOSE = lscov(GLUCOSE,expData.GLUCOSE_2p.response_n);
    if i ==1
        GLUCOSE_n = repmat(scaleGLUCOSE*[sim.dr_gluc_n0.variablevalues(end,idxGLUCOSE); sim.dr_gluc_n6.variablevalues(idxT30,idxGLUCOSE)],1,3);
        GLUCOSE_d = repmat(scaleGLUCOSE*[sim.dr_gluc_d0.variablevalues(end,idxGLUCOSE); sim.dr_gluc_d6.variablevalues(idxT30,idxGLUCOSE)],1,3);
    else
        GLUCOSE_n(:,2)=min(GLUCOSE_n(:,2), scaleGLUCOSE*[sim.dr_gluc_n0.variablevalues(end,idxGLUCOSE); sim.dr_gluc_n6.variablevalues(idxT30,idxGLUCOSE)]);
        GLUCOSE_d(:,2)=min(GLUCOSE_d(:,2), scaleGLUCOSE*[sim.dr_gluc_d0.variablevalues(end,idxGLUCOSE); sim.dr_gluc_d6.variablevalues(idxT30,idxGLUCOSE)]);
        GLUCOSE_n(:,3)=max(GLUCOSE_n(:,3), scaleGLUCOSE*[sim.dr_gluc_n0.variablevalues(end,idxGLUCOSE); sim.dr_gluc_n6.variablevalues(idxT30,idxGLUCOSE)]);
        GLUCOSE_d(:,3)=max(GLUCOSE_d(:,3), scaleGLUCOSE*[sim.dr_gluc_d0.variablevalues(end,idxGLUCOSE); sim.dr_gluc_d6.variablevalues(idxT30,idxGLUCOSE)]);
    end

    % S6K, 10nM, 60min
    % Fold over basal for normal and diabetes data respectively
    % Simdata scaled with lscov for n and d respectively
    if any(idxS6K)
        tIdx=ismember(sim.n_60.time, expData.S6K_time.time);
        S6K_time = sim.n_60.variablevalues(tIdx,idxS6K);
        scaleS6K = lscov(S6K_time, expData.S6K_time.response_n);
        if i ==1
            S6K_time_n = repmat(scaleS6K*sim.n_60.variablevalues(:,idxS6K),1,3);
            S6K_time_d = repmat(scaleS6K*sim.d_60.variablevalues(:,idxS6K),1,3);
        else
            S6K_time_n(:,2)=min(S6K_time_n(:,2), scaleS6K*sim.n_60.variablevalues(:,idxS6K));
            S6K_time_d(:,2)=min(S6K_time_d(:,2), scaleS6K*sim.d_60.variablevalues(:,idxS6K));
            S6K_time_n(:,3)=max(S6K_time_n(:,3), scaleS6K*sim.n_60.variablevalues(:,idxS6K));
            S6K_time_d(:,3)=max(S6K_time_d(:,3), scaleS6K*sim.d_60.variablevalues(:,idxS6K));
        end
    end
    % S6, 10 nM insulin, 60 min
    % Fold over basal for normal and diabetes data respectively
    % Simdata scaled with lscov for n and d respectively
    if any(idxS6)
        tIdx=ismember(sim.n_60.time, expData.S6_time.time);
        S6_time = sim.n_60.variablevalues(tIdx,idxS6);
        scaleS6 = lscov(S6_time, expData.S6_time.response_n);
        if i ==1
            S6_time_n = repmat(scaleS6*sim.n_60.variablevalues(:,idxS6),1,3);
            S6_time_d = repmat(scaleS6*sim.d_60.variablevalues(:,idxS6),1,3);
        else
            S6_time_n(:,2)=min(S6_time_n(:,2), scaleS6*sim.n_60.variablevalues(:,idxS6));
            S6_time_d(:,2)=min(S6_time_d(:,2), scaleS6*sim.d_60.variablevalues(:,idxS6));
            S6_time_n(:,3)=max(S6_time_n(:,3), scaleS6*sim.n_60.variablevalues(:,idxS6));
            S6_time_d(:,3)=max(S6_time_d(:,3), scaleS6*sim.d_60.variablevalues(:,idxS6));
        end
    end

    % ERK
    % 10nm insulin 60min normal (n) and diabetes (d)
    % baseline not subtracted, simulation scaled with lscov
    if any(idxERK)
        tIdx=ismember(sim.n_60.time, expData.ERK_time.time);
        ERK_time = sim.n_60.variablevalues(tIdx,idxERK);
        scaleERK = lscov(ERK_time, expData.ERK_time.mean_n);
        if i ==1
            ERK_time_n = repmat(scaleERK*sim.n_60.variablevalues(:,idxERK),1,3);
            ERK_time_d = repmat(scaleERK*sim.d_60.variablevalues(:,idxERK),1,3);
        else
            ERK_time_n(:,2)=min(ERK_time_n(:,2), scaleERK*sim.n_60.variablevalues(:,idxERK));
            ERK_time_d(:,2)=min(ERK_time_d(:,2), scaleERK*sim.d_60.variablevalues(:,idxERK));
            ERK_time_n(:,3)=max(ERK_time_n(:,3), scaleERK*sim.n_60.variablevalues(:,idxERK));
            ERK_time_d(:,3)=max(ERK_time_d(:,3), scaleERK*sim.d_60.variablevalues(:,idxERK));
        end
    end
    % Elk1
    % 10nm insulin 60min normal (n) and diabetes (d)
    % baseline not subtracted, simulation scaled as % of max
    if any(idxElk1)
        tIdx=ismember(sim.n_60.time, expData.Elk1_time.time);
        Elk1_time = sim.n_60.variablevalues(tIdx,idxElk1);
        scaleElk1 = lscov(Elk1_time, expData.Elk1_time.response_n');
        if i ==1
            Elk1_time_n = repmat(scaleElk1*sim.n_60.variablevalues(:,idxElk1),1,3);
            Elk1_time_d = repmat(scaleElk1*sim.d_60.variablevalues(:,idxElk1),1,3);
        else
            Elk1_time_n(:,2)=min(Elk1_time_n(:,2), scaleElk1*sim.n_60.variablevalues(:,idxElk1));
            Elk1_time_d(:,2)=min(Elk1_time_d(:,2), scaleElk1*sim.d_60.variablevalues(:,idxElk1));
            Elk1_time_n(:,3)=max(Elk1_time_n(:,3), scaleElk1*sim.n_60.variablevalues(:,idxElk1));
            Elk1_time_d(:,3)=max(Elk1_time_d(:,3), scaleElk1*sim.d_60.variablevalues(:,idxElk1));
        end
    end
    % FOXO
    % 10nm insulin 60min normal (n) and diabetes (d)
    % baseline not subtracted, simulation scaled as % of max
    if any(idxFOXO)
        tIdx=ismember(sim.n_60.time, expData.Foxo_time.time);
        FOXO_time = sim.n_60.variablevalues(tIdx,idxFOXO);
        scaleFOXO = lscov(FOXO_time, expData.Foxo_time.response_n);
        if i ==1
            FOXO_time_n = repmat(scaleFOXO*sim.n_60.variablevalues(:,idxFOXO),1,3);
            FOXO_time_d = repmat(scaleFOXO*sim.d_60.variablevalues(:,idxFOXO),1,3);
        else
            FOXO_time_n(:,2)=min(FOXO_time_n(:,2), scaleFOXO*sim.n_60.variablevalues(:,idxFOXO));
            FOXO_time_d(:,2)=min(FOXO_time_d(:,2), scaleFOXO*sim.d_60.variablevalues(:,idxFOXO));
            FOXO_time_n(:,3)=max(FOXO_time_n(:,3), scaleFOXO*sim.n_60.variablevalues(:,idxFOXO));
            FOXO_time_d(:,3)=max(FOXO_time_d(:,3), scaleFOXO*sim.d_60.variablevalues(:,idxFOXO));
        end
    end

    if i ==1
        mTORC1_n = repmat(sim.n_60.variablevalues(:,idxmTORC1),1,3);
        mTORC1_d = repmat(sim.d_60.variablevalues(:,idxmTORC1),1,3);
    else
        mTORC1_n(:,2)=min(mTORC1_n(:,2), sim.n_60.variablevalues(:,idxmTORC1));
        mTORC1_d(:,2)=min(mTORC1_d(:,2), sim.d_60.variablevalues(:,idxmTORC1));
        mTORC1_n(:,3)=max(mTORC1_n(:,3), sim.n_60.variablevalues(:,idxmTORC1));
        mTORC1_d(:,3)=max(mTORC1_d(:,3), sim.d_60.variablevalues(:,idxmTORC1));
    end

    if i ==1
        mTORC2_n = repmat(sim.n_60.variablevalues(:,idxmTORC2),1,3);
        mTORC2_d = repmat(sim.d_60.variablevalues(:,idxmTORC2),1,3);
    else
        mTORC2_n(:,2)=min(mTORC2_n(:,2), sim.n_60.variablevalues(:,idxmTORC2));
        mTORC2_d(:,2)=min(mTORC2_d(:,2), sim.d_60.variablevalues(:,idxmTORC2));
        mTORC2_n(:,3)=max(mTORC2_n(:,3), sim.n_60.variablevalues(:,idxmTORC2));
        mTORC2_d(:,3)=max(mTORC2_d(:,3), sim.d_60.variablevalues(:,idxmTORC2));
    end
    if i==1
        fprintf('%i of %i \n|',i,size(optParams,1))
    elseif mod(i,50)==0
        fprintf(' %i of %i \n',i,size(optParams,1))
    else
        fprintf('|')
    end
end
fprintf('\n')

%% Plot figures
figsize = [450 75 900 600];

conc = 1e-9*[0.001 0.01 0.1 0.3 1 10 100];

set(figure(61),'Position', figsize);
figure(61)

m = 6;
n = 5;
pos = reshape(1:m*n,n,m)'; % Vertical mode

pIdx = 1;

%Dose response data for IRp + IRip (variablevalue 1)
%                               IRS1p (variablevalue 2)
%                               IRS1307 (variablevalue 3)
%                               PKB308 (variablevalue 4)
%                               PKB473 (variablevalue 5)
%                               AS160 (variablevalue 6)
%                               S6 (variablevalue 9)
% The doseresponse curves are normalized to 100 at maximum value
% Baselines removed since data was treated that way
% All dose response curves in this loop have insulin concentrations:
% 1E-11,1E-10,3E-10,1E-9,1E-8,1E-7 M (model is in nM)
% Dose response

subplot(m, n, pos(pIdx))
hold on
hn  = fill([conc flip(conc)]', [dr_n.IR(:,2); flip(dr_n.IR(:,3))], colorNorm); set(hn,'facealpha',0.5,'edgealpha',0);
hd  = fill([conc flip(conc)]', [dr_d.IR(:,2); flip(dr_d.IR(:,3))], colorDiab); set(hd,'facealpha',0.5,'edgealpha',0);
plot(conc, dr_n.IR(:,1),'-', 'color', colorNorm, 'LineWidth',2)
plot(conc, dr_d.IR(:,1),'-', 'color', colorDiab, 'LineWidth',2)
errorbar(conc,expData.IR_dr.response_n,expData.IR_dr.sem_n,'o', 'color', colorNorm,'MarkerFaceColor','auto','MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
errorbar(conc,expData.IR_dr.response_d,expData.IR_dr.sem_d,'o', 'color', colorDiab,'MarkerFaceColor','auto','MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
axis('tight')
set(gca, 'XScale', 'log', 'YTick',[0 50 100],'XTick',[1e-11 1e-9 1e-7]);
box off
xlabel('[Insulin] (M)')
ylabel('% of max')
title('IR-YP')
set(gca,'FontSize',11)
pIdx = pIdx+1;

subplot(m, n, pos(pIdx))
hold on
hn  = fill([conc flip(conc)]', [dr_n.IRS1(:,2); flip(dr_n.IRS1(:,3))], colorNorm); set(hn,'facealpha',0.5,'edgealpha',0);
hd  = fill([conc flip(conc)]', [dr_d.IRS1(:,2); flip(dr_d.IRS1(:,3))], colorDiab); set(hd,'facealpha',0.5,'edgealpha',0);
plot(conc, dr_n.IRS1(:,1),'-', 'color', colorNorm, 'LineWidth',2)
plot(conc, dr_d.IRS1(:,1),'-', 'color', colorDiab, 'LineWidth',2)
errorbar(conc,expData.IRS1_dr.response_n,expData.IRS1_dr.sem_n,'o', 'color', colorNorm,'MarkerFaceColor','auto','MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
errorbar(conc,expData.IRS1_dr.response_d,expData.IRS1_dr.sem_d,'o', 'color', colorDiab,'MarkerFaceColor','auto','MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
axis('tight')
set(gca, 'XScale', 'log','YTick',[0 50 100],'XTick',[1e-11 1e-9 1e-7]);
box off
xlabel('[Insulin] (M)')
ylabel('% of max')
title('IRS1-YP')
set(gca,'FontSize',11)
pIdx = pIdx+1;

subplot(m, n, pos(pIdx))
hold on
hn  = fill([conc flip(conc)]', [dr_n.IRS1307(:,2); flip(dr_n.IRS1307(:,3))], colorNorm); set(hn,'facealpha',0.5,'edgealpha',0);
hd  = fill([conc flip(conc)]', [dr_d.IRS1307(:,2); flip(dr_d.IRS1307(:,3))], colorDiab); set(hd,'facealpha',0.5,'edgealpha',0);
plot(conc, dr_n.IRS1307(:,1),'-', 'color', colorNorm, 'LineWidth',2)
plot(conc, dr_d.IRS1307(:,1),'-', 'color', colorDiab, 'LineWidth',2)
errorbar(conc,expData.IRS1307_dr.response_n,expData.IRS1307_dr.sem_n,'o', 'color', colorNorm,'MarkerFaceColor','auto','MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
errorbar(conc,expData.IRS1307_dr.response_d,expData.IRS1307_dr.sem_d,'o', 'color', colorDiab,'MarkerFaceColor','auto','MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
axis('tight')
set(gca, 'XScale', 'log','YTick',[0 50 100],'XTick',[1e-11 1e-9 1e-7]);
box off
xlabel('[Insulin] (M)')
ylabel('% of max')
title('IRS1-S307P')
set(gca,'FontSize',11)
pIdx = pIdx+1;

subplot(m, n, pos(pIdx))
hold on
hn  = fill([conc flip(conc)]', [dr_n.PKB308(:,2); flip(dr_n.PKB308(:,3))], colorNorm); set(hn,'facealpha',0.5,'edgealpha',0);
hd  = fill([conc flip(conc)]', [dr_d.PKB308(:,2); flip(dr_d.PKB308(:,3))], colorDiab); set(hd,'facealpha',0.5,'edgealpha',0);
plot(conc, dr_n.PKB308(:,1), '-', 'color', colorNorm, 'LineWidth',2)
plot(conc, dr_d.PKB308(:,1), '-', 'color', colorDiab, 'LineWidth',2)
errorbar(conc, expData.PKB308_dr.response_n,expData.PKB308_dr.sem_n,'o', 'color', colorNorm,'MarkerFaceColor','auto','MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
errorbar(conc, expData.PKB308_dr.response_d,expData.PKB308_dr.sem_d,'o', 'color', colorDiab,'MarkerFaceColor','auto','MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
axis('tight')
set(gca, 'XScale', 'log','YTick',[0 50 100],'XTick',[1e-11 1e-9 1e-7]);
box off
xlabel('[Insulin] (M)')
ylabel('% of max')
title('PKB-S308P')
set(gca,'FontSize',11)
pIdx = pIdx+1;

subplot(m, n, pos(pIdx))
hold on
hn  = fill([conc flip(conc)]', [dr_n.PKB473(:,2); flip(dr_n.PKB473(:,3))], colorNorm); set(hn,'facealpha',0.5,'edgealpha',0);
hd  = fill([conc flip(conc)]', [dr_d.PKB473(:,2); flip(dr_d.PKB473(:,3))], colorDiab); set(hd,'facealpha',0.5,'edgealpha',0);
plot(conc, dr_n.PKB473(:,1), '-', 'color', colorNorm, 'LineWidth',2)
plot(conc, dr_d.PKB473(:,1), '-', 'color', colorDiab, 'LineWidth',2)
errorbar(conc, expData.PKB473_dr.response_n,expData.PKB473_dr.sem_n,'o', 'color', colorNorm,'MarkerFaceColor','auto','MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
errorbar(conc, expData.PKB473_dr.response_d,expData.PKB473_dr.sem_d,'o', 'color', colorDiab,'MarkerFaceColor','auto','MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
axis('tight')
set(gca, 'XScale', 'log','YTick',[0 50 100],'XTick',[1e-11 1e-9 1e-7]);
box off
xlabel('[Insulin] (M)')
ylabel('% of max')
title('PKB-T473P')
set(gca,'FontSize',11)
pIdx = pIdx+1;

subplot(m, n, pos(pIdx))
hold on
hn  = fill([conc flip(conc)]', [dr_n.AS160(:,2); flip(dr_n.AS160(:,3))], colorNorm); set(hn,'facealpha',0.5,'edgealpha',0);
hd  = fill([conc flip(conc)]', [dr_d.AS160(:,2); flip(dr_d.AS160(:,3))], colorDiab); set(hd,'facealpha',0.5,'edgealpha',0);
plot(conc, dr_n.AS160(:,1), '-', 'color', colorNorm, 'LineWidth',2)
plot(conc, dr_d.AS160(:,1), '-', 'color', colorDiab, 'LineWidth',2)
errorbar(conc, expData.AS160_dr.response_n,expData.AS160_dr.sem_n,'o', 'color', colorNorm,'MarkerFaceColor','auto','MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
errorbar(conc, expData.AS160_dr.response_d,expData.AS160_dr.sem_d,'o', 'color', colorDiab,'MarkerFaceColor','auto','MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
axis('tight')
set(gca, 'XScale', 'log','YTick',[0 50 100],'XTick',[1e-11 1e-9 1e-7]);
box off
xlabel('[Insulin] (M)')
ylabel('% of max')
title('AS160-T642P')
set(gca,'FontSize',11)
pIdx = pIdx+1;

subplot(m, n, pos(pIdx))
hold on
conc_glucose = expData.GLUCOSE_dr.conc;
hn  = fill([conc_glucose; flip(conc_glucose)], [GLUCOSE_dr_n(:,2); flip(GLUCOSE_dr_n(:,3))], colorNorm); set(hn,'facealpha',0.5,'edgealpha',0);
hd  = fill([conc_glucose; flip(conc_glucose)], [GLUCOSE_dr_d(:,2); flip(GLUCOSE_dr_d(:,3))], colorDiab); set(hd,'facealpha',0.5,'edgealpha',0);
plot(conc_glucose,GLUCOSE_dr_n(:,1),'-', 'color', colorNorm, 'LineWidth',2)
plot(conc_glucose,GLUCOSE_dr_d(:,1),'-', 'color', colorDiab, 'LineWidth',2)
errorbar(expData.GLUCOSE_dr.conc,expData.GLUCOSE_dr.response_n,expData.GLUCOSE_dr.sem_n,'o', 'color', colorNorm,'MarkerFaceColor','auto','MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
errorbar(expData.GLUCOSE_dr.conc,expData.GLUCOSE_dr.response_d,expData.GLUCOSE_dr.sem_d,'o', 'color', colorDiab,'MarkerFaceColor','auto','MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
axis('tight')
set(gca, 'XScale', 'log','YTick',[0 50 100],'XTick',[1e-11 1e-9 1e-7]);
box off
xlabel('[Insulin] (M)')
ylabel('% of max')
title({'Glucose uptake'}) %, 0.05 nM ins.
set(gca,'FontSize',11)
pIdx = pIdx+1;

if any(idxS6)
    subplot(m, n, pos(pIdx))
    hold on
    hn  = fill([conc flip(conc)]', [dr_n.S6(:,2); flip(dr_n.S6(:,3))], colorNorm); set(hn,'facealpha',0.5,'edgealpha',0);
    hd  = fill([conc flip(conc)]', [dr_d.S6(:,2); flip(dr_d.S6(:,3))], colorDiab); set(hd,'facealpha',0.5,'edgealpha',0);
    plot(conc, dr_n.S6(:,1)', '-', 'color', colorNorm, 'LineWidth',2)
    plot(conc, dr_d.S6(:,1)', '-', 'color', colorDiab, 'LineWidth',2)
    errorbar(conc, expData.S6_dr.response_n,expData.S6_dr.sem_n,'o', 'color', colorNorm,'MarkerFaceColor','auto','MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
    axis('tight')
    set(gca, 'XScale', 'log','YTick',[0 50 100],'XTick',[1e-11 1e-9 1e-7]);
    box off
    ylabel('% of max')
    xlabel('[Insulin] (M)')
    title('S6-S235/S236P')
    set(gca,'FontSize',11)
    pIdx = pIdx+1;
end
if any(idxERK)
    subplot(m, n, pos(pIdx))
    hold on
    hn  = fill([conc flip(conc)]', [dr_n.ERK(:,2); flip(dr_n.ERK(:,3))], colorNorm); set(hn,'facealpha',0.5,'edgealpha',0);
    hd  = fill([conc flip(conc)]', [dr_d.ERK(:,2); flip(dr_d.ERK(:,3))], colorDiab); set(hd,'facealpha',0.5,'edgealpha',0);
    plot(conc, dr_n.ERK(:,1)', '-', 'color', colorNorm, 'LineWidth', 2)
    plot(conc, dr_d.ERK(:,1)', '-', 'color', colorDiab, 'LineWidth', 2)
    errorbar(conc, expData.ERK_dr.response_n,expData.ERK_dr.sem_n,'o', 'color', colorNorm,'MarkerFaceColor', colorNorm,'MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
    errorbar(conc, expData.ERK_dr.response_d,expData.ERK_dr.sem_d,'o', 'color', colorDiab,'MarkerFaceColor', colorDiab,'MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
    axis('tight')
    set(gca, 'XScale', 'log','YTick',[0 50 100],'XTick',[1e-11 1e-9 1e-7]);
    box off
    xlabel('[Insulin], M')
    ylabel('% of max')
    title('ERK-T202/Y204P')
    set(gca,'FontSize',11)
    pIdx = pIdx+1;
end

if any(idxFOXO)
    subplot(m, n, pos(pIdx))
    hold on
    concFoxo = [1e-4*1e-9 conc]; % adding the 0.001 ins concentration from the foxo sim.
    hn  = fill([concFoxo flip(concFoxo)]', [dr_n.FOXO(:,2); flip(dr_n.FOXO(:,3))], colorNorm); set(hn,'facealpha',0.5,'edgealpha',0);
    hd  = fill([concFoxo flip(concFoxo)]', [dr_d.FOXO(:,2); flip(dr_d.FOXO(:,3))], colorDiab); set(hd,'facealpha',0.5,'edgealpha',0);
    plot(concFoxo, dr_n.FOXO(:,1), '-', 'color', colorNorm, 'LineWidth', 2)
    plot(concFoxo, dr_d.FOXO(:,1), '-', 'color', colorDiab, 'LineWidth', 2)
    errorbar(concFoxo, expData.Foxo_raw_dr.response_n,expData.Foxo_raw_dr.sem_n,'o', 'color', colorNorm,'MarkerFaceColor', colorNorm,'MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
    errorbar(concFoxo, expData.Foxo_raw_dr.response_d,expData.Foxo_raw_dr.sem_d,'o', 'color', colorDiab,'MarkerFaceColor', colorDiab,'MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
    axis('tight')
    set(gca, 'XScale', 'log','YTick',[0 50 100],'XTick',[1e-11 1e-9 1e-7]);
    box off
    xlabel('[Insulin], M')
    ylabel('% of max')
    title('FOXO-S256P')
    set(gca,'FontSize',11)
    pIdx = pIdx+1;
end

% Time curves
subplot(m, n, pos(pIdx))
hold on
t = sim.n_180.time(1:idxT30)';
hn  = fill([t; flip(t)], [IR_timep_n(1:idxT30,2); flip(IR_timep_n(1:idxT30,3))], colorNorm); set(hn,'facealpha',0.5,'edgealpha',0);
hd  = fill([t; flip(t)], [IR_timep_d(1:idxT30,2); flip(IR_timep_d(1:idxT30,3))], colorDiab); set(hd,'facealpha',0.5,'edgealpha',0);

plot(t,IR_timep_n(1:idxT30,1),'-', 'color', colorNorm, 'LineWidth',2)
plot(t,IR_timep_d(1:idxT30,1),'-', 'color', colorDiab, 'LineWidth',2)

errorbar(expData.IR_time.time,expData.IR_time.response,expData.IR_time.sem,'o', 'color', colorNorm,'MarkerFaceColor','auto','MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
axis('tight')
xlabel('Time (min)')
ylabel('a.u.')
title('IR-YP')
set(gca,'FontSize',11)
pIdx = pIdx+1;

subplot(m, n, pos(pIdx))
hold on
t = sim.n_60.time(1:idxT10)';
hn  = fill([t; flip(t)], [IR_2p_n(1:idxT10,2); flip(IR_2p_n(1:idxT10,3))], colorNorm); set(hn,'facealpha',0.5,'edgealpha',0);
hd  = fill([t; flip(t)], [IR_2p_d(1:idxT10,2); flip(IR_2p_d(1:idxT10,3))], colorDiab); set(hd,'facealpha',0.5,'edgealpha',0);

plot(t,IR_2p_n(1:idxT10,1),'-', 'color', colorNorm, 'LineWidth',2);
plot(t,IR_2p_d(1:idxT10,1),'-', 'color', colorDiab, 'LineWidth',2);

errorbar([0 10],expData.IR_2p.response_n,expData.IR_2p.sem_n,'o', 'color', colorNorm,'MarkerFaceColor','auto','MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
errorbar([0 10],expData.IR_2p.response_d,expData.IR_2p.sem_d,'o', 'color', colorDiab,'MarkerFaceColor','auto','MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
axis('tight')
xlabel('Time (min)')
ylabel('a.u.')
title('IR-YP')
set(gca,'FontSize',11)
pIdx = pIdx+1;

subplot(m, n, pos(pIdx))
hold on
t = sim.n_60.time(1:idxT15)';
hn  = fill([t; flip(t)], [IRint_n(1:idxT15,2); flip(IRint_n(1:idxT15,3))], colorNorm); set(hn,'facealpha',0.5,'edgealpha',0);
hd  = fill([t; flip(t)], [IRint_d(1:idxT15,2); flip(IRint_d(1:idxT15,3))], colorDiab); set(hd,'facealpha',0.5,'edgealpha',0);
plot(t,IRint_n(1:idxT15,1),'-', 'color', colorNorm, 'LineWidth',2)
plot(t,IRint_d(1:idxT15,1),'-', 'color', colorDiab, 'LineWidth',2)
errorbar(expData.IR_int.time,expData.IR_int.response, expData.IR_int.sem,'o', 'color', colorNorm,'MarkerFaceColor','auto','MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
xlabel('Time (min)')
ylabel('a.u.')
title('% internalized IR, 10 nM ins') %.
set(gca,'FontSize',11)
pIdx = pIdx+1;

subplot(m, n, pos(pIdx))
hold on
t = sim.n_60.time(1:idxT3)';
hn  = fill([t; flip(t)], [IRS1_time3p_n(1:idxT3,2); flip(IRS1_time3p_n(1:idxT3,3))], colorNorm); set(hn,'facealpha',0.5,'edgealpha',0);
hd  = fill([t; flip(t)], [IRS1_time3p_d(1:idxT3,2); flip(IRS1_time3p_d(1:idxT3,3))], colorDiab); set(hd,'facealpha',0.5,'edgealpha',0);
plot(t,IRS1_time3p_n(1:idxT3,1),'-', 'color', colorNorm, 'LineWidth',2)
plot(t,IRS1_time3p_d(1:idxT3,1),'-', 'color', colorDiab, 'LineWidth',2)
errorbar(expData.IRS1_time_3.time,expData.IRS1_time_3.response,expData.IRS1_time_3.sem,'o', 'color', colorNorm,'MarkerFaceColor','auto','MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
axis('tight')
xlabel('Time (min)')
ylabel('a.u.')
title('IRS1-YP, 10 nM ins.')
set(gca,'FontSize',11)
pIdx = pIdx+1;

subplot(m, n, pos(pIdx))
hold on
t = sim.n_180.time(1:idxT30)';
hn  = fill([t; flip(t)], [IRS1_time30p_n(1:idxT30,2); flip(IRS1_time30p_n(1:idxT30,3))], colorNorm); set(hn,'facealpha',0.5,'edgealpha',0);
hd  = fill([t; flip(t)], [IRS1_time30p_d(1:idxT30,2); flip(IRS1_time30p_d(1:idxT30,3))], colorDiab); set(hd,'facealpha',0.5,'edgealpha',0);
plot(t,IRS1_time30p_n(1:idxT30,1),'-', 'color', colorNorm, 'LineWidth',2)
plot(t,IRS1_time30p_d(1:idxT30,1),'-', 'color', colorDiab, 'LineWidth',2)
errorbar(expData.IRS1_time_30.time,expData.IRS1_time_30.response,expData.IRS1_time_30.sem,'o', 'color', colorNorm,'MarkerFaceColor','auto','MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
axis('tight')
xlabel('Time (min)')
ylabel('a.u.')
title('IRS1-YP, 100 nM ins.')
set(gca,'FontSize',11)
pIdx = pIdx+1;

subplot(m, n, pos(pIdx))
hold on
t = [sim.ds1_n.time 4+sim.ds2_n.time]';
hn  = fill([t; flip(t)], [IRS1_dsp_n(:,2); flip(IRS1_dsp_n(:,3))], colorNorm); set(hn,'facealpha',0.5,'edgealpha',0);
hd  = fill([t; flip(t)], [IRS1_dsp_d(:,2); flip(IRS1_dsp_d(:,3))], colorDiab); set(hd,'facealpha',0.5,'edgealpha',0);
plot(t,IRS1_dsp_n(:,1),'-', 'color', colorNorm, 'LineWidth',2)
plot(t,IRS1_dsp_d(:,1),'-', 'color', colorDiab, 'LineWidth',2)
errorbar(expData.IRS1_ds.time,expData.IRS1_ds.response,expData.IRS1_ds.sem,'o', 'color', colorNorm,'MarkerFaceColor','auto','MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
axis('tight')
xlabel('Time (min)')
ylabel('a.u.')
title('IRS1-YP, 1.2 + 10 nM ins.')
set(gca,'FontSize',11)
pIdx = pIdx+1;

subplot(m, n, pos(pIdx))
hold on
t = sim.n_60.time(1:idxT10)';
hn  = fill([t; flip(t)], [IRS1_2p_n(1:idxT10,2); flip(IRS1_2p_n(1:idxT10,3))], colorNorm); set(hn,'facealpha',0.5,'edgealpha',0);
hd  = fill([t; flip(t)], [IRS1_2p_d(1:idxT10,2); flip(IRS1_2p_d(1:idxT10,3))], colorDiab); set(hd,'facealpha',0.5,'edgealpha',0);
plot(t,IRS1_2p_n(1:idxT10,1),'-', 'color', colorNorm, 'LineWidth',2);
plot(t,IRS1_2p_d(1:idxT10,1),'-', 'color', colorDiab, 'LineWidth',2);
errorbar([0 10],expData.IRS1_2p.response_n,expData.IRS1_2p.sem_n,'o', 'color', colorNorm,'MarkerFaceColor','auto','MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
errorbar([0 10],expData.IRS1_2p.response_d,expData.IRS1_2p.sem_d,'o', 'color', colorDiab,'MarkerFaceColor','auto','MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
axis('tight')
xlabel('Time (min)')
ylabel('a.u.')
title('IRS1-YP, 10 nM ins.')
set(gca,'FontSize',11)
pIdx = pIdx+1;

subplot(m, n, pos(pIdx))
hold on
t = sim.n_60.time';
hn  = fill([t; flip(t)], [IRS1307_time_n(:,2); flip(IRS1307_time_n(:,3))], colorNorm); set(hn,'facealpha',0.5,'edgealpha',0);
hd  = fill([t; flip(t)], [IRS1307_time_d(:,2); flip(IRS1307_time_d(:,3))], colorDiab); set(hd,'facealpha',0.5,'edgealpha',0);
plot(t,IRS1307_time_n(:,1),'-', 'color', colorNorm, 'LineWidth',2)
plot(t,IRS1307_time_d(:,1),'-', 'color', colorDiab, 'LineWidth',2)
errorbar(expData.IRS1307_time.time,expData.IRS1307_time.response_n,expData.IRS1307_time.sem_n,'o', 'color', colorNorm,'MarkerFaceColor','auto','MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
errorbar(expData.IRS1307_time.time,expData.IRS1307_time.response_d,expData.IRS1307_time.sem_d,'o', 'color', colorDiab,'MarkerFaceColor','auto','MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
axis([0 60 0 5])
xlabel('Time (min)')
ylabel('a.u.')
title('IRS1-S307P, 10 nM ins.')
set(gca,'FontSize',11)
pIdx = pIdx+1;

subplot(m, n, pos(pIdx))
hold on
t = sim.n_60.time';
hn  = fill([t; flip(t)], [PKB308_time_n(:,2); flip(PKB308_time_n(:,3))], colorNorm); set(hn,'facealpha',0.5,'edgealpha',0);
hd  = fill([t; flip(t)], [PKB308_time_d(:,2); flip(PKB308_time_d(:,3))], colorDiab); set(hd,'facealpha',0.5,'edgealpha',0);
plot(t,PKB308_time_n(:,1),'-', 'color', colorNorm, 'LineWidth',2)
plot(t,PKB308_time_d(:,1),'-', 'color', colorDiab, 'LineWidth',2)
errorbar(expData.PKB308_time.time,expData.PKB308_time.response_n,expData.PKB308_time.sem_n,'o', 'color', colorNorm,'MarkerFaceColor','auto','MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
errorbar(expData.PKB308_time.time,expData.PKB308_time.response_d,expData.PKB308_time.sem_d,'o', 'color', colorDiab,'MarkerFaceColor','auto','MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
axis('tight')
title('PKB-S308P, 10 nM ins.')
xlabel('Time (min)')
ylabel('a.u.')
set(gca,'FontSize',11)
pIdx = pIdx+1;

subplot(m, n, pos(pIdx))
hold on
t = sim.n_60.time';
hn  = fill([t; flip(t)], [PKB473_time_n(:,2); flip(PKB473_time_n(:,3))], colorNorm); set(hn,'facealpha',0.5,'edgealpha',0);
hd  = fill([t; flip(t)], [PKB473_time_d(:,2); flip(PKB473_time_d(:,3))], colorDiab); set(hd,'facealpha',0.5,'edgealpha',0);
plot(t,PKB473_time_n(:,1),'-', 'color', colorNorm, 'LineWidth',2)
plot(t,PKB473_time_d(:,1),'-', 'color', colorDiab, 'LineWidth',2)
errorbar(expData.PKB473_time.time,comb_expData.PKB473_time.response_n,comb_expData.PKB473_time.sem_n,'o', 'color', colorNorm,'MarkerFaceColor','auto','MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
errorbar(expData.PKB473_time.time,comb_expData.PKB473_time.response_d,comb_expData.PKB473_time.sem_d,'o', 'color', colorDiab,'MarkerFaceColor','auto','MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
axis('tight')
title('PKB-T473P, 10 nM ins.')
xlabel('Time (min)')
ylabel('a.u.')
set(gca,'FontSize',11)
pIdx = pIdx+1;


subplot(m, n, pos(pIdx))
hold on
t = sim.n_60.time';
hn  = fill([t; flip(t)], [mTORC1_n(:,2); flip(mTORC1_n(:,3))], colorNorm); set(hn,'facealpha',0.5,'edgealpha',0);
hd  = fill([t; flip(t)], [mTORC1_d(:,2); flip(mTORC1_d(:,3))], colorDiab); set(hd,'facealpha',0.5,'edgealpha',0);
plot(t,mTORC1_n(:,1),'-', 'color', colorNorm, 'LineWidth', 2)
plot(t,mTORC1_d(:,1),'-', 'color', colorDiab, 'LineWidth', 2)
axis('tight')
box off
xlabel('Time (min)')
ylabel('a.u.')
title('Active mTORC1')
set(gca,'YTick',[0 50 100],'XTick',[0 10 20 30]);
set(gca,'FontSize',11)
pIdx = pIdx+1;

subplot(m, n, pos(pIdx))
hold on
t = sim.n_60.time';
hn  = fill([t; flip(t)], [AS160_time_n(:,2); flip(AS160_time_n(:,3))], colorNorm); set(hn,'facealpha',0.5,'edgealpha',0);
hd  = fill([t; flip(t)], [AS160_time_d(:,2); flip(AS160_time_d(:,3))], colorDiab); set(hd,'facealpha',0.5,'edgealpha',0);
plot(t,AS160_time_n(:,1),'-', 'color', colorNorm, 'LineWidth',2)
plot(t,AS160_time_d(:,1),'-', 'color', colorDiab, 'LineWidth',2)
errorbar(expData.AS160_time.time,expData.AS160_time.response_n,expData.AS160_time.sem_n,'o', 'color', colorNorm,'MarkerFaceColor','auto','MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
errorbar(expData.AS160_time.time,expData.AS160_time.response_d,expData.AS160_time.sem_d,'o', 'color', colorDiab,'MarkerFaceColor','auto','MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
% axis([0 60 0 2.6])
axis('tight')
title('AS160-T642P, 10 nM ins.')
xlabel('Time (min)')
ylabel('a.u.')
set(gca,'FontSize',11)
pIdx = pIdx+1;

subplot(m, n, pos(pIdx))
bar([1 4 2 5],  [expData.GLUCOSE_2p.response_n' GLUCOSE_n(:,1)'], 'facecolor', colorNorm)
hold on
bar([7 10 8 11],[expData.GLUCOSE_2p.response_d' GLUCOSE_d(:,1)'], 'facecolor', colorDiab)
xticks([1.5 4.5 7.5 10.5])
xticklabels({'N/0 nM','N/100 nM','D/0 nM', 'D/100 nM'})
set(gca,'YTick',[]);
set(gca,'FontSize',11)
title('Glucose uptake')
text(10,1900,'Data','FontSize',8,'FontWeight','Bold','Rotation',60)
text(11,1900,'Sim.','FontSize',8,'FontWeight','Bold','Rotation',60)
box off
set(gca,'FontSize',11)
pIdx = pIdx+1;


subplot(m, n, pos(pIdx))
hold on
t = sim.n_60.time';
hn  = fill([t; flip(t)], [mTORC2_n(:,2); flip(mTORC2_n(:,3))], colorNorm); set(hn,'facealpha',0.5,'edgealpha',0);
hd  = fill([t; flip(t)], [mTORC2_d(:,2); flip(mTORC2_d(:,3))], colorDiab); set(hd,'facealpha',0.5,'edgealpha',0);
plot(t,mTORC2_n(:,1),'-', 'color', colorNorm, 'LineWidth', 2)
plot(t,mTORC2_d(:,1),'-', 'color', colorDiab, 'LineWidth', 2)
axis('tight')
box off
xlabel('Time (min)')
ylabel('a.u.')
title('Active mTORC2')
set(gca,'FontSize',11)
pIdx = pIdx+1;

if any(idxS6K)
    subplot(m, n, pos(pIdx))
    hold on
    t = sim.n_60.time';
    hn  = fill([t; flip(t)], [S6K_time_n(:,2); flip(S6K_time_n(:,3))], colorNorm); set(hn,'facealpha',0.5,'edgealpha',0);
    hd  = fill([t; flip(t)], [S6K_time_d(:,2); flip(S6K_time_d(:,3))], colorDiab); set(hd,'facealpha',0.5,'edgealpha',0);
    plot(t,S6K_time_n(:,1),'-', 'color', colorNorm, 'LineWidth',2)
    plot(t,S6K_time_d(:,1),'-', 'color', colorDiab, 'LineWidth',2)
    errorbar(expData.S6K_time.time,expData.S6K_time.response_n,expData.S6K_time.sem_n,'o', 'color', colorNorm,'MarkerFaceColor','auto','MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
    errorbar(expData.S6K_time.time,expData.S6K_time.response_d,expData.S6K_time.sem_d,'o', 'color', colorDiab,'MarkerFaceColor','auto','MarkerSize',4', 'LineWidth',1.5, 'capsize',8)
    axis('tight')
    xlabel('Time (min)')
    ylabel('a.u.')
    title('S6K-T389P, 10 nM ins.')
    set(gca,'FontSize',11)
    pIdx = pIdx+1;

end
if any(idxS6)
    subplot(m, n, pos(pIdx))
    hold on
    t = sim.n_60.time';
    hn  = fill([t; flip(t)], [S6_time_n(:,2); flip(S6_time_n(:,3))], colorNorm); set(hn,'facealpha',0.5,'edgealpha',0);
    hd  = fill([t; flip(t)], [S6_time_d(:,2); flip(S6_time_d(:,3))], colorDiab); set(hd,'facealpha',0.5,'edgealpha',0);
    plot(t,S6_time_n(:,1),'-', 'color', colorNorm, 'LineWidth',2)
    plot(t,S6_time_d(:,1),'-', 'color', colorDiab, 'LineWidth',2)
    errorbar(expData.S6_time.time,expData.S6_time.response_n,expData.S6_time.sem_n,'o', 'color', colorNorm,'MarkerFaceColor','auto','MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
    errorbar(expData.S6_time.time,expData.S6_time.response_d,expData.S6_time.sem_d,'o', 'color', colorDiab,'MarkerFaceColor','auto','MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
    axis('tight')
    xlabel('Time (min)')
    ylabel('a.u.')
    title('S6-S235/S236P, 10 nM ins.')
    set(gca,'FontSize',11)
    pIdx = pIdx+1;
end

if any(idxERK)
    subplot(m, n, pos(pIdx))
    hold on
    t = sim.n_60.time';
    hn  = fill([t; flip(t)], [ERK_time_n(:,2); flip(ERK_time_n(:,3))], colorNorm); set(hn,'facealpha',0.5,'edgealpha',0);
    hd  = fill([t; flip(t)], [ERK_time_d(:,2); flip(ERK_time_d(:,3))], colorDiab); set(hd,'facealpha',0.5,'edgealpha',0);
    plot(t,ERK_time_n(:,1),'-', 'color', colorNorm, 'LineWidth', 2)
    plot(t,ERK_time_d(:,1),'-', 'color', colorDiab, 'LineWidth', 2)
    errorbar(expData.ERK_time.time,expData.ERK_time.mean_n,expData.ERK_time.sem_n,'o', 'color', colorNorm,'MarkerFaceColor', colorNorm,'MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
    errorbar(expData.ERK_time.time,expData.ERK_time.mean_d,expData.ERK_time.sem_d,'o', 'color', colorDiab,'MarkerFaceColor', colorDiab,'MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
    axis('tight')
    set(gca,'YTick',[0 2 4 6],'XTick',[0 20 40 60]);
    box off
    xlabel('Time (min)')
    ylabel('a.u.')
    title('ERK-T202/Y204P')
    set(gca,'FontSize',11)
    pIdx = pIdx+1;
end

if any(idxElk1)
    subplot(m, n, pos(pIdx))
    hold on
    t = sim.n_60.time';
    hn  = fill([t; flip(t)], [Elk1_time_n(:,2); flip(Elk1_time_n(:,3))], colorNorm); set(hn,'facealpha',0.5,'edgealpha',0);
    hd  = fill([t; flip(t)], [Elk1_time_d(:,2); flip(Elk1_time_d(:,3))], colorDiab); set(hd,'facealpha',0.5,'edgealpha',0);
    plot(t,Elk1_time_n(:,1),'-', 'color', colorNorm, 'LineWidth', 2)
    plot(t,Elk1_time_d(:,1),'-', 'color', colorDiab, 'LineWidth', 2)
    errorbar(expData.Elk1_time.time,expData.Elk1_time.response_n,expData.Elk1_time.sem_n,'o', 'color', colorNorm,'MarkerFaceColor', colorNorm,'MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
    errorbar(expData.Elk1_time.time,expData.Elk1_time.response_d,expData.Elk1_time.sem_d,'o', 'color', colorDiab,'MarkerFaceColor', colorDiab,'MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
    axis('tight')
    set(gca,'YTick',[0 50 100],'XTick',[0 20 40 60]);
    box off
    xlabel('Time (min)')
    ylabel('%a.u.')
    title('Elk1-S383P')
    set(gca,'FontSize',11)
    pIdx = pIdx+1;
end

if any(idxFOXO)
    subplot(m, n, pos(pIdx))
    hold on
    t = sim.n_60.time';
    hn  = fill([t; flip(t)], [FOXO_time_n(:,2); flip(FOXO_time_n(:,3))], colorNorm); set(hn,'facealpha',0.5,'edgealpha',0);
    hd  = fill([t; flip(t)], [FOXO_time_d(:,2); flip(FOXO_time_d(:,3))], colorDiab); set(hd,'facealpha',0.5,'edgealpha',0);
    plot(t,FOXO_time_n(:,1),'-', 'color', colorNorm, 'LineWidth', 2)
    plot(t,FOXO_time_d(:,1),'-', 'color', colorDiab, 'LineWidth', 2)
    errorbar(expData.Foxo_time.time,expData.Foxo_time.response_n,expData.Foxo_time.sem_n,'o', 'color', colorNorm,'MarkerFaceColor', colorNorm,'MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
    errorbar(expData.Foxo_time.time,expData.Foxo_time.response_d,expData.Foxo_time.sem_d,'o', 'color', colorDiab,'MarkerFaceColor', colorDiab,'MarkerSize',4, 'LineWidth',1.5, 'capsize',8)
    axis('tight')
    set(gca,'YTick',[0 1 2 3],'XTick',[0 20 40 60]);
    box off
    xlabel('Time (min)')
    ylabel(' a.u.')
    title('FOXO-S256P')
    set(gca,'FontSize',11)
    pIdx = pIdx+1;
end
set(figure(61), 'outerposition',[249 1 1727 1408], 'PaperType','a4') % Vertical mode

end

