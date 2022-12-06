function [costTotal, sims] = CostNyman(param, model, expInd, expData, limit)
if ~isempty(expInd)
    param(expInd)=exp(param(expInd)); % To adapt based on if running in the combined model or in the individual submodels
end

t=table();
sims.Normal=t;
sims.Diabetes=t;

%% Indices

variables = IQMvariables(model);
idxIR=strcmp(variables,'measuredIR');
idxIRS1=strcmp(variables,'measuredIRS1');
idxIRS1307   = strcmp(variables,'measuredIRS1307');
idxPKB308    = strcmp(variables,'measuredPKB308');
idxPKB473    = strcmp(variables,'measuredPKB473');
idxAS160     = strcmp(variables,'measuredAS160');
idxS6K       = strcmp(variables,'measuredS6K');
idxS6        = strcmp(variables,'measuredS6');
idxGLUT4     = strcmp(variables,'measuredGLUT4');
idxGLUCOSE   = strcmp(variables,'measuredGLUCOSE');
idxIRi       = strcmp(variables,'fractionIRi');
idxERK       = strcmp(variables,'measuredERK');
idxElk1      = strcmp(variables,'measuredElk1');
idxFOXO      = strcmp(variables,'measuredFOXO');

time = expData.AllTimes;

[sim, tmpCost] = SimulateGlucose(param, model, time);
if tmpCost > 0
    costTotal=tmpCost;
    return
end
idxT10=find(sim.n_60.time==10);
idxT30=find(sim.n_60.time==30);

%% Cost calculations
costs=zeros(1,24);

% Extract dose response simulations
dr_n = sim.dr_n;
dr_d = sim.dr_d;

costs(1) = sum((expData.IR_dr.response_n-dr_n.IR).^2./expData.IR_dr.sem_n.^2);
costs(1) = costs(1) + sum((expData.IR_dr.response_d-dr_d.IR).^2./expData.IR_dr.sem_d.^2); %Note: t1 has small SEM
sims.Normal.IR_dr=[dr_n.IR'; expData.IR_dr.conc'];
sims.Diabetes.IR_dr=[dr_d.IR'; expData.IR_dr.conc'];
costs(2) = sum((expData.IRS1_dr.response_n-dr_n.IRS1).^2./expData.IRS1_dr.sem_n.^2); %Note: t1 has sigma=0
costs(2) = costs(2) + sum((expData.IRS1_dr.response_d-dr_d.IRS1).^2./expData.IRS1_dr.sem_d.^2);%Note: t1 has sigma=0
sims.Normal.IRS1_dr=[dr_n.IRS1'; expData.IRS1_dr.conc'];
sims.Diabetes.IRS1_dr=[dr_d.IRS1'; expData.IRS1_dr.conc'];
costs(3) = sum((expData.IRS1307_dr.response_n-dr_n.IRS1307).^2./expData.IRS1307_dr.sem_n.^2); %Note: t1 has sigma=0
costs(3) = costs(3) + sum((expData.IRS1307_dr.response_d-dr_d.IRS1307).^2./expData.IRS1307_dr.sem_d.^2); %Note: t1 has sigma=0
sims.Normal.IRS1307_dr=[dr_n.IRS1307'; expData.IRS1307_dr.conc'];
sims.Diabetes.IRS1307_dr=[dr_d.IRS1307'; expData.IRS1307_dr.conc'];
costs(4) = sum(((expData.PKB308_dr.response_n-dr_n.PKB308).^2)./expData.PKB308_dr.sem_n.^2);
costs(4) = costs(4) + sum(((expData.PKB308_dr.response_d-dr_d.PKB308).^2)./expData.PKB308_dr.sem_d.^2);
sims.Normal.PKB308_dr=[dr_n.PKB308'; expData.PKB308_dr.conc'];
sims.Diabetes.PKB308_dr=[dr_d.PKB308'; expData.PKB308_dr.conc'];
costs(5) = sum(((expData.PKB473_dr.response_n'-dr_n.PKB473).^2)./expData.PKB473_dr.sem_n'.^2);
costs(5) = costs(5) + sum(((expData.PKB473_dr.response_d'-dr_d.PKB473).^2)./expData.PKB473_dr.sem_d'.^2);
sims.Normal.PKB473_dr=[dr_n.PKB473'; expData.PKB473_dr.conc];
sims.Diabetes.PKB473_dr=[dr_d.PKB473'; expData.PKB473_dr.conc];
costs(6) = sum(((expData.AS160_dr.response_n-dr_n.AS160).^2)./expData.AS160_dr.sem_n.^2);
costs(6) = costs(6) + sum(((expData.AS160_dr.response_d-dr_d.AS160).^2)./expData.AS160_dr.sem_d.^2);
sims.Normal.AS160_dr=[dr_n.AS160'; expData.AS160_dr.conc'];
sims.Diabetes.AS160_dr=[dr_d.AS160'; expData.AS160_dr.conc'];

if any(idxS6)
    costs(7) = sum(((expData.S6_dr.response_n-dr_n.S6).^2)./expData.S6_dr.sem_n.^2);
    sims.Normal.S6_dr=[dr_n.S6'; expData.S6_dr.conc'];
    sims.Diabetes.S6_dr=[dr_d.S6'; expData.S6_dr.conc'];
end

if any(idxERK)
    costs(22) = sum((expData.ERK_dr.response_n-dr_n.ERK).^2./expData.ERK_dr.sem_n.^2);
    costs(22) = costs(22) + sum((expData.ERK_dr.response_d-dr_d.ERK).^2./expData.ERK_dr.sem_d.^2); %Note: t1 has small SEM
    sims.Normal.ERK_dr=[dr_n.ERK'; expData.ERK_dr.conc'];
    sims.Diabetes.ERK_dr=[dr_d.ERK'; expData.ERK_dr.conc'];
end

if any(idxFOXO)
    costs(23) = sum((expData.Foxo_raw_dr.response_n-dr_n.FOXO).^2./expData.Foxo_raw_dr.sem_n.^2);
    costs(23) = costs(23) + sum((expData.Foxo_raw_dr.response_d-dr_d.FOXO).^2./expData.Foxo_raw_dr.sem_d.^2); %Note: t1 has small SEM
    sims.Normal.Foxo_raw_dr=[dr_n.FOXO'; expData.Foxo_raw_dr.insulin'];
    sims.Diabetes.Foxo_raw_dr=[dr_d.FOXO'; expData.Foxo_raw_dr.insulin'];
end

costs(8) = sum(((expData.GLUCOSE_dr.response_n - sim.GLUCOSE_dr_n').^2)./expData.GLUCOSE_dr.sem_n.^2);
costs(8) = costs(8) + sum(((expData.GLUCOSE_dr.response_d - sim.GLUCOSE_dr_d').^2)./expData.GLUCOSE_dr.sem_d.^2);
sims.Normal.GLUCOSE_dr=[sim.GLUCOSE_dr_n; expData.GLUCOSE_dr.conc'];
sims.Diabetes.GLUCOSE_dr=[sim.GLUCOSE_dr_d; expData.GLUCOSE_dr.conc'];

costs(9) =sum((1e9*expData.EC50.n(1:5) - sim.EC50_n(1:5)).^2);
costs(9) = costs(9) + sum((1e9*expData.EC50.d(1:5) - sim.EC50_d(1:5)).^2);
sims.Normal.EC50=[sim.EC50_d(1:5)'; 1:5];
sims.Diabetes.EC50=[sim.EC50_d(1:5)'; 1:5];

%% Cost for time curves

% IRp + IRip
% IR baseline/insulin stimulation 10 nM 10 min normal/diabetes
IR_2p_n = [sim.n_0.variablevalues(end,idxIR);sim.n_60.variablevalues(idxT10,idxIR)];
IR_2p_d = [sim.d_0.variablevalues(end,idxIR);sim.d_60.variablevalues(idxT10,idxIR)];
scaleIR_2p = lscov(IR_2p_n,expData.IR_2p.response_n, expData.IR_2p.sem_n.^-2);
costs(10) = sum((expData.IR_2p.response_n-scaleIR_2p*IR_2p_n).^2./expData.IR_2p.sem_n.^2);
costs(10) = costs(10) + sum((expData.IR_2p.response_d-scaleIR_2p*IR_2p_d).^2./expData.IR_2p.sem_d.^2);
sims.Normal.IR_2p=[scaleIR_2p*IR_2p_n'; [0 10]];
sims.Diabetes.IR_2p=[scaleIR_2p*IR_2p_d'; [0 10]];

% IRp 100nm insulin, 30 min
% baseline subtracted
% normalised max to 100
tIdx=ismember(sim.n_180.time, expData.IR_time.time);
IR_time = sim.n_180.variablevalues(tIdx,idxIR)-sim.n_0.variablevalues(end,idxIR);
scaleIR = lscov(IR_time,expData.IR_time.response, expData.IR_time.sem.^-2);
costs(11) = sum((expData.IR_time.response-scaleIR*IR_time).^2./expData.IR_time.sem.^2);
sims.Normal.IR_time=[scaleIR*IR_time'; expData.IR_time.time'];
sims.Diabetes.IR_time=[scaleIR*IR_time'; expData.IR_time.time'];

% Internalized IR nM insulin, 10 min
% percent of total IR
IRint = sim.n_60.variablevalues(idxT10,idxIRi)*100;
costs(11) = costs(11) + ((expData.IR_int.response-IRint)/expData.IR_int.sem)^2;
sims.Normal.IRint=[IRint'; expData.IR_int.time'];

% IRS1
% IRS1 baseline/insulin stimulation 10 nM 10 min normal/diabetes
IRS1_2p_n = [sim.n_0.variablevalues(end,idxIRS1);sim.n_60.variablevalues(idxT10,idxIRS1)];
IRS1_2p_d = [sim.d_0.variablevalues(end,idxIRS1);sim.d_60.variablevalues(idxT10,idxIRS1)];
scaleIRS1_2p = lscov(IRS1_2p_n,expData.IRS1_2p.response_n, expData.IRS1_2p.sem_n.^-2);
costs(12) = sum((expData.IRS1_2p.response_n-scaleIRS1_2p*IRS1_2p_n).^2./expData.IRS1_2p.sem_n.^2);
costs(12) = costs(12) + sum((expData.IRS1_2p.response_d-scaleIRS1_2p*IRS1_2p_d).^2./expData.IRS1_2p.sem_d.^2);
sims.Normal.IRS1_2p=[scaleIRS1_2p*IRS1_2p_n'; [0 10]];
sims.Diabetes.IRS1_2p=[scaleIRS1_2p*IRS1_2p_d'; [0 10]];

% IRS1 10nm insulin 3 min
% baseline not subtracted
% data normalised to max 100, simdata scaled with lscov
tIdx=ismember(sim.n_60.time, expData.IRS1_time_3.time);
IRS1_time3 = sim.n_60.variablevalues(tIdx,idxIRS1);
scaleIRS1_10 = lscov(IRS1_time3,expData.IRS1_time_3.response, expData.IRS1_time_3.sem.^-2);
costs(13) = sum((expData.IRS1_time_3.response-scaleIRS1_10*IRS1_time3).^2./expData.IRS1_time_3.sem.^2);
sims.Normal.IRS1_time3=[scaleIRS1_10*IRS1_time3'; expData.IRS1_time_3.time'];

% IRS1p 100nm insulin 30 min
% baseline subtracted
% data normalised to max 100, simdata scaled with lscov
tIdx=ismember(sim.n_180.time, expData.IRS1_time_30.time);
IRS1_time30 = sim.n_180.variablevalues(tIdx,idxIRS1)-sim.n_0.variablevalues(end,idxIRS1);
scaleIRS1_30 = lscov(IRS1_time30,expData.IRS1_time_30.response, expData.IRS1_time_30.sem.^-2);
costs(13) = costs(13) + sum((expData.IRS1_time_30.response-scaleIRS1_30*IRS1_time30).^2./expData.IRS1_time_30.sem.^2);
sims.Normal.IRS1_time30=[scaleIRS1_30*IRS1_time30'; expData.IRS1_time_30.time'];

% Cost for IRS1p double step (1.2 + 10 nm)
% baseline not subtracted
% simdata scaled with lscov
tIdx_n=ismember(sim.ds1_n.time, expData.IRS1_ds.time);
tIdx_d=ismember(sim.ds2_n.time(2:end)+4, expData.IRS1_ds.time);
IRS1_ds = [sim.ds1_n.variablevalues(tIdx_n,idxIRS1); sim.ds2_n.variablevalues(tIdx_d,idxIRS1)];
scaleIRS1ds = lscov(IRS1_ds,expData.IRS1_ds.response, expData.IRS1_ds.sem.^-2);
costs(13) = costs(13) + sum((expData.IRS1_ds.response-scaleIRS1ds*IRS1_ds).^2./expData.IRS1_ds.sem.^2);
sims.Normal.IRS1_ds=[scaleIRS1ds*IRS1_ds'; expData.IRS1_ds.time'];

% IRS1307
% Controls (n) and diabetes (d)
% IRS1307 10nm insulin 60min
% baseline not subtracted
% simulation scaled with lscov
tIdx=ismember(sim.n_60.time, expData.IRS1307_time.time);
IRS1307_time_n = sim.n_60.variablevalues(tIdx,idxIRS1307);
IRS1307_time_d = sim.d_60.variablevalues(tIdx,idxIRS1307);
scaleIRS1307 = lscov(IRS1307_time_n,expData.IRS1307_time.response_n, expData.IRS1307_time.sem_n.^-2); %scale parameter for IRS1307
costs(14) = sum(((expData.IRS1307_time.response_n - scaleIRS1307*IRS1307_time_n).^2)./expData.IRS1307_time.sem_n.^2);
costs(14) = costs(14) + sum(((expData.IRS1307_time.response_d - scaleIRS1307*IRS1307_time_d).^2)./expData.IRS1307_time.sem_d.^2);
%Note: data_d at t5 has very small sigma (<0.052)
sims.Normal.IRS1307_time=[scaleIRS1307*IRS1307_time_n'; expData.IRS1307_time.time'];
sims.Diabetes.IRS1307_time=[scaleIRS1307*IRS1307_time_d'; expData.IRS1307_time.time'];

% PKB308
tIdx=ismember(sim.n_60.time,expData.PKB308_time.time);
PKB308_time_n = sim.n_60.variablevalues(tIdx,idxPKB308);
PKB308_time_d = sim.d_60.variablevalues(tIdx,idxPKB308);
scalePKB308 = lscov(PKB308_time_n,expData.PKB308_time.response_n, expData.PKB308_time.sem_n.^-2);
costs(15) = sum(((expData.PKB308_time.response_n - scalePKB308.*PKB308_time_n).^2)./expData.PKB308_time.sem_n.^2);
costs(15) = costs(15) + sum(((expData.PKB308_time.response_d - scalePKB308.*PKB308_time_d).^2)./expData.PKB308_time.sem_d.^2);
%Note: data_d at tEnd has very small sigma
sims.Normal.PKB308_time=[scalePKB308*PKB308_time_n'; expData.PKB308_time.time'];
sims.Diabetes.PKB308_time=[scalePKB308*PKB308_time_d'; expData.PKB308_time.time'];

% PKB473
% 10nm insulin 60min normal (n) and diabetes (d)
% baseline not subtracted, ExpData normalised to max 100. Normalised to 10 point in 2p ExpData to get relationship between n and d
% simulation scaled with lscov
% Combination: Multiply expData.time with expData.2p divided by expData.time at 10 min => expData.time(9)=1.
comb_expData.PKB473_time.response_n = expData.PKB473_time_renorm.response_n*expData.PKB473_2p.response_n(2); % normalising the 60 min data to the 10 min 2p data point
comb_expData.PKB473_time.sem_n = expData.PKB473_time_renorm.sem_n*expData.PKB473_2p.response_n(2);
comb_expData.PKB473_time.response_d = expData.PKB473_time_renorm.response_d*expData.PKB473_2p.response_d(2); % normalising the 60 min data to the 10 min 2p data point
comb_expData.PKB473_time.sem_d = expData.PKB473_time_renorm.sem_d*expData.PKB473_2p.response_d(2);

tIdx = ismember(sim.n_60.time, expData.PKB473_time_renorm.time);
PKB473_time_n = sim.n_60.variablevalues(tIdx,idxPKB473);
PKB473_time_d = sim.d_60.variablevalues(tIdx,idxPKB473);
scalePKB473 = lscov(PKB473_time_n,comb_expData.PKB473_time.response_n, comb_expData.PKB473_time.sem_n.^-2);
costs(16)  = sum(((comb_expData.PKB473_time.response_n - scalePKB473.*PKB473_time_n).^2)./comb_expData.PKB473_time.sem_n.^2);
costs(16)  = costs(16)  + sum(((comb_expData.PKB473_time.response_d- scalePKB473.*PKB473_time_d).^2)./comb_expData.PKB473_time.sem_d.^2);
%Note: data at t9 has sigma=0
sims.Normal.PKB473_time=[scalePKB473*PKB473_time_n'; expData.PKB473_time.time];
sims.Diabetes.PKB473_time=[scalePKB473*PKB473_time_d'; expData.PKB473_time.time];

% AS160
tIdx=ismember(sim.n_60.time, expData.AS160_time.time);
AS160_time_n = sim.n_60.variablevalues(tIdx,idxAS160);
AS160_time_d = sim.d_60.variablevalues(tIdx,idxAS160);
scaleAS160 = lscov(AS160_time_n,expData.AS160_time.response_n, expData.AS160_time.sem_n.^-2);
costs(17) =  sum(((expData.AS160_time.response_n - scaleAS160*AS160_time_n).^2)./expData.AS160_time.sem_n.^2);
costs(17) = costs(17) + sum(((expData.AS160_time.response_d - scaleAS160*AS160_time_d).^2)./expData.AS160_time.sem_d.^2);
sims.Normal.AS160_time=[scaleAS160*AS160_time_n'; expData.AS160_time.time'];
sims.Diabetes.AS160_time=[scaleAS160*AS160_time_d'; expData.AS160_time.time'];

% Glucose uptake two point data
% 0 and 100 nM insulin 15 min
% then 0.05 nM glucose 30 min
GLUCOSE_n = [sim.dr_gluc_n0.variablevalues(end,idxGLUCOSE); sim.dr_gluc_n6.variablevalues(idxT30,idxGLUCOSE)];
GLUCOSE_d = [sim.dr_gluc_d0.variablevalues(end,idxGLUCOSE); sim.dr_gluc_d6.variablevalues(idxT30,idxGLUCOSE)];
scaleGLUCOSE = lscov(GLUCOSE_n,expData.GLUCOSE_2p.response_n, expData.GLUCOSE_2p.sem_n.^-2);
costs(18) =  sum(((expData.GLUCOSE_2p.response_n - scaleGLUCOSE.*GLUCOSE_n).^2)./expData.GLUCOSE_2p.sem_n.^2);
costs(18) = costs(18) + sum(((expData.GLUCOSE_2p.response_d - scaleGLUCOSE.*GLUCOSE_d).^2)./expData.GLUCOSE_2p.sem_d.^2);
sims.Normal.GLUCOSE_2p=[scaleGLUCOSE*GLUCOSE_n'; [0 30]];
sims.Diabetes.GLUCOSE_2p=[scaleGLUCOSE*GLUCOSE_d'; [0 30]];

% S6K, 10nM, 60min
% Fold over basal for normal and diabetes data respectively
% Simdata scaled with lscov for n and d respectively
if any(idxS6K)
    tIdx=ismember(sim.n_60.time, expData.S6K_time.time);
    S6K_time_n = sim.n_60.variablevalues(tIdx,idxS6K);
    S6K_time_d = sim.d_60.variablevalues(tIdx,idxS6K);
    scaleS6K = lscov(S6K_time_n, expData.S6K_time.response_n, expData.S6K_time.sem_n.^-2);
    costs(19) = sum(((expData.S6K_time.response_n - scaleS6K*S6K_time_n).^2)./expData.S6K_time.sem_n.^2);
    costs(19) = costs(19) + sum(((expData.S6K_time.response_d - scaleS6K*S6K_time_d).^2)./expData.S6K_time.sem_d.^2);
    sims.Normal.S6K_time=[scaleS6K*S6K_time_n'; expData.S6K_time.time'];
    sims.Diabetes.S6K_time=[scaleS6K*S6K_time_d'; expData.S6K_time.time'];
end

% S6, 10 nM insulin, 60 min
% Fold over basal for normal and diabetes data respectively
% Simdata scaled with lscov for n and d respectively
if any(idxS6)
    tIdx=ismember(sim.n_60.time, expData.S6_time.time);
    S6_time_n = sim.n_60.variablevalues(tIdx,idxS6);
    S6_time_d = sim.d_60.variablevalues(tIdx,idxS6);
    scaleS6 = lscov(S6_time_n, expData.S6_time.response_n, expData.S6_time.sem_n.^-2);
    costs(20) = sum(((expData.S6_time.response_n - scaleS6*S6_time_n).^2)./expData.S6_time.sem_n.^2);
    costs(20) = costs(20) + sum(((expData.S6_time.response_d - scaleS6*S6_time_d).^2)./expData.S6_time.sem_d.^2);
    sims.Normal.S6_time=[scaleS6*S6_time_n'; expData.S6_time.time'];
    sims.Diabetes.S6_time=[scaleS6*S6_time_d'; expData.S6_time.time'];
end

if any(idxERK)
    tIdx=ismember(sim.n_60.time, expData.ERK_time.time);
    ERK_time_n = sim.n_60.variablevalues(tIdx,idxERK);
    ERK_time_d = sim.d_60.variablevalues(tIdx,idxERK);
    scaleERK = lscov(ERK_time_n, expData.ERK_time.mean_n, expData.ERK_time.sem_n.^-2);
    costs(24) = sum(((expData.ERK_time.mean_n - scaleERK*ERK_time_n).^2)./expData.ERK_time.sem_n.^2);
    costs(24) = costs(24) + sum(((expData.ERK_time.mean_d - scaleERK*ERK_time_d).^2)./expData.ERK_time.sem_d.^2);
    sims.Normal.ERK_time=[scaleERK*ERK_time_n'; expData.ERK_time.time'];
    sims.Diabetes.ERK_time=[scaleERK*ERK_time_d'; expData.ERK_time.time'];
end
if any(idxElk1)
    tIdx=ismember(sim.n_60.time, expData.Elk1_time.time);
    Elk1_time_n = sim.n_60.variablevalues(tIdx,idxElk1);
    Elk1_time_d = sim.d_60.variablevalues(tIdx,idxElk1);
    scaleElk1 = lscov(Elk1_time_n, expData.Elk1_time.response_n', expData.Elk1_time.sem_n.^-2);
    costs(25) = sum(((expData.Elk1_time.response_n - scaleElk1*Elk1_time_n').^2)./expData.Elk1_time.sem_n.^2);
    costs(25) = costs(25) + sum(((expData.Elk1_time.response_d - scaleElk1*Elk1_time_d').^2)./expData.Elk1_time.sem_d.^2);
    sims.Normal.Elk1_time=[scaleElk1*Elk1_time_n'; expData.Elk1_time.time'];
    sims.Diabetes.Elk1_time=[scaleElk1*Elk1_time_d'; expData.Elk1_time.time'];
end
if any(idxFOXO)
    tIdx=ismember(sim.n_60.time, expData.Foxo_time.time);
    FOXO_time_n = sim.n_60.variablevalues(tIdx,idxFOXO);
    FOXO_time_d = sim.d_60.variablevalues(tIdx,idxFOXO);
    scaleFOXO = lscov(FOXO_time_n, expData.Foxo_time.response_n, expData.Foxo_time.sem_n.^-2);
    costs(26) = sum(((expData.Foxo_time.response_n - scaleFOXO*FOXO_time_n).^2)./expData.Foxo_time.sem_n.^2);
    costs(26) = costs(26) + sum(((expData.Foxo_time.response_d - scaleFOXO*FOXO_time_d).^2)./expData.Foxo_time.sem_d.^2);
    sims.Normal.FOXO_time=[scaleFOXO*FOXO_time_n'; expData.Foxo_time.time'];
    sims.Diabetes.FOXO_time=[scaleFOXO*FOXO_time_d'; expData.Foxo_time.time'];
end


%% Ad hoc
costs(21)=0;
% GLUT 4 normal data: Basal 5-30 %, 10 nM insulin 10 min > 30-80 %
if sim.n_0.variablevalues(end,idxGLUT4) < 5
    costs(21) = costs(21) + (sim.n_0.variablevalues(end,idxPKB473)-5)^2;
end
if sim.n_0.variablevalues(end,idxGLUT4) > 30
    costs(21) = costs(21) + (sim.n_0.variablevalues(end,idxPKB473)-30)^2;
end
if sim.n_60.variablevalues(idxT10,idxGLUT4) < 30
    costs(21) = costs(21) + (sim.n_60.variablevalues(idxT10,idxPKB473)-30)^2;
end
if sim.n_60.variablevalues(idxT10,idxGLUT4) > 80
    costs(21) = costs(21) + (sim.n_60.variablevalues(idxT10,idxPKB473)-80)^2;
end

costs(21) = tmpCost;
% end Ad hoc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    FINAL CHECKS AND RETURN COMMAND
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sims.Normal.Properties.RowNames={'Sim','Point'};
sims.Diabetes.Properties.RowNames={'Sim','Point'};

costTotal = nansum(costs)+sum(isnan(costs))*limit;

end

