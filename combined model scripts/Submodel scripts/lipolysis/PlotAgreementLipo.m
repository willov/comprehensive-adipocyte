function [] = PlotAgreementLipo(optParam, modelName, res, useHSL, predictHSL, baseFolder)
clear mex

if nargin<2, modelName='lipolysis'; end
if nargin<3, res=1; disp('Running in low-res mode. Set res=0.01 for the same resolution as in the original scripts'), end
if nargin<4 , useHSL=1; end
if nargin<5 || isempty(predictHSL), predictHSL=0; end
if nargin<6, baseFolder=''; end

[model,data, ~, ~, nParams, expInd] = Init_lipo(modelName, useHSL);

stimulusHighRes=table();
ins=log10(unique([data.InVitro.FA.Ins; data.InVitro.Glycerol.Ins]));
stimulusHighRes.Ins=[0 10.^(ins(2):res:ins(end)) 0]';
stimulusHighRes.Iso=[0.01*ones(height(stimulusHighRes)-1,1); 0]; %10 nM = 0.01 µM

pNames=IQMparameters(model);
if length(optParam)==nParams-1 % If one diabetes parameter is missing
    optParam=[optParam(1:end-1) 1 optParam(end)]; %assumes that the missing parameter is the one from the glucose model
elseif length(optParam)==nParams-2 % if two diabetes parameters is missing
    optParam=[optParam 1 1];
end
diabInd = contains(pNames,'diab');
optParam(expInd) = exp(optParam(expInd));
diab = optParam(diabInd); %if the vector is shorter than the index of diab_reest, use the last value
optParam(diabInd)=[];

Best=simulateInVitro(model, optParam, diab, stimulusHighRes, 0);
if strcmp(modelName, 'lipolysis') && ~useHSL
    in = load(FindBestParametersFile(baseFolder, 1, 'HSL, cost, min'), 'optParam');
    HSLParams = in.optParam;

    if length(HSLParams)==nParams-1 % If one diabetes parameter is missing
        HSLParams=[HSLParams(1:end-1) 1 HSLParams(end)];
    elseif length(HSLParams)==nParams-2 % if two diabetes parameters is missing
        HSLParams=[HSLParams 1 1];
    end

    HSLParams(expInd) = exp(HSLParams(expInd));
    diab=HSLParams(diabInd);
    HSLParams(diabInd)=[];
    BestHSL=simulateInVitro(model, HSLParams, diab, stimulusHighRes, 0);
    Best.Normal.HSL = BestHSL.Normal.HSL;
end

InVitrotmp.high=Best.Normal;
InVitrotmp.low=Best.Normal;

InVitrotmpD.high=Best.Diabetes;
InVitrotmpD.low=Best.Diabetes;

noDiab = [1 1];
TSBest = SimulateInVivo(optParam, model, data.InVivo, noDiab, 0);
TSBest.FA{:,2:end}=TSBest.FA{:,2:end}./TSBest.FA{1,2:end}; 
TStmp.Gly.low=TSBest.Gly;
TStmp.Gly.high=TSBest.Gly;
TStmp.FA.low=TSBest.FA;
TStmp.FA.high=TSBest.FA;

TSBestD=TSBest;
TSBestD.Gly{:,2:end}=nan;
TSBestD.FA{:,2:end}=nan;
TStmpD.Gly.low=TSBestD.Gly;
TStmpD.Gly.high=TSBestD.Gly;
TStmpD.FA.low=TSBestD.FA;
TStmpD.FA.high=TSBestD.FA;

baseFolder=strrep(baseFolder, 'Results','Results-PPL');
files=dir(sprintf('%s/**/*%s*.mat', baseFolder, modelName));

optParams=[];
for i = fliplr(1:length(files))
    load([files(i).folder '/' files(i).name],'optParam');
    if length(optParam)==nParams-1 % If one diabetes parameter is missing
        optParam=[optParam(1:end-1) 1 optParam(end)];
    elseif length(optParam)==nParams-2 % if two diabetes parameters is missing
        optParam=[optParam 1 1];
    end
    optParams(i,:)=optParam;
end
optParams=unique(optParams,'rows');
if size(optParams,1)>0
    fprintf('\nSimulating the uncertainty for the experiments of the lipolysis submodel:\n')
end
for i = 1:size(optParams,1) 
        optParam=optParams(i,:);
        optParam(expInd)=exp(optParam(expInd));
        diab=optParam(diabInd);
        optParam(diabInd)=[];

        [Tmp]=simulateInVitro(model, optParam, diab, stimulusHighRes);
        InVitrotmp.low{:,:}=min(InVitrotmp.low{:,:}, Tmp.Normal{:,:});
        InVitrotmp.high{:,:}=max(InVitrotmp.high{:,:}, Tmp.Normal{:,:});
        
        TmpTS = SimulateInVivo(optParam, model, data.InVivo, noDiab, 0);
        if isfield(Tmp,'Diabetes')
            InVitrotmpD.low{:,:}=min(InVitrotmpD.low{:,:}, Tmp.Diabetes{:,:});
            InVitrotmpD.high{:,:}=max(InVitrotmpD.high{:,:}, Tmp.Diabetes{:,:});
            
            TmpTSD = SimulateInVivo(optParam, model, data.InVivo, diab, 0);
            TmpTSD.FA{:,2:end}=TmpTSD.FA{:,2:end}./TmpTS.FA{1,2:end};
            
            TStmpD.Gly.low{:,:}=min(TStmpD.Gly.low{:,:}, TmpTSD.Gly{:,:});
            TStmpD.Gly.high{:,:}=max(TStmpD.Gly.high{:,:}, TmpTSD.Gly{:,:});
            TStmpD.FA.low{:,:}=min(TStmpD.FA.low{:,:}, TmpTSD.FA{:,:});
            TStmpD.FA.high{:,:}=max(TStmpD.FA.high{:,:}, TmpTSD.FA{:,:});
        end
        TmpTS.FA{:,2:end}=TmpTS.FA{:,2:end}./TmpTS.FA{1,2:end};
        TStmp.Gly.low{:,:}=min(TStmp.Gly.low{:,:}, TmpTS.Gly{:,:});
        TStmp.Gly.high{:,:}=max(TStmp.Gly.high{:,:}, TmpTS.Gly{:,:});
        TStmp.FA.low{:,:}=min(TStmp.FA.low{:,:}, TmpTS.FA{:,:});
        TStmp.FA.high{:,:}=max(TStmp.FA.high{:,:}, TmpTS.FA{:,:});
        
    if i==1
        fprintf('%i of %i \n|',i,size(optParams,1))
    elseif mod(i,50)==0
        fprintf(' %i of %i \n',i,size(optParams,1))
    else
        fprintf('|')
    end
end
fprintf('\n')

%% Do the plotting
%Setup simulations and data (in vitro)
allInVitroData=struct();
allInVitroData.Normal=data.InVitro;

InVitro.Normal=ConcatenateTableColums(Best.Normal, InVitrotmp.low(:,3:end));
InVitro.Normal=ConcatenateTableColums(InVitro.Normal, InVitrotmp.high(:,3:end));
allInVitroData.Diabetes=data.InVitro_diabetes;

InVitro.Diabetes=ConcatenateTableColums(Best.Diabetes, InVitrotmpD.low(:,3:end));
InVitro.Diabetes=ConcatenateTableColums(InVitro.Diabetes, InVitrotmpD.high(:,3:end));

TS.FA.Diabetes=ConcatenateTableColums(TSBestD.FA, TStmpD.FA.low(:,2:end));
TS.FA.Diabetes=ConcatenateTableColums(TS.FA.Diabetes, TStmpD.FA.high(:,2:end));

%Setup simulations and data (in vivo)
TS.Gly.Normal=ConcatenateTableColums(TSBest.Gly, TStmp.Gly.low(:,2:end));
TS.Gly.Normal=ConcatenateTableColums(TS.Gly.Normal, TStmp.Gly.high(:,2:end));

TS.FA.Normal=ConcatenateTableColums(TSBest.FA, TStmp.FA.low(:,2:end));
TS.FA.Normal=ConcatenateTableColums(TS.FA.Normal, TStmp.FA.high(:,2:end));

InVitro.Normal(:,~ismember(InVitro.Normal.Properties.VariableNames,{'Ins','Iso','HSL', 'FA', 'Glycerol','PKB473','PKB308','Reesterification'}))=[];
InVitro.Diabetes(:,~ismember(InVitro.Diabetes.Properties.VariableNames,{'Ins','Iso','HSL', 'FA', 'Glycerol','PKB473','PKB308', 'Reesterification'}))=[];
Best.Normal(:,~ismember(Best.Normal.Properties.VariableNames,{'Ins','Iso','FA', 'HSL', 'Glycerol','PKB473','PKB308', 'Reesterification'}))=[];

%% Plot
if predictHSL
    InVitroHSL=InVitro;
    InVitroHSL.Normal(:,~ismember(InVitroHSL.Normal.Properties.VariableNames,{'Ins','Iso','HSL'}))=[];
    InVitroHSL.Diabetes(:,~ismember(InVitroHSL.Diabetes.Properties.VariableNames,{'Ins','Iso','HSL'}))=[];
    InVitroHSL.Normal.HSL(:,2:end)=[];
    InVitroHSL.Diabetes.HSL(:,2:end)=[];

    PlotInVitro(InVitroHSL, allInVitroData, {'Normal', 'Diabetes'}, 53)
    set(figure(53), 'outerposition',[37 461 2424 716], 'PaperType','a4')
    set(gca,'FontSize', 15)
else
    if ~useHSL
        InVitro.Normal.HSL=[];
        InVitro.Diabetes.HSL=[];
    end

    data.InVivo.Fig3Epi=[];
    TS.Gly.Normal.Fig3Epi=[];
    pos = [1 4 7 nan 10];

    PlotInVivo(data, TS.Gly.Normal, 11,'Glycerol',[], pos, [4 3])
    PlotInVitro(InVitro, allInVitroData, {'Normal', 'Diabetes'}, 11) % This is the non-overlapping plot

    if height(stimulusHighRes)<800
        fprintf('\n\nNote: running with a low resolution might yield slightly incorrect bounds for the reesterification\n')
    end
end
end

