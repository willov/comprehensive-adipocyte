
function [] = PlotInhibitorsGlu(param, model, expInd, expData)

colorNorm = [2 64 167]/256;
colorRapa = [0.3010 0.7450 0.9330];
colorPD = [.9290 0.6940 0.1250];

figure(62)
time=unique([expData.AllTimes 0:0.01:60]);

param(expInd)=exp(param(expInd));

variables = IQMvariables(model);
idxIR=strcmp(variables,'measuredIR');
idxIRS1307   = strcmp(variables,'measuredIRS1307');
idxPKB473    = strcmp(variables,'measuredPKB473');
idxAS160     = strcmp(variables,'measuredAS160');
idxmTORC1    = strcmp(variables,'measuredmTORC1');
idxS6        = strcmp(variables,'measuredS6');
idxmTORC2    = strcmp(variables,'measuredmTORC2');
idxERK       = strcmp(variables,'measuredERK');
idxFOXO      = strcmp(variables,'measuredFOXO');

[simNormal, tmpCost] = SimulateGlucose(param, model, time,1,1);
if tmpCost > 0
    disp('Simulation has crashed.');
    return
end

%% Simulate rapamycin inhibition
RapamycinInhib=table();
RapamycinInhib.mTORC1=simNormal.n_60.variablevalues(:,idxmTORC1)';
RapamycinInhib.PKB473=simNormal.n_60.variablevalues(:,idxPKB473)';
RapamycinInhib.FOXO=simNormal.n_60.variablevalues(:,idxFOXO)';
RapamycinInhib.S6=simNormal.n_60.variablevalues(:,idxS6)';
RapamycinInhib.IR=simNormal.n_60.variablevalues(:,idxIR)';
RapamycinInhib.IRS1307=simNormal.n_60.variablevalues(:,idxIRS1307)';
RapamycinInhib.AS160=simNormal.n_60.variablevalues(:,idxAS160)';

RapamycinDR_IRS1=simNormal.dr_n.IRS1';

pNames=IQMparameters(model);
idxInhib=ismember(pNames, {'kG5a1', 'kG5a2'});
fprintf("Simulating rapamycin inhibition: \n")
for i=1:7
    [sim] = SimulateGlucose(param, model, time, idxInhib,(i*2));
    RapamycinInhib(i+1,:)=table(sim.n_60.variablevalues(:,idxmTORC1)',...
        sim.n_60.variablevalues(:,idxPKB473)',...
        sim.n_60.variablevalues(:,idxFOXO)',...
        sim.n_60.variablevalues(:,idxS6)',...
        sim.n_60.variablevalues(:,idxIR)',...
        sim.n_60.variablevalues(:,idxIRS1307)',...
        sim.n_60.variablevalues(:,idxAS160)');

    RapamycinDR_IRS1=[RapamycinDR_IRS1; sim.dr_n.IRS1'];
    if i==1
        fprintf('%i of %i \n|', i, 7)
    elseif mod(i,50)==0
        fprintf(' %i of %i \n', i, 7)
    else
        fprintf('|')
    end
end
fprintf('\n')

RapamycinInhib(end+1,:)=table(sim.n_60.time,...
    sim.n_60.time,...
    sim.n_60.time,...
    sim.n_60.time,...
    sim.n_60.time,...
    sim.n_60.time,...
    sim.n_60.time);

color=[colorNorm; repmat(colorRapa,height(RapamycinInhib)-1,1)];
subplot(4,3,1)

for i =fliplr(1:height(RapamycinInhib)-1)
    plot(RapamycinInhib.mTORC1(end,:), RapamycinInhib.mTORC1(i,:),'color', color(i,:))
    hold on
end
title('mTORC1a')
ylabel('a.u.')
xlabel('Time (min)')

subplot(4,3,7)
for i =fliplr(1:height(RapamycinInhib)-1)
    plot(RapamycinInhib.PKB473(end,:), RapamycinInhib.PKB473(i,:),'color', color(i,:))
    hold on
end
ylabel({'PKB-473P'; 'a.u.'})
xlabel('Time (min)')

subplot(4,3,10)
for i =fliplr(1:height(RapamycinInhib)-1)
    plot(RapamycinInhib.FOXO(end,:), RapamycinInhib.FOXO(i,:),'color', color(i,:))
    hold on
end
ylabel({'FOXO1-S256P'; 'a.u.'})
xlabel('Time (min)')

figure(63) 
sgtitle('Brännmark - Rapamycin plots')

subplot(7,1,1)
for i =fliplr(1:height(RapamycinInhib)-1)
    plot(RapamycinInhib.S6(end,:), RapamycinInhib.S6(i,:),'color', color(i,:))
    hold on
end
ylabel({'S6-S235/236P, a.u.'})
xlabel('Time (min)')

subplot(7,1,2)
for i =fliplr(1:height(RapamycinInhib)-1)
    plot(RapamycinInhib.IR(end,:), RapamycinInhib.IR(i,:),'color', color(i,:))
    hold on
end
ylabel({'IR-YP, a.u.'})
xlabel('Time (min)')
 
subplot(7,1,3)
conc = 1e-9*[0.001 0.01 0.1 0.3 1 10 100]; % Corresponds to ins_conc2 in simulation function;
for i=1:size(RapamycinDR_IRS1,1)
    semilogx(conc, RapamycinDR_IRS1(i,:),'color',color(i,:))
    hold on
end
ylabel({'IR-YP, %-of max'})
xlabel('[insulin] M')

subplot(7,1,4)
for i =fliplr(1:height(RapamycinInhib)-1)
    plot(RapamycinInhib.IRS1307(end,:), RapamycinInhib.IRS1307(i,:),'color', color(i,:))
    hold on
end
ylabel({'IRS1-307P, a.u.'})
xlabel('Time (min)')

subplot(7,1,5)
for i =fliplr(1:height(RapamycinInhib)-1)
    plot(RapamycinInhib.AS160(end,:), RapamycinInhib.AS160(i,:),'color', color(i,:))
    hold on
end
ylabel({'AS160-T642, a.u.'})
xlabel('Time (min)')

subplot(7,1,7)
for i =fliplr(1:height(RapamycinInhib)-1)
    plot(RapamycinInhib.PKB473(end,:), RapamycinInhib.PKB473(i,:),'color', color(i,:))
    hold on
end
ylabel({'PKB-S473, a.u.'})
xlabel('Time (min)')

%% Simulate torin inhibition
TorinInhib=table();
TorinInhib.mTORC1=simNormal.n_60.variablevalues(:,idxmTORC1)';
TorinInhib.mTORC2=simNormal.n_60.variablevalues(:,idxmTORC2)';
TorinInhib.PKB473=simNormal.n_60.variablevalues(:,idxPKB473)';
TorinInhib.FOXO=simNormal.n_60.variablevalues(:,idxFOXO)';

pNames=IQMparameters(model);
idxInhib=ismember(pNames, {'kG5a1', 'kG5a2', 'kG5c'});

fprintf("Simulating torin inhibition: \n")
for i=1:7
    [sim] = SimulateGlucose(param, model, time, idxInhib,[(i*2)*2 (i*2)*2 (i*2)*1000]); 
    TorinInhib(i+1,:)=table(sim.n_60.variablevalues(:,idxmTORC1)',...
        sim.n_60.variablevalues(:,idxmTORC2)',...
        sim.n_60.variablevalues(:,idxPKB473)',...
        sim.n_60.variablevalues(:,idxFOXO)');
        if i==1
            fprintf('%i of %i \n|', i, 7)
        elseif mod(i,50)==0
            fprintf(' %i of %i \n', i, 7)
        else
            fprintf('|')
        end
end
fprintf('\n')

figure(62)
TorinInhib(end+1,:)=table(sim.n_60.time,...
    sim.n_60.time,...
    sim.n_60.time,...
    sim.n_60.time);

color=[colorNorm; repmat([0.4940 0.1840 0.5560],height(TorinInhib)-1,1)];
subplot(4,3,2)
for i =fliplr(1:height(TorinInhib)-1)
    plot(TorinInhib.mTORC1(end,:), TorinInhib.mTORC1(i,:),'color', color(i,:))
    hold on
end
title('mTORC1a')
ylabel('a.u.')
xlabel('Time (min)')

subplot(4,3,5)
for i =fliplr(1:height(TorinInhib)-1)
    plot(TorinInhib.mTORC2(end,:), TorinInhib.mTORC2(i,:),'color', color(i,:))
    hold on
end
title('mTORC2a')
ylabel('a.u.')
xlabel('Time (min)')

subplot(4,3,8)
for i =fliplr(1:height(TorinInhib)-1)
    plot(TorinInhib.PKB473(end,:), TorinInhib.PKB473(i,:),'color', color(i,:))
    hold on
end
ylabel('a.u.')
xlabel('Time (min)')

subplot(4,3,11)
for i =fliplr(1:height(TorinInhib)-1)
    plot(TorinInhib.FOXO(end,:), TorinInhib.FOXO(i,:),'color', color(i,:))
    hold on
end
ylabel('a.u.')
xlabel('Time (min)')

%% Simulate PD184352 inhibition
PDinhib=table();
PDinhib.ERK=simNormal.n_60.variablevalues(:,idxERK)';
PDinhib.PKB473=simNormal.n_60.variablevalues(:,idxPKB473)';
PDinhib.FOXO=simNormal.n_60.variablevalues(:,idxFOXO)';

pNames=IQMparameters(model);
idxInhib=ismember(pNames, {'kG10basal', 'kG10a1', 'kG10a2'});
fprintf("Simulating PD184352 inhibition: \n")

for i=1:7
    [sim] = SimulateGlucose(param, model, time, idxInhib,(i*2));
    PDinhib(i+1,:)=table(sim.n_60.variablevalues(:,idxERK)',...
        sim.n_60.variablevalues(:,idxPKB473)',...
        sim.n_60.variablevalues(:,idxFOXO)');
        if i==1
            fprintf('%i of %i \n|', i, 7)
        elseif mod(i,50)==0
            fprintf(' %i of %i \n', i, 7)
        else
            fprintf('|')
        end
end
fprintf('\n')

PDinhib(end+1,:)=table(sim.n_60.time,...
    sim.n_60.time,...
    sim.n_60.time);

color=[colorNorm; repmat(colorPD,height(PDinhib)-1,1)];
subplot(4,3,3)
for i =fliplr(1:height(PDinhib)-1)
    plot(PDinhib.ERK(end,:), PDinhib.ERK(i,:),'color', color(i,:))
    hold on
end
title('ERK1/2')
ylabel('a.u.')
xlabel('Time (min)')

subplot(4,3,9)
for i =fliplr(1:height(PDinhib)-1)
    plot(PDinhib.PKB473(end,:), PDinhib.PKB473(i,:),'color', color(i,:))
    hold on
end
ylabel('a.u.')
xlabel('Time (min)')

subplot(4,3,12)
for i =fliplr(1:height(PDinhib)-1)
    plot(PDinhib.FOXO(end,:), PDinhib.FOXO(i,:),'color', color(i,:))
    hold on
end
ylabel('a.u.')
xlabel('Time (min)')

set(figure(62), 'outerposition',[1641 240 900 800], 'PaperType','a4')
set(figure(63), 'outerposition',[0 0 300 1440], 'PaperType','a4')
end