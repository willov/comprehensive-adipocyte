function [] = PlotInVitro(responses, expData, experiments, id, posOverride)

if nargin<5, posOverride=0; end
if nargin<3 || isempty(experiments), experiments = fieldnames(responses); end

colorNorm = [2 64 167]/256;
colorRapa = [0.3010 0.7450 0.9330];
colorWort = [0.2344  0.6484    0.2891];
colorTorin = [129, 25, 133]/256;
colorAkti = [168, 88, 42]/256;
colorPD = [.9290 0.6940 0.1250];
colorDiab = [0.6350 0.0780 0.1840];

figure(id);

if posOverride ~=0
    m=3; n=2; pos = posOverride;
elseif width(responses.Normal)==3
    m=1; n=2; pos=2;
    hold on
elseif ismember('Diabetes', experiments) && any(contains(responses.Normal.Properties.VariableNames, 'PKB473'))
    m=3; n=3; pos=[2 3 5 6 8 9];
elseif any(ismember({'Rapamycin','Torin', 'PD'}, experiments))
    m=2; n=3; pos=1:6;
elseif ismember('Diabetes', experiments) && ismember('Reesterification', responses.Diabetes.Properties.VariableNames) % For plotting the original diabetes plot
    m=2; n=2; pos=[1 2 3 4];
elseif any(contains(responses.Normal.Properties.VariableNames, 'Ins_')) % for plotting only the input
    [m,n] = CloseToSquare(width(responses.Normal)-2);
    pos=1:width(responses.Normal)-2;
else
    m=3; n=2; pos=[2 4 6];
end

for i =1:length(experiments)
    switch experiments{i}
        case 'Normal'
            c = colorNorm;
        case 'PD'
            c = colorPD;
        case 'Rapamycin'
            c = colorRapa;
        case 'Wort'
            c = colorWort; 
        case 'Torin'
            c = colorTorin;
        case 'Akti'
            c = colorAkti;
        case 'Diabetes'
            c = colorDiab;
        otherwise
            c = [0 0 0];
    end

    response=responses.(experiments{i});

    if isfield(expData,experiments{i})
        data=expData.(experiments{i});
    else
        data=[];
    end

    PlotExperiment(m,n,pos,response, data, c, colorDiab)
end
PlotInput(m,n,pos,responses.Normal)
set(figure(id), 'outerposition',[0 0 2560 1440], 'PaperType','a4')
end

function []=PlotExperiment(m,n,pos,response, data, c, colorDiab)

set(0,'DefaultLineLineWidth',2)

for j=3:width(response)
    variable=response.Properties.VariableNames{j};
    if strcmp(response.Properties.VariableNames{j},'u_C')
        subplot(m,n,pos(j-2))
        x=response{:,'u_cAMP'};
        y=response{:,'u_C'};
        x=reshape(x,numel(x),1);
        y=reshape(y,numel(y),1);

        [~,ind]=sort(x);
        plot(x(ind), y(ind),'o')
    else

        if isequal(c, colorDiab) && ismember(response.Properties.VariableNames{j},{'PKB473','PKB308', 'HSL','Glycerol','FA'})
            l='--';
        elseif size(response,2)==3 && strcmp(response.Properties.VariableNames{3},'HSL')
            l='--';
        elseif strcmp(response.Properties.VariableNames{j},'PKB308')
            l='--';
        else
            l='-';
        end

        if isfield(data,variable)
            PlotSubplot(data.(variable), response(:,[1 2 j]),m,n,pos(j-2), c, l)
        else
            PlotSubplot([], response(:,[1 2 j]),m,n,pos(j-2), c, l)
        end
    end
end
end

function []=PlotSubplot(data, sim, m,n,ind,c, l)
x=sim.Ins*1e-9;
if x(1)==0
    x(1)=x(2)/4;
end
if x(end)==0
    x(end)=x(end-1)*10;
end

subplot(m,n,ind)

y=sim{:,end};
if size(y,2)>1 && isequal(y(:,2), y(:,3))
    y=y(:,2);
end

if size(y,2)==3
    h  = fill([x(1)/1.5 x(1)/1.5 x(1)*1.5 x(1)*1.5],[y(1, 2) y(1, 3) y(1, 3) y(1, 2)  ], c);
    set(h,'facealpha',0.5,'edgealpha',0);
    hold on
    h2 = fill([x(2:end-1); flipud(x(2:end-1))], [y(2:end-1,2); flipud(y(2:end-1,3))],c);
    set(h2,'facealpha',0.5,'edgealpha',0);
    h3  = fill([x(end)/1.5 x(end)/1.5 x(end)*1.5 x(end)*1.5],[y(end, 2) y(end, 3) y(end, 3) y(end, 2)  ], c);
    set(h3,'facealpha',0.5,'edgealpha',0);
elseif size(y,2)>3
    plot([x(1)/1.5 x(1)*1.5],[y(1, :); y(1, :)],'-','Color',c,'linewidth', 2.5);
    hold on
    semilogx(x(2:end-1),y(2:end-1,:),l,'Color',c)
    plot([x(end)/1.5 x(end)*1.5],[y(end, :); y(end, :)],'-','Color',c,'linewidth', 2.5);
end

if any(ismember(size(y,2),[1 3]))
    plot([x(1)/1.5 x(1)*1.5],[y(1, 1) y(1, 1)],'-','Color',c,'linewidth', 2.5);
    hold on
    plot(x(2:end-1), y(2:end-1,1),l, 'color', c,'linewidth', 2.5);
    plot([x(end)/1.5 x(end)*1.5],[y(end, 1) y(end, 1)], '-','Color',c,'linewidth', 2.5);
end

if ~isempty(data)
    insData=data.Ins*1e-9;
    insData([1 end])=[x(1) x(end)]; % adjust for zeros in logscale (use same as in simulation)
    if strcmp(l,'-')
        errorbar(insData, data.Mean, data.SEM,'o','MarkerFaceColor', 'auto', 'linewidth', 2.5,'capsize',12,'Color',c)
    elseif strcmp(l, '--')
        errorbar(insData, data.Mean, data.SEM,'^','MarkerFaceColor', 'white', 'linewidth', 2.5,'capsize',12,'Color',c)
    end
end
end
function [] = PlotInput(m,n,pos, response)

for i = 1:width(response)-2
    subplot(m,n,pos(i))
    axis('tight');

    if contains (response.Properties.VariableNames{i+2},{'u_C'})
        title(response.Properties.VariableNames{i+2})
        set(gca,'XScale','log', 'XMinorTick','on','FontSize', 20)
        axis('tight');
        box off
        return
    elseif contains (response.Properties.VariableNames{i+2},{'u_'})
        ylabel({sprintf('%s,',response.Properties.VariableNames{i+2})})
    elseif contains(response.Properties.VariableNames{i+2},'HSL')
        ylabel({'HSL-S552P'; '% of iso only'}) 
    elseif contains(response.Properties.VariableNames{i+2},'PKB473')
        ylabel({'PKB-T473P'; '% of iso only'}) 
    elseif contains(response.Properties.VariableNames{i+2},'PKB308')
        ylabel({'PKB-S308P'; '% of iso only'}) 
    elseif  contains(response.Properties.VariableNames{i+2},'IR_')
        ylabel(['Response, ' response.Properties.VariableNames{i+2}])
    elseif  contains (response.Properties.VariableNames{i+2},{'Reesterification'})
        ylabel('% Reesterification')
    elseif  contains (response.Properties.VariableNames{i+2},{'Glycerol'})
        ylabel({sprintf('Released %s,',lower(response.Properties.VariableNames{i+2})); '% of iso only'})
    else
        ylabel({sprintf('Released %s,',response.Properties.VariableNames{i+2}); '% of iso only'})
    end
    xlabel('[Insulin] (M)')

    axsPatch=findobj(gca,'type','patch');
    axsLine=findobj(gca,'type','line');
    axsErrorbar = findobj(gca, 'type', 'errorbar');

    if ~isempty(axsLine)
        YData=horzcat(axsLine.YData);
    end
    if ~isempty(axsPatch)
        YData=vertcat(axsPatch.YData)';
    end
    if ~isempty(axsErrorbar)
        YData=[horzcat(axsErrorbar.YData)-horzcat(axsErrorbar.YNegativeDelta) horzcat(axsErrorbar.YData)+horzcat(axsErrorbar.YPositiveDelta)];
    end

    if  contains (response.Properties.VariableNames{i+2},{'Reesterification'})
        minimum = 10;
    else
        minimum = min(YData);
    end
    maximum = max(YData);

    x=response.Ins*1e-9;
    x(1)=x(2)/4/1.5;
    plot(x([1 end-1]), [minimum-0.1*(maximum-minimum) minimum-0.1*(maximum-minimum)],'-','color', [0.7 0.7 0.7],'linewidth', 5) % iso stimulus bar
    set(gca,'XScale','log', 'XMinorTick','on','FontSize', 20)

    if  contains (response.Properties.VariableNames{i+2},{'Reesterification'})
        ylim([0, 100])
    else
        axis('tight');
    end
    a=axis;
    axis([a(1), a(2), minimum-0.175*(maximum-minimum), a(4)])
    xticks([1e-12 1e-9 1e-6])
    box off
end
end