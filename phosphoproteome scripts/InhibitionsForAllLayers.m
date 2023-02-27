function []=InhibitionsForAllLayers(coreModel)

choice = input('Do you want to rerun the calculations, or only plot the figure? \nRerunning the calculations will be slow. \nOnly plot? \nPress ctrl+C to cancel \n(y/n):  ', 's');
if lower(choice)=='n'
    files = [dir('Results/empty*'); dir('Results/layer*')];
    layers=double(regexprep(string({files.name}),{'[a-zA-Z\.]'},{''}));
    [~,idx]=sort(layers);
    files=files(idx); % Sort files based on layer

    x=[0 15/60 30/60 1 2 5 10 20 60];

    simulateDiabetes=0;
    plotNormal=0;
    plotDiabetes=0;
    plotInhib=0;
    plotAllInhib=1;

    % Initialize results tables
    inhibitions = strcat('I', cellstr(string(round(100-100./2.^(1:5)))));
    direction_MK=array2table(zeros(length(files),length(inhibitions)+3), 'VariableNames', ['Layer', 'Rel', 'Additions', inhibitions]);
    direction_MK.DataSelection(:)={'none'}; % Add a column for which data selection have been used
    direction_MK=direction_MK(:,[1:2 end 3:end-1]); %Rearrange the table
    direction_LY=direction_MK;
    direction_MK_layer = direction_MK;
    direction_LY_layer = direction_LY;

    disp('Calculating prediction accuracy for each layer, this may take a while')

    for i = 1:length(files)
        fprintf('Done with %i of %i layers\n', i, length(files))
        try
            load([files(i).folder '/' files(i).name],'list','expData','structure','parameters', 'layer', 'rel', 'layerStructure')

            if min(list.Conf)<0 && rel>=0 %If the layer contains datadriven interactions
                rel=-1;
            end

            direction_MK{i,1:2} = [layer,rel];
            direction_LY{i,1:2} = [layer,rel];
            direction_MK_layer{i,1:2} = [layer,rel];
            direction_LY_layer{i,1:2} = [layer,rel];

            dataSelection = unique(list.DataSelection); % Only works for the specific type of data division
            dataSelection = dataSelection(1);
            direction_MK.DataSelection(i)=dataSelection; direction_LY.DataSelection(i)=dataSelection; direction_LY_layer.DataSelection(i)=dataSelection; direction_MK_layer.DataSelection(i)=dataSelection;

            if ~contains(files(i).name,'empty') && ~isempty(layerStructure) % If the current layer has additions
                modelName=strrep(files(i).name, '.mat','');
                modelName=regexprep(modelName, ' .+','');

                if ~isempty(structure)
                    GenerateModel(structure,modelName,pwd,100);
                    IQMmakeMEXmodel(IQMmodel([modelName '.txt']))
                end

                [~, totalInhib, inhibitionValuesDirection]= PlotFinalModel(parameters, x, expData, modelName, coreModel, list, simulateDiabetes, plotNormal,plotDiabetes,plotInhib, plotAllInhib);

                additions =unique(regexprep(layerStructure.Target,{'_p','_U'},{'',''}));
                direction_MK.Additions(i)=length(additions); direction_LY.Additions(i)=length(additions);
                direction_MK{i,5:end} = totalInhib{'direction_MK',:};
                direction_LY{i,5:end} = totalInhib{'direction_LY',:};

                direction_MK_layer(i,:)=direction_MK(i,:);
                direction_LY_layer(i,:)=direction_LY(i,:);

                nInhibs = length(totalInhib.Properties.VariableNames);
                inhibitionValuesDirection(~ismember(inhibitionValuesDirection.Protein,additions),:)=[];
                direction_MK_layer{i,5:end}=nansum(inhibitionValuesDirection{:,2:nInhibs+1})./sum(~isnan(inhibitionValuesDirection{:,2:nInhibs+1}));
                direction_LY_layer{i,5:end}=nansum(inhibitionValuesDirection{:,nInhibs+2:end})./sum(~isnan(inhibitionValuesDirection{:,nInhibs+2:end}));

                delete([modelName '.txt'])
                delete([modelName '.mex*'])

            else
                disp('Empty layer, skipping')
                direction_MK{i,5:end} = nan;
                direction_LY{i,5:end} = nan;
                direction_MK_layer{i,5:end} = nan;
                direction_LY_layer{i,5:end} = nan;
            end
        catch err
            disp(getReport(err))
        end
    end
    save('./Results/post-all-inhibs.mat') % Saving for the possibility to run the plotting without rerunning the calculations
else
    load('./Results/post-all-inhibs-precalculated.mat', 'direction_MK', 'direction_LY' )
end

%% Plot the results
PlotResults(direction_MK, direction_LY, 9,1) % Cumulative

%% Final cleanup
figFiles = dir('./**/*.fig');
for i =1:length(figFiles)
    delete([figFiles(i).folder '/' figFiles(i).name])
end
end

function [] = PlotResults(direction_MK, direction_LY, fig, m)
direction_MK.Additions=cumsum(direction_MK.Additions);
direction_LY.Additions=cumsum(direction_LY.Additions);

direction_MK(~any(direction_MK{:,[1:2 5:end]}~=0,2),:)=[];
direction_LY(~any(direction_LY{:,[1:2 5:end]}~=0,2),:)=[];

direction_MK = sortrows(direction_MK,'Layer','ascend');
direction_LY = sortrows(direction_LY,'Layer','ascend');

% Handle empty layers (removing,renaming layers)
if any(direction_MK.Additions==0)
    direction_MK(direction_MK.Additions==0,:)=[];
    direction_LY(direction_LY.Additions==0,:)=[];
else
    direction_MK(all(isnan(direction_MK{:,5:end}),2),:)=[];
    direction_LY(all(isnan(direction_LY{:,5:end}),2),:)=[];
end

direction_MK.Layer=(1:height(direction_MK))';
direction_LY.Layer=(1:height(direction_LY))';

maxLayer=max(direction_LY.Layer);
cStart = [8 69 148]/256;
cEnd = [230 238 247]/256;
layerGradient = [linspace(cStart(1),cEnd(1),maxLayer)', linspace(cStart(2),cEnd(2),maxLayer)', linspace(cStart(3),cEnd(3),maxLayer)'];

cStart=[0.3 0.3 0.3];
cEnd=[0.9 0.9 0.9];
relGradient = [0 0 0; [linspace(cStart(1),cEnd(1),21)', linspace(cStart(2),cEnd(2),21)', linspace(cStart(3),cEnd(3),21)']];

entries=cellstr(strcat(strrep(direction_MK.Properties.VariableNames(end-4:end),'I',''), '% inhibition'));

firstNonResponder=direction_MK.Layer(find(strcmp(direction_MK.DataSelection,'Nonresponders'),1));

figure(fig)
clf
colorNorm = [2 64 167]/256;
colorInhib = [32 237 178]/256; %[0.4660    0.6740    0.1880]

cStart = colorNorm; %[8 69 148]/256;
cEnd = colorInhib;
inhibColors = [linspace(cStart(1),cEnd(1),6)', linspace(cStart(2),cEnd(2),6)', linspace(cStart(3),cEnd(3),6)'];
inhibColors(1,:)=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VERSION WITH DATA IN SAME GRAPH

ylower = min(45, min(movmean(direction_LY{:,5:end}*100,m,1,'omitnan')-6,[], 'all'));

subplot(2,2,1)

respondersMK=direction_MK(strcmp(direction_MK.DataSelection, 'Responders'),:);
nonrespondersMK=direction_MK(strcmp(direction_MK.DataSelection, 'Nonresponders'),:);

direction_MK.Rel(strcmp(direction_MK.DataSelection, 'Nonresponders'))=direction_MK.Rel(strcmp(direction_MK.DataSelection, 'Nonresponders'))-21;

[~,ind]=unique(direction_MK.Rel);
for i = ind' % Plot lines corresponding to gradual decrease in interaction reliability
    plot([direction_MK.Layer(i), direction_MK.Layer(i)]-0.5, [-5 105],'color', [0.9 0.9 0.9], 'LineWidth',1, 'HandleVisibility','off')
    hold on
end

[~,ind]=unique(respondersMK.Rel);
for i = ind' % Plot lines corresponding to gradual decrease in interaction reliability
    area([respondersMK.Layer(i)-0.5 max(respondersMK.Layer(respondersMK.Rel==respondersMK.Rel(i)))+0.5],[ylower+4  ylower+4],'FaceColor',relGradient(respondersMK.Rel(i)+2,:),'EdgeColor','none', 'HandleVisibility','off')
end
[~,ind]=unique(nonrespondersMK.Rel);
for i = ind' % Plot lines corresponding to gradual decrease in interaction reliability
    area([nonrespondersMK.Layer(i)-0.5 max(nonrespondersMK.Layer(nonrespondersMK.Rel==nonrespondersMK.Rel(i)))+0.5],[ylower+4 ylower+4],'FaceColor',relGradient(nonrespondersMK.Rel(i)+2,:),'EdgeColor','none', 'HandleVisibility','off')
end

plot([firstNonResponder, firstNonResponder]-0.5, [-5 105],'color', 'k', 'LineWidth',2, 'HandleVisibility','off')

[~,ind]=unique(direction_MK);
for i = ind' % Plot lines corresponding to gradual decrease in interaction reliability
    area([direction_MK.Layer(i)-0.5 max(direction_MK.Layer(direction_MK.Rel==direction_MK.Rel(i)))+0.5],[ylower+2 ylower+2],'FaceColor',layerGradient(direction_MK.Layer(i),:),'EdgeColor','none', 'HandleVisibility','off')
end

p = direction_MK{:,5:end}*100;
pmm=movmean(p,m,1, 'omitnan');
h = plot(direction_MK.Layer, pmm, 'linewidth',2); % Plot percentages of correct prediction
set(h, {'color'}, num2cell(inhibColors,2))

% Set labels etc
title('PKB inhibitor')
ylabel('% correct direction')
xlabel('Layer')
set(gca,'FontSize', 15)
axis([direction_MK.Layer(1), direction_MK.Layer(end), ylower,100])
% yticks((0:5)*10+50)
yticks(0:25:100)
xt = ((0:4)*ceil(((max(direction_MK.Layer)-min(direction_MK.Layer))/6)/10)*10+ceil(min(direction_MK.Layer)/10)*10);
xticks([min(direction_MK.Layer) xt max(direction_MK.Layer)])
box off
set(gca, 'FontSize',20)

legend(entries,'location','best')

% LY-inhibitions
subplot(2,2,2)

respondersLY=direction_LY(strcmp(direction_LY.DataSelection, 'Responders'),:);
nonrespondersLY=direction_LY(strcmp(direction_LY.DataSelection, 'Nonresponders'),:);

direction_LY.Rel(strcmp(direction_LY.DataSelection, 'Nonresponders'))=direction_LY.Rel(strcmp(direction_LY.DataSelection, 'Nonresponders'))-21;

[~,ind]=unique(direction_LY.Rel);
for i = ind' % Plot lines corresponding to gradual decrease in interaction reliability
    plot([direction_LY.Layer(i), direction_LY.Layer(i)]-0.5, [-5 105],'color', [0.9 0.9 0.9], 'LineWidth',1, 'HandleVisibility','off')
    hold on
end
[~,ind]=unique(respondersLY.Rel);
for i = ind' % Plot areas corresponding to gradual decrease in interaction reliability
    area([respondersLY.Layer(i)-0.5 max(respondersLY.Layer(respondersLY.Rel==respondersLY.Rel(i)))+0.5],[ylower+4  ylower+4],'FaceColor',relGradient(respondersLY.Rel(i)+2,:),'EdgeColor','none', 'HandleVisibility','off')
end
[~,ind]=unique(nonrespondersLY.Rel);
for i = ind' % Plot lines corresponding to gradual decrease in interaction reliability
    area([nonrespondersLY.Layer(i)-0.5 max(nonrespondersLY.Layer(nonrespondersLY.Rel==nonrespondersLY.Rel(i)))+0.5],[ylower+4 ylower+4],'FaceColor',relGradient(nonrespondersLY.Rel(i)+2,:),'EdgeColor','none', 'HandleVisibility','off')
end
plot([firstNonResponder, firstNonResponder]-0.5, [-5 105],'color', 'k', 'LineWidth',2, 'HandleVisibility','off')
[~,ind]=unique(direction_LY);
for i = ind' % Plot areas corresponding to expansion iteration
    area([direction_LY.Layer(i)-0.5 max(direction_LY.Layer(direction_LY.Rel==direction_LY.Rel(i)))+0.5],[ylower+2 ylower+2],'FaceColor',layerGradient(direction_LY.Layer(i),:),'EdgeColor','none', 'HandleVisibility','off')
end

p = direction_LY{:,5:end}*100;
pmm=movmean(p,m,1, 'omitnan');
h = plot(direction_LY.Layer, pmm,'linewidth',2); % Plot percentages of correct prediction
set(h, {'color'}, num2cell(inhibColors,2))

% Set labels etc
title('PI3K inhibitor')
ylabel('% correct direction')
xlabel('Layer')
set(gca,'FontSize', 15)
axis([direction_LY.Layer(1), direction_LY.Layer(end), ylower,100])
yticks(0:25:100)
xt = ((0:4)*ceil(((max(direction_LY.Layer)-min(direction_LY.Layer))/6)/10)*10+ceil(min(direction_LY.Layer)/10)*10);
xticks([min(direction_LY.Layer) xt max(direction_LY.Layer)])
box off
set(gca, 'FontSize',20)

subplot(2,2,3)
maxAdditions=max(direction_LY.Additions);

[~,ind]=unique(direction_LY.Rel);
for i = ind' % Plot lines corresponding to gradual decrease in interaction reliability
    plot([direction_LY.Layer(i), direction_LY.Layer(i)]-0.5, [1 maxAdditions],'color', [0.9 0.9 0.9], 'LineWidth',1, 'HandleVisibility','off')
    hold on
end
plot([firstNonResponder, firstNonResponder]-0.5, [0.1 maxAdditions],'color', 'k', 'LineWidth',2, 'HandleVisibility','off')
plot(direction_LY.Layer, direction_LY.Additions, 'k', 'linewidth',2)
title('Total additions (cumulative)')
xlabel('Layer')
box off
ylabel('Number of additions (log)')
set(gca, 'YScale', 'log')
set(gca, 'FontSize',20)
axis([direction_LY.Layer(1), direction_LY.Layer(end), 0.5, maxAdditions])
xt = ((0:4)*ceil(((max(direction_LY.Layer)-min(direction_LY.Layer))/6)/10)*10+ceil(min(direction_LY.Layer)/10)*10);
xticks([min(direction_LY.Layer) xt max(direction_LY.Layer)])
yticks(10.^(0:3))

%%% PLOT LEGENDS
subplot(2,4,7)
relGradientFlip = flipud(relGradient);
for i = 1:size(relGradientFlip,1) % Plot lines corresponding to gradual decrease in interaction reliability
    area([i-0.5 i+0.5]-2,[ylower+2 ylower+2],'FaceColor',relGradientFlip(i,:),'EdgeColor','none', 'HandleVisibility','off')
    hold on
end
axis tight
box off
ylim([ylower,100])
xt = ((0:3)*ceil(((size(relGradient,1)-1)/6)/5)*5+5);
xticks([-1 xt-1 20])
xticklabels([flip(xt) 0 -1])

title('Legend for number of primary sources')

subplot(2,4,8)
[~,ind]=unique(direction_MK);
for i = ind' % Plot lines corresponding to gradual decrease in interaction reliability
    area([direction_MK.Layer(i)-0.5 max(direction_MK.Layer(direction_MK.Rel==direction_MK.Rel(i)))+0.5],[ylower+2 ylower+2],'FaceColor',layerGradient(direction_MK.Layer(i),:),'EdgeColor','none', 'HandleVisibility','off')
    hold on
end
box off

axis([direction_LY.Layer(1), direction_LY.Layer(end), ylower,100])
xt = ((0:4)*ceil(((max(direction_LY.Layer)-min(direction_LY.Layer))/6)/10)*10+ceil(min(direction_LY.Layer)/10)*10);
xticks([min(direction_LY.Layer) xt max(direction_LY.Layer)])

title('Legend for layer')

% Set final properties of the figure
set(figure(fig), 'outerposition',[166 1 1722 1396], 'PaperType','a4')
exportgraphics(figure(fig), '../Fig. 10 Inhibitors.pdf', 'ContentType','vector')
close all
end