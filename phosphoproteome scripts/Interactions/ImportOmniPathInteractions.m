interactionsInput = readtable('OmniPathInteractions.csv');
load('../../Data/MS data/uniprotTranslation.mat','uniprotTranslation')
uniprotTranslation=uniprotTranslation(:,{'Entry','Genenamesprimary'});

interactionsOmni=innerjoin(interactionsInput, uniprotTranslation,"LeftKeys","Source","RightKeys","Entry");
interactionsOmni.Source=interactionsOmni.Genenamesprimary;
interactionsOmni.Genenamesprimary=[];
interactionsOmni=innerjoin(interactionsOmni, uniprotTranslation,"LeftKeys","Target","RightKeys","Entry");
interactionsOmni.Target=interactionsOmni.Genenamesprimary;
interactionsOmni.Genenamesprimary=[];
save('interactionsOmni_all', 'interactionsOmni')

interactionsOmni = interactionsOmni(:,{'Source', 'Target','Conf'});
save('interactionsOmni', 'interactionsOmni')