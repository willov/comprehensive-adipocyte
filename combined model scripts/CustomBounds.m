function [lb, ub] = CustomBounds(pNames,lb,ub)

%% Adiponectin special bounds
lb(ismember(pNames,'kAVpip')) = 20e-6;
lb(ismember(pNames,'kAVcell')) = (10.8e-15)*0.8;
ub(ismember(pNames,'kAVpip')) = 60e-6;
ub(ismember(pNames,'kAVcell')) = (214.6e-15)*0.8;

%% Lipolysis special bounds
lb(strcmp(pNames,'phe_effect'))=0.6;
lb(ismember(pNames,{'min1','min2','min3'}))=0;
lb(ismember(pNames,{'n1','n2','n3'}))=0.5;
lb(strcmp(pNames,'isoscale'))=8;
lb(strcmp(pNames,'EC501'))=log10(0.5);
lb(strcmp(pNames,'EC502'))=-5;
lb(strcmp(pNames,'EC503'))=-5;

ub(strcmp(pNames,'phe_effect'))=1;
ub(ismember(pNames,{'min1','min2','min3'}))=100;
ub(ismember(pNames,{'n1','n2','n3'}))=3;
ub(ismember(pNames,'isoscale'))=12;
ub(strcmp(pNames,'EC501'))=log10(1.1);
ub(strcmp(pNames,'EC502'))=3;
ub(strcmp(pNames,'EC503'))=3;

end