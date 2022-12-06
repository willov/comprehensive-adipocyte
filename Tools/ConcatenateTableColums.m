
function [t]=ConcatenateTableColums(t1,t2)
% FUNCTION [t]=CONCATENATETABLECOLUMNS(t1,t2)
%The fuction takes two tables, and combines them columnwise.
% Adds shared variables from t2 to t1. In practice, all columns in t2 which is also t1 is added. 
% Any column not in t1 will be ignored. 
t=t1;
for i = 1:length(t1.Properties.VariableNames)
    var=t.Properties.VariableNames{i};
    if any(ismember(t2.Properties.VariableNames,var))
        t.(var)=[t1.(var) t2.(var)];
    end
end
end