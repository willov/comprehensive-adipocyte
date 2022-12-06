function [S] = CatSims(structs)
%Concatenates multiple structs. If concatenating simulation structs from IQM, you probably want to use CatSim instead.

    S = structs(1);
    for i = 2:length(structs)
        fields = fieldnames(S);
        for k = 1:numel(fields)
            aField     = fields{k}; % EDIT: changed to {}
            if strcmp(aField,'time') % If field is a horizontal 1D array
                S.(aField) = horzcat(S.(aField), structs(i).(aField));
            elseif contains(aField, 'values')%if field is a vertical 1D array
                S.(aField) = vertcat(S.(aField), structs(i).(aField));
            end
        end
    end
end
