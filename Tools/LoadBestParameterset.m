function [optParam, valid]=LoadBestParameterset(baseFolder, modelName, files,p, costFun, limit)

valid=0;
names=string({files.name});
files(~contains(names, modelName)| contains(names, 'tmp'))=[]; %removes all files not containing "key"

if p==1
    pos = 1;
else
    pos=length(files);
end

if ~isempty(files)
    matches = regexp({files.name}, '.+opt\(([\+\-0-9\.eE]+)\)', 'tokens');
    values = double(string(matches));
    [sorted_values, ind] = sort(values);
    files=files(ind);
    for n=0:length(files)-1 %Some numerical differences might yield nonvalid start guesses.
        if ~isempty(who('-file',[files(pos+p*n).folder '/' files(pos+p*n).name]','optParam'))
            load([files(pos+p*n).folder '/' files(pos+p*n).name],'optParam'); %ascend or descend if min or max
            
            [~,cost]=costFun(optParam);
            if floor(cost*1e4)<=round(limit*1e4) % We accept slightly worse solutions to cope with numerical differences in differnt computers and operating systems.
                valid=1;
                break;
            else
                if ~exist(strrep(files(pos+p*n).folder,'PL','PL_notValid'),'dir')
                    mkdir(strrep(files(pos+p*n).folder,'PL','PL_notValid'))
                end
                movefile([files(pos+p*n).folder '/' files(pos+p*n).name],...
                    [strrep(files(pos+p*n).folder,'PL','PL_notValid') '/' files(pos+p*n).name]);
            end
        end
    end
end
if isempty(files)  || valid==0
    load(FindBestParametersFile(strrep(baseFolder,'-PPL',''), 1, [modelName ', opt']), 'optParam')
end
end