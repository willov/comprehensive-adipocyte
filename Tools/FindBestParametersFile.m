function [path,dirStruct] = FindBestParametersFile(pathToSearch, useSub, key)
%FINDBESTPARAMETERSFILE finds the file containing the best parameter set
%based on the name. Assumes that the cost is given in parenthesis. 
%   [path,name] = FindBestParametersFile(path)

if nargin<3, key='.'; end
if nargin<2, useSub=0; end

%%
assert(exist(pathToSearch,'dir')>0,'Path does not exist. Are you in the right directory?');
if useSub
    pathToSearch=[pathToSearch '/**/*.mat'];
end
files=dir(pathToSearch);
assert(length(files)>0,'No files to search in the provided path')
names=string({files.name});
if length(files)>1
    files(cellfun(@isempty,regexp(names, key)))=[]; %removes all files not containing "key"
    
    names=string({files.name});
    values=double(regexprep(names,{'.*\(','\).*'},{'',''}));
    assert(any(~isnan(values)),'No files contains valid names. Make sure the files has the cost in parenthesis, e.g opt(12.34)');
    [~,bestInd]=min(values);
    name=files(bestInd).name;
    path=[files(bestInd).folder '/' name];
    dirStruct=files(bestInd);
else
    assert(contains(names,key), "File does not contain provided key")
    path =[files.folder '/' files.name];
    dirStruct=files;
end
end

