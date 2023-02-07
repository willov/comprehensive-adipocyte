function [ SEM ] = nanSEM(x, dim)
%NANSEM Returns the answer  to SEM = std(x)/sqrt(n), given a matrix x, either as a single
%value or column/row specified by 'dim'. 
if nargin<2
    dim=0;
end


if dim>0
    SEM=nanstd(x,[],dim);
    n=sum(~isnan(x),dim);
else
    SEM=nanstd(reshape(x,1,numel(x)));
    n=sum(sum(~isnan(x)));
end

n(n<=0)=nan;
SEM=SEM./sqrt(n);
end

