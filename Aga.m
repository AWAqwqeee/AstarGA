function [nchrset,fitAstarmax,maxfittra,meanfittra,varargout] = Aga(chrset,minval,maxval,T,varargin)
%GENALG
% function [nchrset,fitAstarmax,varargout] = genalg(chrset,minval,maxval,T,varargin)
% chrset是染色体集，maxval和minval是上下限统一为单个数，T是杂交代数，最后可选变异因子，默认为0.2。此函数需预先写好limit(chrset)和fitAstar(chrset)
%   此处提供详细说明
    tmpchrset=chrset;
    maxfittra=zeros(T,1);
    meanfittra=zeros(T,1);
    for i=1:T
        tmpchrset=Arou(tmpchrset);
        [~,I]=max(fitAstar(tmpchrset));
        record=tmpchrset(:,I);
        tmpchrset=tmpchrset(:,[1:I-1,I+1:end]);
        if(size(varargin)==0)
            tmpchrset=motate(tmpchrset,minval,maxval);
        else
            tmpchrset=motate(tmpchrset,minval,maxval,varargin);
        end
        tmpchrset=[tmpchrset,record];
        [~,b]=fitAstar(tmpchrset);
        maxfittra(i,1)=min(b);
        meanfittra(i,1)=sum(b)./length(b);
    end
    [fitAstarmax,I]=max(fitAstar(chrset));
    varargout=mat2cell(chrset(:,I),size(chrset,1));
    nchrset=tmpchrset;
    %route(chrset(:,I));
end

