function [fitval,distlist] = fitAstar(chrset)
%FITNESS 示例函数
%   此处显示详细说明
    s=size(chrset);
    c=[ones(1,s(2));chrset;400.*ones(1,s(2))];
    len=s(1)+2;
    distlist=zeros(1,s(2));
    fit=zeros(1,s(2));
    for index = 1:len-1
        test1=c(index,:);
        test2=c(index+1,:);
        for num = 1:s(2)    %index是列索引
            t1xy=[floor((test1(num)-1)/20)+1,round(rem(test1(num)-1,20))+1];
            t2xy=[floor((test2(num)-1)/20)+1,round(rem(test2(num)-1,20))+1];
            distlist(num)=distlist(num)+pdist([t1xy;t2xy]);
            fit(num)=fit(num)+abs((t1xy(1)-t1xy(2))./sqrt(2));
        end
    end
    fitval = 1./(distlist-20.*sqrt(2))+1./fit;
end

