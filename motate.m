function nchrset = motate(chrset,minval,maxval,varargin)
%   pm是变异因子（每条染色体的变异概率），需提前写好约束函数
%   此处提供详细说明
    pm=0.2;
    if(size(varargin)~=0)
        pm=cell2mat(varargin{1});
    end
    [m,n]=size(chrset);
    tmpchrset=chrset;
    for i=1:n
        rdyn=rand();
        if(rdyn<=pm)
            while(rdyn<=pm)
                rd=randi([1 m]);
                jug=randi([0 1]);
                if(jug==1)%完全随机
                    tmpchrset(rd,i)=randi([minval maxval]);
                else
                    tmpchrset(rd,i)=rem(abs(tmpchrset(rd,i)+randi([-10 10])+randi([-10 10]).*20)-1,400)+1;
                    if(tmpchrset(rd,i)==0)
                        tmpchrset(rd,i)=1;
                    end
                end
                if(~limit(tmpchrset(:,i)))
                   tmpchrset(:,i)=chrset(:,i);%变异失败
                   continue;
                end
                rdyn=rand();
            end

        end
    end
    nchrset=tmpchrset;
end
