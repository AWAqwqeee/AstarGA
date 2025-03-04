function nchrset = roulette(chrset)
%ROULETTE 轮盘赌选择杂交父代并进行杂交并选择存活，需提前写好适应度函数和约束函数
%   此处显示详细说明
    [rows,chrsum]=size(chrset);%获取染色体、基因个数
    fitval=fitness(chrset);%计算适应度
    fitrate=fitval./sum(fitval);%计算适应度比例（归一化）
    childsum=floor(chrsum.*0.6);%杂交子代个数    可设置
    newchr=zeros(rows,childsum);%子代染色体集初始化
    survchr=zeros(rows,chrsum-childsum);%存活的父代染色体集初始化
    childnow=1;%目前在产生的子代序号
    crochrs=zeros(rows,2);%初始化准备杂交的染色体
    cur=1;%初始化查询光标
    while(childnow<=childsum)%轮盘赌
        r=rand();
        while(r>0)
            r=r-fitrate(cur);
            cur=rem(cur,chrsum)+1;
        end
        crochrs(:,1)=chrset(:,cur);
        r=rand();
        while(r>0)
            r=r-fitrate(cur);
            cur=rem(cur,chrsum)+1;
        end
        crochrs(:,2)=chrset(:,cur);
        crofac=randi([0 1],rows,1);
        newchr(:,childnow)=crochrs(:,1).*crofac+crochrs(:,2).*~crofac;
        if(limit(newchr(:,childnow)))
            childnow=childnow+1;
        end
    end
    [B,I]=sort(fitval,'descend');
    survchr=(chrset(:,I(1:chrsum-childsum)));
    nchrset=[newchr survchr];
end
