function [nchrset,fitnessmax,mintra,meantra,Anchrset,Afitnessmax,Amintra,Ameantra] = compare(chrset,minval,maxval,T,Achrset,Aminval,Amaxval,AT)
%COMPARE 对比任意两个染色体的genalg结果并输出、绘制平均和最短路径长度图像
%   此处显示详细
    N=1000;   %虚拟并发线程数
    minmat=zeros(N,T);
    meanmat=zeros(N,T);
    Aminmat=zeros(N,AT);
    Ameanmat=zeros(N,AT);
    for thread=1:N
        [nchrset,fitnessmax,mintra,meantra]=genalg(chrset,minval,maxval,T);
        [Anchrset,Afitnessmax,Amintra,Ameantra]=Aga(Achrset,Aminval,Amaxval,AT);
        minmat(thread,:)=mintra;
        meanmat(thread,:)=meantra;
        Aminmat(thread,:)=Amintra;
        Ameanmat(thread,:)=Ameantra;
        thread
    end
    min=mean(minmat);
    means=mean(meanmat);
    Amin=mean(Aminmat);
    Amean=mean(Ameanmat);
    hold off;
    plot(1:T,min,'-.',1:T,means,'-.',1:AT,Amin,'-',1:AT,Amean,'-');
    hold on;
    legend('经典GA最短路径','经典GA平均路径','A*GA最短路径','A*GA平均路径');
    xlabel('迭代次数');
    ylabel('路径长度');
%   hold off;
%   plot(1:T,maxfittra,1:T,meanfittra,1:AT,Amaxfittra,'h-',1:AT,Ameanfittra,'h-');
%   hold on;
end

