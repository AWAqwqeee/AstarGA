function suitable = limit(chro)
%LIMIT chro是单条染色体列向量 suitable为0是失败
%   此处显示详细说明
    suit=0;

map_up=[
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	1	1	1	0	0	1	0	0	0	0	1	1	1	1	1	0	0	0	0
0	1	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	1	0	0	0	0	0	0	0	0	0	0	1	1	0	0	0	0	0	0
0	1	0	0	1	1	1	0	0	0	0	0	1	1	0	0	0	0	0	0
0	1	0	0	1	1	1	0	0	1	1	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	1	1	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	1	1	0	0	0	0	0	0
0	0	0	1	0	1	0	0	0	0	0	0	1	1	0	0	0	0	0	0
0	0	0	0	0	0	0	0	1	1	1	0	0	0	0	0	1	1	1	1
0	0	1	1	0	0	0	0	1	1	1	0	0	0	0	0	1	1	1	1
0	0	1	1	0	0	0	0	1	1	1	0	1	1	0	0	1	1	1	1
0	0	1	1	0	1	0	0	0	0	0	0	1	1	0	0	0	0	0	0
0	0	1	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	1	1	0	0
0	0	0	0	0	0	0	0	1	1	0	0	1	1	0	0	1	1	0	0
0	0	0	0	0	0	0	0	1	1	0	0	1	1	0	0	0	0	0	0
0	0	0	0	0	0	0	0	1	1	0	0	0	0	0	0	0	0	0	0];

    chro2=[1;chro;400];
    len=length(chro2);
    for index = 1:len-1
        test1=chro2(index);
        test2=chro2(index+1);
        t1xy=[floor((test1-1)/20)+1,round(rem(test1-1,20))+1];
        t2xy=[floor((test2-1)/20)+1,round(rem(test2-1,20))+1];
        dist=pdist([t1xy;t2xy]);
        cur=t1xy;
        vec=t2xy-t1xy;
        curv=vec./dist;

        suit=suit+map_up(round(t1xy(2)),round(t1xy(1)));
        suit=suit+map_up(round(t2xy(2)),round(t2xy(1)));
        for curi = 1:0.5:dist
            curnew=cur+curi.*curv;
            suit=suit+map_up(round(curnew(2)),round(curnew(1)));
            %map_up(round(curnew(2)),round(curnew(1)))=map_up(round(curnew(2)),round(curnew(1)))+0.5;
        end
    end
    suitable=~suit;
    %mpup=map_up;
end

