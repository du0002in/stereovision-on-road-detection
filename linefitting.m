function [sampled_xtmp, sampled_ytmp, xtmp, ytmp]=linefitting(kcsditmp,th)
% tic;
ImgSize2=size(kcsditmp,2);
% remove right part of the image
% kcsditmp(:,ImgSize2/2:end)=0;
[indvr, indur]=find(kcsditmp>4.5*th);
indvr=-indvr;
Datapoint=[indvr,indur]';
SubDatapoint=Datapoint;
% toc;
% [line, bedf, tmp]=RANSAC(Datapoint, SubDatapoint);
% toc;
[line, bedf]=P_RANSAC(Datapoint, SubDatapoint);
% sum(line)
% toc;
if sum(line)==0
    xtmp=[];
    ytmp=[];
    return;
end
if sum(abs(bedf(1:4)))==0
    Ytmp=[2*indur(line==1) ones(sum(line),1)];
    d=-2*indvr(line==1);
    cf=(Ytmp'*Ytmp)\(Ytmp')*d;
    findInf=find(cf==Inf);
    if sum(findInf)~=0
        c=1;f=-2*Fullconpoint(1,1);e=0;
    else
        c=cf(1);f=cf(2);e=1;
    end
    ytmp=-21:-1:-240;
    xtmp=(-f-2*e*ytmp)/(2*c);
%     tmpD=0;tmpDO=0;
% % % % %     xtmp2=(xtmp(1:end-1)+xtmp(2:end))/2.0;
% % % % %     ytmp=41:1:479;
% % % % %     xtmp3=zeros(1,length(ytmp));
% % % % %     xtmp3(1:2:end)=xtmp;
% % % % %     xtmp3(2:2:end)=xtmp2;
% % % % %     GoodPoint=find(xtmp3<=ImgSize2 & xtmp3>=1);
% % % % %     xtmp=xtmp3(GoodPoint);
% % % % %     ytmp=ytmp(GoodPoint);
else
    Ytmp=[indvr(line==1).*indvr(line==1) indvr(line==1) indur(line==1) ones(sum(line),1)];
    d=indvr(line==1).*indur(line==1);
    bedftmp=(Ytmp'*Ytmp)\(Ytmp')*d;
    if bedftmp(3)<max(indvr(line==1))
        bedftmp=bedf(1:4);
    end
%     bedftmp=bedf(1:4);
    B=bedftmp(1);E=bedftmp(2);D=bedftmp(3);F=bedftmp(4);
    C=E+B*D; A=F+C*D;
    ytmp=min(floor(D),-21):-1:-240;
    xtmp=A./(ytmp-D)+B*ytmp+C;
%     tmpD=D;tmpDO=bedf(3);
% % % % %     xtmp2=(xtmp(1:end-1)+xtmp(2:end))/2.0;
% % % % %     ytmp=(-ytmp(1)*2-1):1:(-ytmp(end)*2-1);
% % % % %     xtmp3=zeros(1,length(ytmp));
% % % % %     xtmp3(1:2:end)=xtmp;
% % % % %     xtmp3(2:2:end)=xtmp2;
% % % % % %     xtmp3=round(xtmp3);
% % % % %     GoodPoint=find(xtmp3<=ImgSize2 & xtmp3>=1);
% % % % %     xtmp=xtmp3(GoodPoint);
% % % % %     ytmp=ytmp(GoodPoint);
end
% toc
ytmp(round(xtmp)<=0)=[];
xtmp(round(xtmp)<=0)=[];
ytmp(round(xtmp)>ImgSize2)=[];
xtmp(round(xtmp)>ImgSize2)=[];
% toc
% % % figure;plot(indur, indvr,'.');axis equal;
% % % fl=Datapoint(:,line==1);
% % % hold on;plot(fl(2,:),fl(1,:),'r.');
%select 20 points equally spaced based on y or v
% [~,~,indytmp]=intersect(indvr, ytmp);
ind1=intersectAB(indvr(line==1),ytmp);
indytmp=find(ind1==1);
tmp=floor(length(indytmp)/20);
% toc
if tmp>=1
    equal_ind=(1:20)*tmp;
else
    equal_ind=1:length(indytmp);
end
sampled_xtmp=xtmp(indytmp(equal_ind));
sampled_ytmp=ytmp(indytmp(equal_ind));
% toc;
% tmpv=indvr(line==1);
% tmpu=indur(line==1);
end
% indvr(line==1)=[]; indur(line==1)=[];
% delta_x=0.308*ytmp+7.385;
% xtmp_left=xtmp-delta_x;
% xtmp_right=xtmp+delta_x;