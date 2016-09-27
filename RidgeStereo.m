function [kcsditmpl, kcsditmpr]=RidgeStereo(imgl, imgr, ImgSize1, ImgSize2, ix, xfil2, yfil2, GsiCenter,GsiMiddle,GsiCorner)

% ImgSmoothl=zeros(ImgSize1,ImgSize2);
% ImgSmoothr=zeros(ImgSize1,ImgSize2);
img=[imgl(1:2:ImgSize1,:) zeros(ImgSize1/2, 30) imgr(1:2:ImgSize1,:)];
ImgSmoothl=zeros(size(img));

for i=1:(length(ix)-1)
    ImgSmoothl(ix(i):(ix(i+1)-1),:)=conv2(img(ix(i):(ix(i+1)-1),:),xfil2{i},'same');
end
ImgSmoothl=conv2(ImgSmoothl,yfil2,'same');

% ImgSmoothr=conv2(ImgSmoothr,yfil2,'same');
% % % % % 
% % % % % [Gr_Lsd_x, Gr_Lsd_y]=C_Gradient(ImgSmoothl);
% % % % % %compute the structure tensor field Ssdi
% % % % % ssd11=Gr_Lsd_x.^2;
% % % % % ssd12=Gr_Lsd_x.*Gr_Lsd_y;
% % % % % ssd22=Gr_Lsd_y.^2;
% % % % % % Ssdi12=ssd12*(GsiCenter+GsiCorner)+(ssd11+ssd22)*GsiMiddle;
% % % % % % Ssdi11_22=(ssd11-ssd22).*(GsiCenter-GsiCorner);
% % % % % Ssdi12=ssd12*.6306+(ssd11+ssd22)*0.0838;
% % % % % Ssdi11_22=(ssd11-ssd22).*(0.608);
% % % % % %calculate the max eignvector of Ssdi
% % % % % EigenValue_22D12=0.5*(Ssdi11_22+sqrt((Ssdi11_22).^2+4*Ssdi12.^2))./Ssdi12;
% % % % % EigenValue_22D12(isnan(EigenValue_22D12))=0;
% % % % % wpsdiV=sqrt(1./(1+EigenValue_22D12.^2));
% % % % % wpsdiU=wpsdiV.*EigenValue_22D12;
% % % % % %     wpsdiV(isnan(wpsdiV))=0;
% % % % % %     wpsdiU(isnan(wpsdiU))=0;
% % % % % psdi=sign(wpsdiU.*Gr_Lsd_x+wpsdiV.*Gr_Lsd_y);
% % % % % U=psdi.*wpsdiU;
% % % % % V=psdi.*wpsdiV;
% % % % % %ridgeness measurement calculate negtive divergence kcsdi
% % % % % 
% % % % % % tmp=(Gr_Lsd_x.*Gr_Lsd_x-Gr_Lsd_y.*Gr_Lsd_y)*(GsiCorner-GsiCenter)./((Gr_Lsd_x.*Gr_Lsd_x+Gr_Lsd_y.*Gr_Lsd_y)*GsiMiddle+Gr_Lsd_x.*Gr_Lsd_y*(GsiCorner+GsiCenter));
% % % % % % wpsdiV1=2./sqrt(4+(tmp+sqrt(tmp.*tmp+4)).^2);
% % % % % kcsdil=C_Divergence(U,V);
% % % % % toc;
% % % % % tic;
kcsdil=C_DivergenceNew(ImgSmoothl);

kcsditmp=kcsdil;
% Bl=sort(kcsditmp(:),'descend');
% ridge_thl=Bl(5000);

ridge_thl=C_select2(kcsditmp(:),309400);

kcsditmp(kcsdil<ridge_thl)=0;
kcsditmp(kcsdil>=ridge_thl)=1;

% imgtmp=zeros(ImgSize1,size(img,2));
% imgtmp(1:2:ImgSize1,:)=kcsditmp;
% imgtmp(2:2:ImgSize1,:)=kcsditmp;
% imgtmp=kcsditmp;
kcsditmpl=kcsditmp(:,1:ImgSize2);
kcsditmpr=kcsditmp(:,end-ImgSize2+1:end);

end