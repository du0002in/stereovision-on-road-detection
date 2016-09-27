function [left_img_painted, right_img_painted]=paint_img(left_img,right_img,seq)
load 'night_intensity_distribution';
pdf_lane=night_lanemarking_intensity_distribution(seq,:);
pdf_road=night_road_intensity_distribution(seq,:);
csum_pdf_lane=cumsum(pdf_lane);
csum_pdf_road=cumsum(pdf_road);
%for left image
ind_lane=find(left_img==1); ind_road=find(left_img==0);
assigned_intensity=C_assign_intensity(csum_pdf_lane,length(ind_lane));
left_img(ind_lane)=assigned_intensity;
assigned_intensity=C_assign_intensity(csum_pdf_road,length(ind_road));
left_img(ind_road)=assigned_intensity;
left_img_painted=left_img;

%for right image
ind_lane=find(right_img==1); ind_road=find(right_img==0);
assigned_intensity=C_assign_intensity(csum_pdf_lane,length(ind_lane));
right_img(ind_lane)=assigned_intensity;
assigned_intensity=C_assign_intensity(csum_pdf_road,length(ind_road));
right_img(ind_road)=assigned_intensity;
right_img_painted=right_img;