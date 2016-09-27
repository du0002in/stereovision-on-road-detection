function [xtmpl, ytmpl, xtmpr, ytmpr, refine_para]=refine_model(leftrightlane,para,sel_ytmpr,sel_ytmpl,D_pre,CorepxR,CorepxL)
global_varibles;
ImgSize2=640;
refine_para=para;
top_row=2*floor((cc_y-40)/2)-1; %make sure cc_y-40 is odd
btm_row=2*ceil((cc_y-480)/2); %make sure cc_y-480 is even
if sum(abs(para(:,2)))==0 %model is hypobola
    if leftrightlane==0
        ytmp=sel_ytmpr+cc_y; %convert to camera coordinate
        Atmp=[1./(ytmp'-D_pre) ytmp' ones(length(ytmp), 1)];
        Btmp=CorepxR';
        ABC=pinv(Atmp)*Btmp;
        A=ABC(1); Bright=ABC(2); Cright=ABC(3);
        Atmp=[ytmp' ones(length(ytmp), 1)];
        Btmp=CorepxL'-A./(ytmp'-D_pre);
        BC=pinv(Atmp)*Btmp;
        Bleft=BC(1); Cleft=BC(2);
        refine_para(1)=A; refine_para(2)=Bleft; refine_para(3)=Cleft; 
        refine_para(4)=Bright; refine_para(5)=Cright; 
    else
        ytmp=sel_ytmpl+cc_y; %convert to camera coordinate
        Atmp=[1./(ytmp'-D_pre) ytmp' ones(length(ytmp), 1)];
        Btmp=CorepxL';
        ABC=pinv(Atmp)*Btmp;
        A=ABC(1); Bleft=ABC(2); Cleft=ABC(3);
        Atmp=[ytmp' ones(length(ytmp), 1)];
        Btmp=CorepxR'-A./(ytmp'-D_pre);
        BC=pinv(Atmp)*Btmp;
        Bright=BC(1); Cright=BC(2);
        refine_para(1)=A; refine_para(2)=Bleft; refine_para(3)=Cleft; 
        refine_para(4)=Bright; refine_para(5)=Cright; 
    end
    if mod(floor(D_pre),2)==0
        ytop=floor(D_pre)-1;
    else ytop=floor(D_pre);
    end
    
    ytmpl=min(ytop,top_row):-2:btm_row;
    xtmpl=A./(ytmpl-D_pre)+Bleft*ytmpl+Cleft;
    ytmpr=min(ytop,top_row):-2:btm_row;
    xtmpr=A./(ytmpr-D_pre)+Bright*ytmpl+Cright;
    %convert to 240x640 image coordinate
    ytmpl=-(-ytmpl+cc_y+1)/2;
    xtmpl=-(xtmpl-cc_x_left);
    ytmpr=-(-ytmpr+cc_y+1)/2;
    xtmpr=-(xtmpr-cc_x_right);
    
    ytmpl(round(xtmpl)<=0)=[];
    xtmpl(round(xtmpl)<=0)=[];
    ytmpl(round(xtmpl)>ImgSize2)=[];
    xtmpl(round(xtmpl)>ImgSize2)=[];
    
    ytmpr(round(xtmpr)<=0)=[];
    xtmpr(round(xtmpr)<=0)=[];
    ytmpr(round(xtmpr)>ImgSize2)=[];
    xtmpr(round(xtmpr)>ImgSize2)=[];
else %model is line
    if leftrightlane==0
        ytmp=sel_ytmpr+cc_y; %convert to camera coordinate
        Atmp=[2*ytmp' ones(length(ytmp), 1)];
        Btmp=-2*CorepxR';
        ef=pinv(Atmp)*Btmp;
        eright=ef(1); fright=ef(2);
        Atmp=[2*ytmp' ones(length(ytmp), 1)];
        Btmp=-2*CorepxL';
        ef=pinv(Atmp)*Btmp;
        eleft=ef(1); fleft=ef(2);
        refine_para(7)=eleft; refine_para(8)=fleft; 
        refine_para(9)=eright; refine_para(10)=fright; 
    else
        ytmp=sel_ytmpl+cc_y; %convert to camera coordinate
        Atmp=[2*ytmp' ones(length(ytmp), 1)];
        Btmp=-2*CorepxL';
        ef=pinv(Atmp)*Btmp;
        eleft=ef(1); fleft=ef(2);
        Atmp=[2*ytmp' ones(length(ytmp), 1)];
        Btmp=-2*CorepxR';
        ef=pinv(Atmp)*Btmp;
        eright=ef(1); fright=ef(2);
        refine_para(7)=eleft; refine_para(8)=fleft; 
        refine_para(9)=eright; refine_para(10)=fright; 
    end
    ytmpl=top_row:-2:btm_row;
    xtmpl=-eleft*ytmpl-fleft/2;
    ytmpr=top_row:-2:btm_row;
    xtmpr=-eright*ytmpr-fright/2;
    %convert to -240x640 image coordinate
    ytmpl=-(-ytmpl+cc_y+1)/2;
    xtmpl=-(xtmpl-cc_x_left);
    ytmpr=-(-ytmpr+cc_y+1)/2;
    xtmpr=-(xtmpr-cc_x_right);
    
    ytmpl(round(xtmpl)<=0)=[];
    xtmpl(round(xtmpl)<=0)=[];
    ytmpl(round(xtmpl)>ImgSize2)=[];
    xtmpl(round(xtmpl)>ImgSize2)=[];
    
    ytmpr(round(xtmpr)<=0)=[];
    xtmpr(round(xtmpr)<=0)=[];
    ytmpr(round(xtmpr)>ImgSize2)=[];
    xtmpr(round(xtmpr)>ImgSize2)=[];
end