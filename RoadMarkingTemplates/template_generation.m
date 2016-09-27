clear all;

del_w=-30:5:30;
theta=asin(del_w/65);
del_h=floor(26*cos(theta));
for i=1:length(del_w)
    log_ind{i}=[round((65*sin(theta(i))-1)/(26*cos(theta(i))-1)*((1:del_h(i))-1)+1); 1:del_h(i)];
    x=[log_ind{i}(1,:) log_ind{i}(1,:)+1 log_ind{i}(1,:)+2];
    y=[log_ind{i}(2,:) log_ind{i}(2,:) log_ind{i}(2,:)];
    min_x=min(x);
    if min_x<=0
        x=x-min_x+1;
    end
    mask=zeros(max(y), max(x));
    ind=sub2ind(size(mask),y,x);
    mask(ind)=1;
    log_mask{i}=mask;
end

load 'roadmarkingtemplates.mat';
total_counter=1; dupl_char_ind=[4 9 17 22 25];
log_char_all_theta=cell(13,27);
log_rr_all_theta=[]; log_cc_all_theta=[];
for iii=[1,3:7]
% for iii=7
    temp=roadmarkingtemplates{iii};
    sum_col=sum(temp);
    i=1;
    while (sum_col(i)==0)
        i=i+1;
    end
    char=[];j=1;
    while (i<=length(sum_col))
        if sum_col(i)~=0
            char=[char temp(:,i)];
            i=i+1;
            while(sum_col(i)~=0)
                char=[char temp(:,i)];
                i=i+1;
                if (i>length(sum_col))
                    break;
                end
            end
            log_char{j}=char; char=[];
            %         figure;imshow(log_char{j},[]);
            j=j+1;
        else
            i=i+1;
        end
    end
    
    for ii=1:(j-1)
% % % %         if sum(dupl_char_ind==total_counter)~=0
% % % %             total_counter=total_counter+1;
% % % %             continue;
% % % %         end
        [ximg,yimg]=meshgrid(1:size(log_char{ii},2),1:size(log_char{ii},1));

        for i=1:length(theta)
            R=[cos(theta(i)) -sin(theta(i)); sin(theta(i)) cos(theta(i))]';
            xy=round(R*[ximg(:)' 1; yimg(:)' size(log_char{ii},1)/2]);
            xx=xy(1,:); yy=xy(2,:);
            if min(xy(1,:))<=0
                xx=xx-min(xy(1,:))+1;
            end
            if min(xy(2,:))<=0
                yy=yy-min(xy(2,:))+1;
            end
            newimg=zeros(max(yy),max(xx))-1;
            xxyy=sub2ind(size(newimg),yy,xx);
            newimg(xxyy(1:(end-1)))=log_char{ii}(:);
            convnew=conv2(newimg,[1 1 1],'same');
            convnew(convnew>=1)=1; 
            convnew(convnew<=0 & convnew>=-1)=0;
            convnew(convnew<=-2)=-1;
            convnew(newimg==1)=1;
            convnew(newimg==0)=0;
            convnew=convnew(2:5:end,1:2:end);
            convnew(ceil(yy(end)/5)+1,ceil(xx(end)/2))=2;
            sum_conv=sum(convnew);
            convnew(:,sum_conv==-size(convnew,1))=[];
            sum_conv=sum(convnew,2);
            convnew(sum_conv==-size(convnew,2),:)=[];
            [log_rr(i), log_cc(i)]=find(convnew==2);
            convnew(convnew==2)=1;
%             log_ahead{i}=convnew;
            log_char_all_theta{i,total_counter}=double(convnew);
%             figure;imshow(convnew,[]);
        end
%         log_char_all_theta{total_counter}=log_ahead;
        log_rr_all_theta=[log_rr_all_theta log_rr'];
        log_cc_all_theta=[log_cc_all_theta log_cc'];
        log_char_width(total_counter)=size(log_char{ii},2);
        total_counter=total_counter+1;
    end
end
