%Adjust the index, so that the rectification is applied to odd rows only
%% Left Image
a1_left(mod(ind_new_left,2)==0)=[];
a2_left(mod(ind_new_left,2)==0)=[];
a3_left(mod(ind_new_left,2)==0)=[];
a4_left(mod(ind_new_left,2)==0)=[];

ind_1_left(mod(ind_new_left,2)==0)=[];
ind_2_left(mod(ind_new_left,2)==0)=[];
ind_3_left(mod(ind_new_left,2)==0)=[];
ind_4_left(mod(ind_new_left,2)==0)=[];
ind_new_left(mod(ind_new_left,2)==0)=[];

%because in C index always starts from 0.
ind_1_left=ind_1_left-1;
ind_2_left=ind_2_left-1;
ind_3_left=ind_3_left-1;
ind_4_left=ind_4_left-1;
ind_new_left=ind_new_left-1;
ind_new_left=ind_new_left';

%% Right Image
a1_right(mod(ind_new_right,2)==0)=[];
a2_right(mod(ind_new_right,2)==0)=[];
a3_right(mod(ind_new_right,2)==0)=[];
a4_right(mod(ind_new_right,2)==0)=[];

ind_1_right(mod(ind_new_right,2)==0)=[];
ind_2_right(mod(ind_new_right,2)==0)=[];
ind_3_right(mod(ind_new_right,2)==0)=[];
ind_4_right(mod(ind_new_right,2)==0)=[];
ind_new_right(mod(ind_new_right,2)==0)=[];

%because in C index always starts from 0.
ind_1_right=ind_1_right-1;
ind_2_right=ind_2_right-1;
ind_3_right=ind_3_right-1;
ind_4_right=ind_4_right-1;
ind_new_right=ind_new_right-1;
ind_new_right=ind_new_right';