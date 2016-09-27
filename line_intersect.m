function [x, y]=line_intersect(abc1, abc2)
a1=abc1(1); b1=abc1(2); c1=abc1(3);
a2=abc2(1); b2=abc2(2); c2=abc2(3);
if a1*b2==a2*b1
    x=[]; y=[];
else
    xy=([a1 b1; a2 b2])\[-c1; -c2];
    x=xy(1); y=xy(2);
end
end