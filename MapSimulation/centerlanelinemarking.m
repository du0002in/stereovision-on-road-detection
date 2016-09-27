load 'map_cordinate4_full';
Lane_marking_length=2.0;
Lane_marking_gap=4.0;
i=0;
while i<length(map_x)
    sumlen=0;
    while sumlen<Lane_marking_length
        i=i+1;
        if i>=length(map_x)
            break;
        end
        d=sqrt((map_x(i+1)-map_x(i))^2+(map_y(i+1)-map_y(i))^2);
        sumlen=sumlen+d;
    end
    sumlen=0;
    while sumlen<Lane_marking_gap
        i=i+1;
        if i>=length(map_x)
            break;
        end
        d=sqrt((map_x(i+1)-map_x(i))^2+(map_y(i+1)-map_y(i))^2);
        sumlen=sumlen+d;
        Lane_marking_point_x(i,3)=nan;
        Lane_marking_point_x(i,6)=nan;
    end
end