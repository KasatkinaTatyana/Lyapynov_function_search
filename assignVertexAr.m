function F = assignVertexAr(T)
%F = [x y] - координаты узлов триангуляции
size = length(T(:,1));
F = [];

for i=1:size
    x1 = T(i,2);
    y1 = T(i,3);
    x2 = T(i,4);
    y2 = T(i,5);
    x3 = T(i,6);
    y3 = T(i,7);
    
    if ( i > 1)
        [minVal, ind] = min( (F(:,1) - x1).^2 + (F(:,2) - y1).^2 );
        if (minVal > 1e-10)
            F = [F; x1 y1];
        end
    else
        F = [F; x1 y1];
    end
    
    [minVal, ind] = min( (F(:,1) - x2).^2 + (F(:,2) - y2).^2 );
    if (minVal > 1e-10)
        F = [F; x2 y2];
    end
    
    [minVal, ind] = min( (F(:,1) - x3).^2 + (F(:,2) - y3).^2 );
    if (minVal > 1e-10)
        F = [F; x3 y3];
    end
end

end