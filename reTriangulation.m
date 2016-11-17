function T_new = reTriangulation(T_new, T_old)
T_new = [];
trCnt = length(T_old(:,1));

for i=1:trCnt
    x1 = 0.75 * T_old(i,2);
    y1 = 0.75 * T_old(i,3);
    x2 = 0.75 * T_old(i,4);
    y2 = 0.75 * T_old(i,5);
    x3 = 0.75 * T_old(i,6);
    y3 = 0.75 * T_old(i,7);
    
% debug    
%     x1 = T_old(i,2);
%     y1 = T_old(i,3);
%     x2 = T_old(i,4);
%     y2 = T_old(i,5);
%     x3 = T_old(i,6);
%     y3 = T_old(i,7);
    
    % если одна из вершин - точка (0; 0)
    if ((abs(x1) < 1e-10) && (abs(y1) < 1e-10))
        T_new = [T_new; i x1 y1 x2 y2 0.5*(x2+x3) 0.5*(y2+y3)];
        T_new = [T_new; i x1 y1 x3 y3 0.5*(x2+x3) 0.5*(y2+y3)];
    elseif ((abs(x2) < 1e-10) && (abs(y2) < 1e-10))
        T_new = [T_new; i x2 y2 x1 y1 0.5*(x1+x3) 0.5*(y1+y3)];
        T_new = [T_new; i x2 y2 x3 y3 0.5*(x1+x3) 0.5*(y1+y3)];
    elseif ((abs(x3) < 1e-10) && (abs(y3) < 1e-10))
        T_new = [T_new; i x3 y3 x1 y1 0.5*(x1+x2) 0.5*(y1+y2)];
        T_new = [T_new; i x3 y3 x2 y2 0.5*(x1+x2) 0.5*(y1+y2)];
    % если точка (0; 0) не является вершиной треугольника
    else
        T_new = [T_new; i x1 y1 0.5*(x1+x2) 0.5*(y1+y2) 0.5*(x1+x3) 0.5*(y1+y3)];
        T_new = [T_new; i x2 y2 0.5*(x1+x2) 0.5*(y1+y2) 0.5*(x2+x3) 0.5*(y2+y3)];
        T_new = [T_new; i x3 y3 0.5*(x3+x2) 0.5*(y3+y2) 0.5*(x1+x3) 0.5*(y1+y3)];
        T_new = [T_new; i 0.5*(x2+x3) 0.5*(y2+y3) 0.5*(x1+x2) 0.5*(y1+y2) 0.5*(x1+x3) 0.5*(y1+y3)];
    end
    
end

end