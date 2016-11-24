function F = verificate(T, L_fAr)

vertexAr = assignVertexAr(T);
% переменные Vx0, Vx1, Vx2 располагаются в том же порядке, что и вершины в 
% массиве vertexAr
vertexCount = length(vertexAr(:,1));
trCount = length(T(:,1));

derSysAr = [];
for i=1:trCount
    x1 = T(i,2);
    x2 = T(i,4);
    x3 = T(i,6);
    
    y1 = T(i,3);
    y2 = T(i,5);
    y3 = T(i,7);
    
    [minVal, ind_1] = min( (vertexAr(:,1) - x1).^2 + ...
        (vertexAr(:,2) - y1).^2 );
    [minVal, ind_2] = min( (vertexAr(:,1) - x2).^2 + ...
        (vertexAr(:,2) - y2).^2 );
    [minVal, ind_3] = min( (vertexAr(:,1) - x3).^2 + ...
        (vertexAr(:,2) - y3).^2 );
    
    V1 = L_fAr(ind_1);
    V2 = L_fAr(ind_2);
    V3 = L_fAr(ind_3);
    
    A = [x1 y1 1;
         x2 y2 1;
         x3 y3 1
         ];
    B = [V1; V2; V3];
    
    X = inv(A)*B;
    
    w = [X(1) X(2)];
    a = X(3);
    
    derSysAr = [derSysAr; dot(w, f_sys([x1 y1]))];
    derSysAr = [derSysAr; dot(w, f_sys([x2 y2]))];
    derSysAr = [derSysAr; dot(w, f_sys([x2 y2]))];
    
    p = [0.25*x1 + 0.25*x2 + 0.5*x3, 0.25*y1 + 0.25*y2 + 0.5*y3];
    derSysAr = [derSysAr; dot(w, f_sys(p)) ];
end

F = derSysAr;