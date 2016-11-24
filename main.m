clear all
close all
clc

global r
r = 0.02;

%условие выхода за границы
global boundCond 
boundCond = [0 0 0 0];
% boundCond[0] = 1: вышли за левую границу области (по x)
% boundCond[1] = 1: вышли за правую границу области (по x)
% boundCond[2] = 1: вышли за нижнюю границу области (по y)
% boundCond[3] = 1: вышли за верхнюю границу области (по y)

% границы области
global xL xR yL yR
xL = -3.1;
xR = 3.1;
yL = -3.1;
yR = 3.1;
global trCnt
trCnt = 0;

n1 = 0;
n2 = 0;

% итоговый массив треугольников
T = [];

K = 0;
while true
    n1 = 0;
    while (n1 <= K)
        n2 = 0;
        while (n2 <= K)           
            if ((n1 + n2) >= K)
                trAr = assignTriangleAr(n1, n2);
                T = appendTriangleAr(T, trAr);
            end                       
            n2 = n2 + 1;          
        end % цикл по n2 
        n1 = n1 + 1;
    end % цикл по n1 
    
    % если ни в какую сторону двигаться нельзя, то триангуляция построена
    if ( boundCond(1) && boundCond(2) && boundCond(3) && boundCond(4) )
        break;
    end
    
    K = K + 1;
end


T1 = [];
T1 = reTriangulation(T1, T);

if (trCnt > 0)
    figure(1);
    hold on
    grid on
    axis([xL-1 xR+1 yL-1 yR+1])
    for i = 1:trCnt
        x1 = T(i,2);
        x2 = T(i,4);
        x3 = T(i,6);
        
        y1 = T(i,3);
        y2 = T(i,5);
        y3 = T(i,7);
        line([x1 x2], [y1 y2]);
        line([x1 x3], [y1 y3]);
        line([x2 x3], [y2 y3]);
    end
end


if (length(T1(:,1)) > 0) 
    figure(2);
    hold on
    grid on
    for i = 1 : length(T1(:,1))
        x1 = T1(i,2);
        x2 = T1(i,4);
        x3 = T1(i,6);
        
        y1 = T1(i,3);
        y2 = T1(i,5);
        y3 = T1(i,7);
        line([x1 x2], [y1 y2]);
        line([x1 x3], [y1 y3]);
        line([x2 x3], [y2 y3]);
    end
end

[exitflag, L_fAr] = solveLinProgProblem_2(T1);
% нумерация вершин в L_fAr совпадает с нумерацией вершин в vertexAr

% визуализация
if exitflag
    vertexAr = assignVertexAr(T1);
    figure(3);
    hold on
    grid on
    axis([xL-1 xR+1 yL-1 yR+1])
    for i = 1 : length(T1(:,1))
        x1 = T1(i,2);
        x2 = T1(i,4);
        x3 = T1(i,6);
        
        y1 = T1(i,3);
        y2 = T1(i,5);
        y3 = T1(i,7);
        
        [minVal, ind_1] = min( (vertexAr(:,1) - x1).^2 + ...
                               (vertexAr(:,2) - y1).^2 );
        [minVal, ind_2] = min( (vertexAr(:,1) - x2).^2 + ...
                               (vertexAr(:,2) - y2).^2 );
        [minVal, ind_3] = min( (vertexAr(:,1) - x3).^2 + ...
                               (vertexAr(:,2) - y3).^2 );
                           
        z1 = L_fAr(ind_1);
        z2 = L_fAr(ind_2);
        z3 = L_fAr(ind_3);
                           
        line([x1 x2], [y1 y2], [z1 z2]);
        line([x1 x3], [y1 y3], [z1 z3]);
        line([x2 x3], [y2 y3], [z2 z3]);
        
    end

    
%     for i=1 : length(L_fAr(:,1))
%     % for i=1 : 40
%         x1 = L_fAr(i,4);
%         x2 = L_fAr(i,6);
%         x3 = L_fAr(i,8);
%         
%         y1 = L_fAr(i,5);
%         y2 = L_fAr(i,7);
%         y3 = L_fAr(i,9);
%         
%         line([x1 x2], [y1 y2], [L_fAr(i,1) L_fAr(i,2)]);
%         line([x1 x3], [y1 y3], [L_fAr(i,1) L_fAr(i,3)]);
%         line([x2 x3], [y2 y3], [L_fAr(i,2) L_fAr(i,3)]);
%     end
end

derSysAr = verificate(T1, L_fAr);