clear all
clc

%условие выхода за границы
global boundCond 
boundCond = [0 0 0 0];
% boundCond[0] = 1: вышли за левую границу области (по x)
% boundCond[1] = 1: вышли за правую границу области (по x)
% boundCond[2] = 1: вышли за нижнюю границу области (по y)
% boundCond[3] = 1: вышли за верхнюю границу области (по y)

% границы области
global xL xR yL yR
xL = -2.1;
xR = 2.1;
yL = -2.1;
yR = 2.1;
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
    % figure(2);
    % hold on
    % grid on
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