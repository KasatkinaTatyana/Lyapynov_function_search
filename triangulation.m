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
xL = -2.1;
xR = 2.1;
yL = -2.1;
yR = 2.1;

n1 = 0;
n2 = 0;
trAr = [];
% итоговый массив треугольников
T = [];
while true
    assignTriangleAr(n1, n2, trAr);
    appendTriangleAr(T, trAr);
    
    n1 = n1 + 1;
    
    assignTriangleAr(n1, n2, trAr);
    appendTriangleAr(T, trAr);
    
    n1 = n1 - 1;
    n2 = n2 + 1;
    
    assignTriangleAr(n1, n2, trAr);
    appendTriangleAr(T, trAr);
    
    n1 = n1 + 1;
    
    % если ни в какую сторону двигаться нельзя, то триангуляция построена
    if ( boundCond(1) && boundCond(2) && boundCond(3) && boundCond(4) )
        break;
    end
end