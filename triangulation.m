clear all
clc

%������� ������ �� �������
global boundCond 
boundCond = [0 0 0 0];
% boundCond[0] = 1: ����� �� ����� ������� ������� (�� x)
% boundCond[1] = 1: ����� �� ������ ������� ������� (�� x)
% boundCond[2] = 1: ����� �� ������ ������� ������� (�� y)
% boundCond[3] = 1: ����� �� ������� ������� ������� (�� y)

% ������� �������
global xL xR yL yR
xL = -2.1;
xR = 2.1;
yL = -2.1;
yR = 2.1;
global trCnt
trCnt = 0;

n1 = 0;
n2 = 0;

% �������� ������ �������������
T = [];
while true
    trAr = assignTriangleAr(n1, n2);
    T = appendTriangleAr(T, trAr);
    
    n1 = n1 + 1;
    
    trAr = assignTriangleAr(n1, n2);
    T = appendTriangleAr(T, trAr);
    
    n1 = n1 - 1;
    n2 = n2 + 1;
    
    trAr = assignTriangleAr(n1, n2);
    T = appendTriangleAr(T, trAr);
    
    n1 = n1 + 1;
    
    % ���� �� � ����� ������� ��������� ������, �� ������������ ���������
    if ( boundCond(1) && boundCond(2) && boundCond(3) && boundCond(4) )
        break;
    end
end

if (trCnt > 0)
    figure(1);
    hold on
    grid on
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