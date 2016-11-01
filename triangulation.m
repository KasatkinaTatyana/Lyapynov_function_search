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
xL = -2.1;
xR = 2.1;
yL = -2.1;
yR = 2.1;

n1 = 0;
n2 = 0;
trAr = [];
% �������� ������ �������������
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
    
    % ���� �� � ����� ������� ��������� ������, �� ������������ ���������
    if ( boundCond(1) && boundCond(2) && boundCond(3) && boundCond(4) )
        break;
    end
end