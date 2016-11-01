function F = assignTriangleAr(n1, n2, trAr)
% 1� ������� ������ - ����� ������������
% 2� � 3� ������� - ���������� 1� ������� (�� x, �� y)
% 4� � 5� ������� - ���������� 2� ������� (�� x, �� y)
% 6� � 7� ������� - ���������� 3� ������� (�� x, �� y)
trAr = [trAr; 0 n1  n2  n1+1  n2     n1+1  n2+1];
trAr = [trAr; 1 n1  n2  n1    n2+1   n1+1  n2+1];
trAr = [trAr; 2 -n1 n2  -n1-1 n2     -n1-1 n2+1];
trAr = [trAr; 3 -n1 n2  -n1   n2+1   -n1-1 n2+1];
trAr = [trAr; 4 -n1 -n2 -n1-1 -n2    -n1-1 -n2-1];
trAr = [trAr; 5 -n1 -n2 -n1   -n2-1  -n1-1 -n2-1];
trAr = [trAr; 6 n1  -n2 n1+1  -n2    n1+1  -n2-1];
trAr = [trAr; 7 n1  -n2 n1    -n2-1  n1+1  -n2-1];
end