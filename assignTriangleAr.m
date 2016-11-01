function F = assignTriangleAr(n1, n2, trAr)
% 1й элемент строки - номер треугольника
% 2й и 3й элемент - координаты 1й вершины (по x, по y)
% 4й и 5й элемент - координаты 2й вершины (по x, по y)
% 6й и 7й элемент - координаты 3й вершины (по x, по y)
trAr = [trAr; 0 n1  n2  n1+1  n2     n1+1  n2+1];
trAr = [trAr; 1 n1  n2  n1    n2+1   n1+1  n2+1];
trAr = [trAr; 2 -n1 n2  -n1-1 n2     -n1-1 n2+1];
trAr = [trAr; 3 -n1 n2  -n1   n2+1   -n1-1 n2+1];
trAr = [trAr; 4 -n1 -n2 -n1-1 -n2    -n1-1 -n2-1];
trAr = [trAr; 5 -n1 -n2 -n1   -n2-1  -n1-1 -n2-1];
trAr = [trAr; 6 n1  -n2 n1+1  -n2    n1+1  -n2-1];
trAr = [trAr; 7 n1  -n2 n1    -n2-1  n1+1  -n2-1];
end