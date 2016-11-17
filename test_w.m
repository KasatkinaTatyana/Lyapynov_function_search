Vx0 = 4;
Vx1 = 5;
Vx2 = 1;

x0 = [0 0];
x1 = [1 0];
x2 = [1 1];

w1 = Vx0*(x1(2) - x2(2)) - Vx1*(x0(2) - x2(2)) + Vx2*(x0(2) - x1(2)) / ... 
     (x1(1)*x2(2) - x1(2)*x2(1) - (x0(1)*x2(2) - x2(1)*x0(2)) + ...
     x0(1)*x1(2) - x1(1)*x0(2));
 
w2 = -Vx0*(x1(1) - x2(1)) + Vx1*(x0(1) - x2(1)) - Vx2*(x0(1) - x1(1)) / ...
     (x1(1)*x2(2) - x1(2)*x2(1) - (x0(1)*x2(2) - x2(1)*x0(2)) + ...
     x0(1)*x1(2) - x1(1)*x0(2));
 
 a = [w1*x0(1) + w2*x0(2) - Vx0
      w1*x1(1) + w2*x1(2) - Vx1
      w1*x2(1) + w2*x2(2) - Vx2];
 