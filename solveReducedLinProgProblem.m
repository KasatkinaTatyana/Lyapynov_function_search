function F = solveReducedLinProgProblem(L_fAr)
global r
f = [1 1];

x0 = [L_fAr(4) L_fAr(5)];
x1 = [L_fAr(6) L_fAr(7)];
x2 = [L_fAr(8) L_fAr(9)];

delta = x1(1)*x2(2) - x1(2)*x2(1) - (x0(1)*x2(2) - x2(1)*x0(2)) + ...
    x0(1)*x1(2) - x1(1)*x0(2);

k1_1 = (x1(2) - x2(2)) / delta;
k2_1 = -(x0(2) - x2(2)) / delta;
k3_1 = (x0(2) - x1(2)) / delta;

k1_2 = -(x1(1) - x2(1)) / delta;
k2_2 = (x0(1) - x2(1)) / delta;
k3_2 = -(x0(1) - x1(1)) / delta;

w1 = k1_1*L_fAr(1) + k2_1*L_fAr(2) + k3_1*L_fAr(3);
w2 = k1_2*L_fAr(1) + k2_2*L_fAr(2) + k3_2*L_fAr(3);

A = [-1 0;
     0 -1
     ];
 
b = [min(w1, -w1);
     min(w2, -w2)
     ];
 
 % блок неравенств на w * f(x)
 B = 2*r;
 f0 = f_sys(x0);
 f1 = f_sys(x1);
 f2 = f_sys(x2);
 
 E0 = 0;
 E1 = norm(x1 - x0,2)*(max(norm(x1 - x0,2), norm(x2 - x0,2)) + ...
      norm(x1 - x0,2))*B;
 E2 = norm(x2 - x0,2)*(max(norm(x1 - x0,2), norm(x2 - x0,2)) + ...
      norm(x2 - x0,2))*B;
  
 w1*f0(1) + w2*f0(2)
 
 A = [A;
      E1 E1;
      E2 E2
      ];
  
b = [b;
     -norm(x1,2) - w1*f1(1) - w2*f1(2);
     -norm(x2,2) - w1*f2(1) - w2*f2(2)
     ];
 
[x, fval, exitflag, output] = linprog(f, A, b);
F = exitflag;
end