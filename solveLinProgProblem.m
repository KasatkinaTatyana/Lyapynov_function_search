function [exitflag, L_fAr] = solveLinProgProblem(T)
global r
% T - массив с треугольниками
% L_fAr - массив вида Vx0 Vx1 Vx2 x0_x x0_y x1_x x1_y x2_x x2_y

L_fAr = [];

% Vx0 Vx1 Vx2 C1 C2

trCnt = length(T(:,1));
for i=1 : trCnt
    xVertex_1 = T(i,2);
    yVertex_1 = T(i,3);
    xVertex_2 = T(i,4);
    yVertex_2 = T(i,5);
    xVertex_3 = T(i,6);
    yVertex_3 = T(i,7);
    
    if (abs(xVertex_1) < 1e-10) && (abs(yVertex_1) < 1e-10)
        % сохраним обозначения вершин x0, x1, x2 из статьи
        % в x0 записываем координаты нулевой вершины, если она принадлежит 
        % текущему треугольнику
        x0 = [xVertex_1 yVertex_1];
        x1 = [xVertex_2 yVertex_2];
        x2 = [xVertex_3 yVertex_3];
    elseif (abs(xVertex_2) < 1e-10) && (abs(yVertex_2) < 1e-10)
        x0 = [xVertex_2 yVertex_2];
        x1 = [xVertex_1 yVertex_1];
        x2 = [xVertex_3 yVertex_3];
    elseif (abs(xVertex_3) < 1e-10) && (abs(yVertex_3) < 1e-10)
        x0 = [xVertex_3 yVertex_3];
        x1 = [xVertex_2 yVertex_2];
        x2 = [xVertex_1 yVertex_1];
    else
        % неважно, какую вершину выбрать за x0
        x0 = [xVertex_1 yVertex_1];
        x1 = [xVertex_2 yVertex_2];
        x2 = [xVertex_3 yVertex_3];
    end
        
    % неравенства Ax <= b
    % Vx0 >= ||x0||
    % Vx1 >= ||x1||
    % Vx2 >= ||x2||
    A = [-1 0 0 0 0;
          0 -1 0 0 0;
          0 0 -1 0 0];
      
    b = [-norm(x0, 2);
         -norm(x1, 2);
         -norm(x2, 2)];
    % |w1| <= C1
    % |w2| <= C2
    delta = x1(1)*x2(2) - x1(2)*x2(1) - (x0(1)*x2(2) - x2(1)*x0(2)) + ...
     x0(1)*x1(2) - x1(1)*x0(2);
 
    k1_1 = (x1(2) - x2(2)) / delta;
    k2_1 = -(x0(2) - x2(2)) / delta;
    k3_1 = (x0(2) - x1(2)) / delta;
    
    k1_2 = -(x1(1) - x2(1)) / delta;
    k2_2 = (x0(1) - x2(1)) / delta;
    k3_2 = -(x0(1) - x1(1)) / delta;
    
    
    
    A = [A;
         -k1_1  -k2_1  -k3_1  -1 0;
          k1_1  k2_1  k3_1  -1 0;
         -k1_2  -k2_2  -k3_2  0 -1;
          k1_2  k2_2  k3_2  0 -1];
      
    b = [b;
         0;
         0;
         0;
         0];
     
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
   
   A = [A;
        k1_1*f0(1)+k1_2*f0(2)  k2_1*f0(1)+k2_2*f0(2)  k3_1*f0(1)+k3_2*f0(2)  E0 E0;
        k1_1*f1(1)+k1_2*f1(2)  k2_1*f1(1)+k2_2*f1(2)  k3_1*f1(1)+k3_2*f1(2)  E1 E1;
        k1_1*f2(1)+k1_2*f2(2)  k2_1*f2(1)+k2_2*f2(2)  k3_1*f2(1)+k3_2*f2(2)  E2 E2;
        ];
   b = [b;
        -norm(x0, 2);
        -norm(x1, 2);
        -norm(x2, 2)
        ];
  
   % целевая функция
   f = [1 1 1 0 0];
   % проверяем, что значение функции в этой вершине уже найдено
%    Aeq = [];
%    beq = [];
%    L_fAr_aux = zeros(1,9);
%    if (length(L_fAr) > 1)
%        % x0
%        [minVal, ind] = min( (L_fAr(:,4) - x0(1)).^2 + (L_fAr(:,5) - x0(2)).^2 );
%        if (minVal < 1e-10)
%            Aeq = [Aeq; 1 0 0 0 0];
%            beq = [beq; L_fAr(ind, 1)];
%            L_fAr_aux(1) = L_fAr(ind, 1);
%            L_fAr_aux(4) = L_fAr(ind, 4);
%            L_fAr_aux(5) = L_fAr(ind, 5);
%        end
%        [minVal, ind] = min( (L_fAr(:,6) - x0(1)).^2 + (L_fAr(:,7) - x0(2)).^2 );
%        if (minVal < 1e-10)
%            Aeq = [Aeq; 1 0 0 0 0];
%            beq = [beq; L_fAr(ind, 2)];
%            L_fAr_aux(1) = L_fAr(ind, 2);
%            L_fAr_aux(4) = L_fAr(ind, 6);
%            L_fAr_aux(5) = L_fAr(ind, 7);
%        end
%        [minVal, ind] = min( (L_fAr(:,8) - x0(1)).^2 + (L_fAr(:,9) - x0(2)).^2 );
%        if (minVal < 1e-10)
%            Aeq = [Aeq; 1 0 0 0 0];
%            beq = [beq; L_fAr(ind, 3)];
%            L_fAr_aux(1) = L_fAr(ind, 3);
%            L_fAr_aux(4) = L_fAr(ind, 8);
%            L_fAr_aux(5) = L_fAr(ind, 9);
%        end
%        
%        % x1
%        [minVal, ind] = min( (L_fAr(:,4) - x1(1)).^2 + (L_fAr(:,5) - x1(2)).^2 );
%        if (minVal < 1e-10)
%            Aeq = [Aeq; 0 1 0 0 0];
%            beq = [beq; L_fAr(ind, 1)];
%            L_fAr_aux(2) = L_fAr(ind, 1);
%            L_fAr_aux(6) = L_fAr(ind, 4);
%            L_fAr_aux(7) = L_fAr(ind, 5);
%        end
%        [minVal, ind] = min( (L_fAr(:,6) - x1(1)).^2 + (L_fAr(:,7) - x1(2)).^2 );
%        if (minVal < 1e-10)
%            Aeq = [Aeq; 0 1 0 0 0];
%            beq = [beq; L_fAr(ind, 2)];
%            L_fAr_aux(2) = L_fAr(ind, 2);
%            L_fAr_aux(6) = L_fAr(ind, 6);
%            L_fAr_aux(7) = L_fAr(ind, 7);
%        end
%        [minVal, ind] = min( (L_fAr(:,8) - x1(1)).^2 + (L_fAr(:,9) - x1(2)).^2 );
%        if (minVal < 1e-10)
%            Aeq = [Aeq; 0 1 0 0 0];
%            beq = [beq; L_fAr(ind, 3)];
%            L_fAr_aux(2) = L_fAr(ind, 3);
%            L_fAr_aux(6) = L_fAr(ind, 8);
%            L_fAr_aux(7) = L_fAr(ind, 9);
%        end
%        
%        %x2
%        [minVal, ind] = min( (L_fAr(:,4) - x2(1)).^2 + (L_fAr(:,5) - x2(2)).^2 );
%        if (minVal < 1e-10)
%            Aeq = [Aeq; 0 0 1 0 0];
%            beq = [beq; L_fAr(ind, 1)];
%            L_fAr_aux(3) = L_fAr(ind, 1);
%            L_fAr_aux(8) = L_fAr(ind, 4);
%            L_fAr_aux(9) = L_fAr(ind, 5);
%        end
%        [minVal, ind] = min( (L_fAr(:,6) - x2(1)).^2 + (L_fAr(:,7) - x2(2)).^2 );
%        if (minVal < 1e-10)
%            Aeq = [Aeq; 0 0 1 0 0];
%            beq = [beq; L_fAr(ind, 2)];
%            L_fAr_aux(3) = L_fAr(ind, 2);
%            L_fAr_aux(8) = L_fAr(ind, 6);
%            L_fAr_aux(9) = L_fAr(ind, 7);
%        end
%        [minVal, ind] = min( (L_fAr(:,8) - x2(1)).^2 + (L_fAr(:,9) - x2(2)).^2 );
%        if (minVal < 1e-10)
%            Aeq = [Aeq; 0 0 1 0 0];
%            beq = [beq; L_fAr(ind, 3)];
%            L_fAr_aux(3) = L_fAr(ind, 3);
%            L_fAr_aux(8) = L_fAr(ind, 8);
%            L_fAr_aux(9) = L_fAr(ind, 9);
%        end
%    end
%    
%    if (length(Aeq) > 1)
%        if (length(Aeq(:,1)) == 3)
%            exitflag = solveReducedLinProgProblem(L_fAr_aux);
%        end
%    end
   
%   [x, fval, exitflag, output] = linprog(f, A, b, Aeq, beq);
   [x, fval, exitflag, output] = linprog(f, A, b);
   
   [(k1_1*x(1) + k2_1*x(2) + k3_1*x(3))*x0(1) + ...
       (k1_2*x(1) + k2_2*x(2) + k3_2*x(3))*x0(2) - x(1);
       (k1_1*x(1) + k2_1*x(2) + k3_1*x(3))*x1(1) + ...
       (k1_2*x(1) + k2_2*x(2) + k3_2*x(3))*x1(2) - x(2);
       (k1_1*x(1) + k2_1*x(2) + k3_1*x(3))*x2(1) + ...
       (k1_2*x(1) + k2_2*x(2) + k3_2*x(3))*x2(2) - x(3);
       ]
   
   if (exitflag ~= 1)
       break
   else
       L_fAr = [L_fAr; x(1) x(2) x(3) x0 x1 x2];
   end
end

end