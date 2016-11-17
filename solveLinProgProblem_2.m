function [exitflag, L_fAr] = solveLinProgProblem_2(T)
global r

vertexAr = assignVertexAr(T);
% переменные Vx0, Vx1, Vx2 располагаются в том же порядке, что и вершины в 
% массиве vertexAr
vertexCount = length(vertexAr(:,1));
trCount = length(T(:,1));


% матрица A имеет структуру Vx0, Vx1, Vx2, ..., Vxn, C_1_1, C_2_1, C_1_2,
% C_2_2, C_1_3, C_2_3, ... где C_1_i, C_2_i - константы i-го треугольника
A = zeros(vertexCount + 7*trCount, vertexCount + 2*trCount);
b = zeros(vertexCount + 7*trCount, 1);

% Vx0 >= ||x0||
% Vx1 >= ||x1||
% Vx2 >= ||x2||
for i=1:vertexCount
    A(i,i) = -1;
    b(i) = -norm(vertexAr(i,:), 2);
end


for i=1:trCount
    % текущий индекс в матрицах A, b - такой, чтобы начать с 1
    cur_ind = vertexCount + 7*i - 6;
    
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
    
    delta = x1(1)*x2(2) - x1(2)*x2(1) - (x0(1)*x2(2) - x2(1)*x0(2)) + ...
            x0(1)*x1(2) - x1(1)*x0(2);
 
    k1_1 = (x1(2) - x2(2)) / delta;
    k2_1 = -(x0(2) - x2(2)) / delta;
    k3_1 = (x0(2) - x1(2)) / delta;
    
    k1_2 = -(x1(1) - x2(1)) / delta;
    k2_2 = (x0(1) - x2(1)) / delta;
    k3_2 = -(x0(1) - x1(1)) / delta;
    
   B = 2*r;
   f0 = f_sys(x0);
   f1 = f_sys(x1);
   f2 = f_sys(x2);
   
   E0 = 0;
   E1 = norm(x1 - x0,2)*(max(norm(x1 - x0,2), norm(x2 - x0,2)) + ...
                         norm(x1 - x0,2))*B;
   E2 = norm(x2 - x0,2)*(max(norm(x1 - x0,2), norm(x2 - x0,2)) + ...
                         norm(x2 - x0,2))*B;
   
   % определяем, как соотносятся вершины треугольника и вершины массива
   % vertexAr
   [minVal, ind_0] = min( (vertexAr(:,1) - x0(1)).^2 + ...
                          (vertexAr(:,2) - x0(2)).^2 );
   [minVal, ind_1] = min( (vertexAr(:,1) - x1(1)).^2 + ...
                          (vertexAr(:,2) - x1(2)).^2 );
   [minVal, ind_2] = min( (vertexAr(:,1) - x2(1)).^2 + ...
                          (vertexAr(:,2) - x2(2)).^2 );
   
   % |w1| <= C1
   % |w2| <= C2
   % -w1 - C1 <= 0
   A(cur_ind, ind_0) = -k1_1;
   A(cur_ind, ind_1) = -k2_1;
   A(cur_ind, ind_2) = -k3_1;
   A(cur_ind, vertexCount + (i-1)*2 + 1) = -1;
   
   % w1 - C1 <= 0
   A(cur_ind + 1, ind_0) = k1_1;
   A(cur_ind + 1, ind_1) = k2_1;
   A(cur_ind + 1, ind_2) = k3_1;
   A(cur_ind + 1,  vertexCount + (i-1)*2 + 1) = -1;
   
   % -w2 - C2 <= 0
   A(cur_ind + 2, ind_0) = -k1_2;
   A(cur_ind + 2, ind_1) = -k2_2;
   A(cur_ind + 2, ind_2) = -k3_2;
   A(cur_ind + 2,  vertexCount + (i-1)*2 + 2) = -1;
   
   % w2 - C2 <= 0
   A(cur_ind + 3, ind_0) = k1_2;
   A(cur_ind + 3, ind_1) = k2_2;
   A(cur_ind + 3, ind_2) = k3_2;
   A(cur_ind + 3,  vertexCount + (i-1)*2 + 2) = -1;
   
   % w*f(x0) + E0*(C1 + C2) <= -||x0||_2
   A(cur_ind + 4, ind_0) = k1_1*f0(1) + k1_2*f0(2);
   A(cur_ind + 4, ind_1) = k2_1*f0(1) + k2_2*f0(2);
   A(cur_ind + 4, ind_2) = k3_1*f0(1) + k3_2*f0(2);
   A(cur_ind + 4,  vertexCount + (i-1)*2 + 1) = E0;
   A(cur_ind + 4,  vertexCount + (i-1)*2 + 2) = E0;
   b(cur_ind + 4) = -norm(x0,2);
   
   % w*f(x1) + E1*(C1 + C2) <= -||x1||_2
   A(cur_ind + 5, ind_0) = k1_1*f1(1) + k1_2*f1(2);
   A(cur_ind + 5, ind_1) = k2_1*f1(1) + k2_2*f1(2);
   A(cur_ind + 5, ind_2) = k3_1*f1(1) + k3_2*f1(2);
   A(cur_ind + 5,  vertexCount + (i-1)*2 + 1) = E1;
   A(cur_ind + 5,  vertexCount + (i-1)*2 + 2) = E1;
   b(cur_ind + 5) = -norm(x1,2);
   
   % w*f(x0) + E0*(C1 + C2) <= -||x0||_2
   A(cur_ind + 6, ind_0) = k1_1*f2(1) + k1_2*f2(2);
   A(cur_ind + 6, ind_1) = k2_1*f2(1) + k2_2*f2(2);
   A(cur_ind + 6, ind_2) = k3_1*f2(1) + k3_2*f2(2);
   A(cur_ind + 6,  vertexCount + (i-1)*2 + 1) = E2;
   A(cur_ind + 6,  vertexCount + (i-1)*2 + 2) = E2;
   b(cur_ind + 6) = -norm(x2,2);
end

f = zeros(vertexCount + 2*trCount,1);
for i=1:vertexCount
    f(i) = 1;
end

[L_fAr, fval, exitflag, output] = linprog(f, A, b);

exitflag
end