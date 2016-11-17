function F = f_sys(x)
global r

F = zeros(1,2);
F(1) = -3*x(1) - 4*x(2) + r*(x(1)^2 - x(2)^2);
F(2) = x(1) + x(2);
end