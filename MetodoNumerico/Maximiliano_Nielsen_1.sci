clc
clear
// Entrada : a = vector de coef . del polinomio 'p' ordenados de menor grado a mayor
//         x_0 = numero real
// Salida : b = evaluacion p(x_0)
//          c = evaluacion p'(x_0)
function [b,c] = gHorner(a, x_0)
    n = length(a);
    b = a(n);
    c = b;
    for i = 1:(n-1) do
        b = a(n-i) + x_0*b;
        if n-i > 1 then
            c = b + x_0*c;
        end
    end
endfunction

a = [1 , -2 ,3 ,1];
x_0 = 2;
[p_x0, dp_x0] = gHorner(a, x_0);
disp("p_x0 = " + string(p_x0));
disp("p_x0_exacto = 17");
disp("error_p_x0 = "+string(abs(p_x0-17)/17));
disp("dp_x0 = " + string(dp_x0));
disp("dp_x0_exacto = 22");
disp("error_p_x0 = "+string(abs(p_x0-22)/22));
//Salida
//p_x0 = 17
//p_x0_exacto = 17
//error_p_x0 = 0
//dp_x0 = 22
//dp_x0_exacto = 22
//error_p_x0 = 0.2272727
