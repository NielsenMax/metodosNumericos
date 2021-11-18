clc()
clear()
xdel()
funcprot(0); 
//======= E1 
// A
//Funcion que aplica regla del Trapecio para calcular la integral definida de a a b de f
// Entrada: a,b = extremos de integración; f = función de Scilab
// Salida: w = aproximación de la integral de f en [a,b] por la Regla del Trapecio
function w = Trapecio(a,b,f)
    w = ((b-a)/2)*(f(a)+f(b))
endfunction

//Funcion que aplica la regla de Simpson para calcular la integral definida de a a b de f
// Entrada: a,b = extremos de integración; f = función de Scilab
// Salida: w = aproximación de la integral de f en [a,b] por la Regla de Simpson
function w = Simpson(a,b,f)
    h = (b-a)/2
    w = (h/3)*(f(a)+4*f(a+h)+f(b))
endfunction

T_1 = Trapecio(1, 2,log);
S_1 = Simpson(1, 2, log);
real_1 = intg(1, 2, log);
err_T1 = abs(T_1-real_1)/real_1;
err_S1 = abs(S_1-real_1)/real_1;
/*
cot_T = -(h^3**f''(c))/12 = -(1/12)*(-1/c^2) <= 1/12
cot_S = -(h^5*f^(4)(c))/90 = -(0.5^5/90)*(-6/c^4) <= 1/480
*/


deff("y=f1_2(x)", "y=x^(1/3)");
T_2 = Trapecio(0,0.1, f1_2);
S_2 = Simpson(0,0.1, f1_2);
real_2 = intg(0,0.1, f1_2);
err_T2 = abs(T_2-real_2)/real_2;
err_S2 = abs(S_2-real_2)/real_2;
/*
cot_T = -(h^3*f''(c))/12 = -(0.1^3/12)*(-2/(9*c^(5/3))) <= 0.0008596
cot_S = -(h^5*f^(4)(c))/90 = -(0.05^5/90)*(-80/(81*c^(11/3))) <= 0.0000159
*/

deff("y=f1_3(x)", "y=sin(x)^2");
T_3 = Trapecio(0, %pi/3, f1_3);
S_3 = Simpson(0,%pi/3, f1_3);
real_3 = intg(0,%pi/3, f1_3);
err_T3 = abs(T_3-real_3)/real_3;
err_S3 = abs(S_3-real_3)/real_3;
/*
cot_T = -(h^3*f''(c))/12 = -((%pi/3)^3/12)*(2*cos(2*c)) <= 0 ???
cot_S = -(h^5*f^(4)(c))/90 = -((%pi/6)^5/90)*(-8*cos(2*c)) <= 0 ???
*/

//================== E2
//Funcion que aplica regla del Trapecio compuesta para calcular la integral definida de a a b de f
// Entrada: a,b = extremos de integración; n = cantidad de subintervalos; f = función de Scilab
// Salida: w = aproximación de la integral de f en [a,b] por la Regla del Trapecio
function w = TrapecioCompuesto(a, b, n, fx)
    h = (b-a)/n;
    w = 0.5*(fx(a)+fx(b));
    for i=1:n-1
        xi = a+i*h;
        w = w+fx(xi);
    end
    w = w*h;
endfunction
deff("y=f2_1(x)", "y=1/x");
T2_1 = TrapecioCompuesto(1,3,4,f2_1);
R2_1 = intg(1,3, f2_1);
deff("y=f2_2(x)", "y=x^3");
T2_2 = TrapecioCompuesto(0,2,4,f2_2);
R2_2 = intg(0,2, f2_2);
deff("y=f2_3(x)", "y=x*(1+x^2)^(1/2)");
T2_3 = TrapecioCompuesto(0,3,6,f2_3);
R2_3 = intg(0,3, f2_3);
deff("y=f2_4(x)", "y=sin(%pi*x)");
T2_4 = TrapecioCompuesto(0,1,8,f2_4);
R2_4 = intg(0,1, f2_4);
deff("y=f2_5(x)", "y=x*sin(x)");
T2_5 = TrapecioCompuesto(0,2*%pi,8,f2_5);
R2_5 = intg(0,2*%pi, f2_5);
deff("y=f2_6(x)", "y=x^2*exp(x)");
T2_6 = TrapecioCompuesto(0,1,8,f2_6);
R2_6 = intg(0,1, f2_6);
//================== E3
//Funcion que aplica la regla de Simpson para calcular la integral definida de a a b de f
// Entrada: a,b = extremos de integración; n = cantidad de intervalos; f = función de Scilab
// Salida: w = aproximación de la integral de f en [a,b] por la Regla de Simpson
function w = SimpsonCompuesto(a, b, n, f)
    h = (b-a)/n;
    w = (f(a)+f(b));
    for i=1:n-1
        w = w+f(a+i*h)*(2*modulo(i,2)+2);
    end
    w = w * (h/3);
endfunction
//deff("y=f2_1(x)", "y=1");
T3_1 = SimpsonCompuesto(1,3,4,f2_1);
R3_1 = intg(1,3, f2_1);
//deff("y=f2_2(x)", "y=x^3");
T3_2 = SimpsonCompuesto(0,2,4,f2_2);
R3_2 = intg(0,2, f2_2);
//deff("y=f2_3(x)", "y=x*(1+x^2)^(1/2)");
T3_3 = SimpsonCompuesto(0,3,6,f2_3);
R3_3 = intg(0,3, f2_3);
//deff("y=f2_4(x)", "y=sin(%pi*x)");
T3_4 = SimpsonCompuesto(0,1,8,f2_4);
R3_4 = intg(0,1, f2_4);
//deff("y=f2_5(x)", "y=x*sin(x)");
T3_5 = SimpsonCompuesto(0,2*%pi,8,f2_5);
R3_5 = intg(0,2*%pi, f2_5);
//deff("y=f2_6(x)", "y=x^2*exp(x)");
T3_6 = SimpsonCompuesto(0,1,8,f2_6);
R3_6 = intg(0,1, f2_6);
disp("aprox",T3_1);
disp("real",R3_1);
disp("aprox",T3_2);
disp("real",R3_2);
disp("aprox",T3_3);
disp("real",R3_3);
disp("aprox",T3_4);
disp("real",R3_4);
disp("aprox",T3_5);
disp("real",R3_5);
disp("aprox",T3_6);
disp("real",R3_6);
//======================= E 4
deff("y=f4(x)", "y=(x+1)^(-1)");
T4 = TrapecioCompuesto(0, 1.5, 10, f4);
S4 = TrapecioCompuesto(0, 1.5, 10, f4);
R4 = 0.9262907;
// ===================== E5
//Funcion que aplica regla del Trapecio compuesta con 2 intervalos para calcular la integral definida de a a b y de c a d de f
// Entrada: a,b = extremos de integración exteriores;
// c,d = extremos de integracion interiores;
// n = cantidad de subintervalos; f = función de Scilab
// Salida: w = aproximación de la integral de la integral de f en [c(y),d(y)], en [a,b] por la Regla del Trapecio
function w = TrapecioExt(a,b,c,d,f)
    h = (b-a)/2;
    h2 = (d-c)/2;
    y = [a, a+h, b];
    x = [c, c+h2, d];
    w = h*h2*(0.5*(0.5*f(x(1),y(1))+f(x(1),y(2))+0.5*f(x(1),y(3)))+(0.5*f(x(2),y(1))+f(x(2),y(2))+0.5*f(x(2),y(3)))+0.5*(0.5*f(x(3),y(1))+f(x(3),y(2))+0.5*f(x(3),y(3))));
endfunction

deff("z=f5(x,y)", "z=sin(x+y)");
T5 = TrapecioExt(0,1,0,2,f5);

//=================== E6
//Funcion que aplica regla del Trapecio compuesta para calcular la integral definida de a a b de f(c,x)
// Entrada: a,b = extremos de integración; n = cantidad de subintervalos; f = función de Scilab; c = real a aplicar como primer parametro de f
// Salida: w = aproximación de la integral de f en [a,b] por la Regla del Trapecio
function w = TrapecioCompuestoFD(a, b, n, fx, c)
    h = (b-a)/n;
    w = 0.5*(fx(c, a)+fx(c, b));
    for i=1:n-1
        xi = a+i*h;
        w = w+fx(xi);
    end
    w = w*h;
endfunction

//Funcion que aplica regla del Trapecio compuesta para calcular la integral definida de a a b y de c(x) a d(x) de f
// Entrada: a,b = extremos de integración exteriores;
// c,d = funciones para los extremos de integracion interiores;
// n = cantidad de subintervalos para la integral externa; 
//m = cantidad de intervalos para la integral interna;
//f = función de Scilab
// Salida: w = aproximación de la integral de la integral de f en [c(y),d(y)], en [a,b] por la Regla del Trapecio
function w = DobleTrap(a,b,c,d,n,m,f)
    h = (b-a)/n;
    //deff("z=fxa(y)","z=f("+string(a)+",y)");
    //deff("z=fxb(y)","z=f("+string(b)+",y)");
    w = 0.5*(TrapecioCompuestoFD(c(a),d(a),m,f,a)+TrapecioCompuestoFD(c(b),d(b),m,f, b));
    for i=1:n-1
        xi = a+i*h;
        //deff("z=fxi(y)","z=f("+string(xi)+",y)");
        w = w + TrapecioCompuestoFD(c(xi),d(xi),m,f, xi);
    end
    w = w*h;
endfunction

//Funcion que aplica regla de Simpson compuesta para calcular la integral definida de a a b de f(c,x)
// Entrada: a,b = extremos de integración; n = cantidad de subintervalos; f = función de Scilab; c = real a aplicar como primer parametro de f
// Salida: w = aproximación de la integral de f en [a,b] por la Regla de Simpson
function w = SimpsonFD(a, b, n, f, c)
    h = (b-a)/n;
    w = (f(c,a)+f(c,b));
    for i=1:n-1
        w = w+f(c,a+i*h)*(2*modulo(i,2)+2);
    end
    w = w * (h/3);
endfunction

//Funcion que aplica regla de Simpson compuesta para calcular la integral definida de a a b y de c(x) a d(x) de f
// Entrada: a,b = extremos de integración exteriores;
// c,d = funciones para los extremos de integracion interiores;
// n = cantidad de subintervalos para la integral externa; 
//m = cantidad de intervalos para la integral interna;
//f = función de Scilab
// Salida: w = aproximación de la integral de la integral de f en [c(y),d(y)], en [a,b] por la Regla de Simpson
function w = SimpsonDoble(a,b,c,d,n,m,f)
    h = (b-a)/n;
    w = (SimpsonFD(c(a),d(a),m,f,a)+SimpsonFD(c(b),d(b),m,f,b));
    for i=1:n-1
        xi = a+i*h;
        w = w+SimpsonFD(c(xi),d(xi),m,f, xi)*(2*modulo(i,2)+2);
    end
    w = w * (h/3);
endfunction

function z = uno(x,y)
    z=1;
endfunction
function y = cx1(x)
    y = -sqrt(2*x-x^2);
endfunction
function y = dx1(x)
    y = sqrt(2*x-x^2);
endfunction
T6 = DobleTrap(0,2, cx1, dx1, 100, 100, uno);
S6 = SimpsonDoble(0,2, cx1, dx1, 100, 100, uno);
