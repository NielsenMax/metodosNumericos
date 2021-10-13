clear
clc
//ejercicio 1
//si b<0: x- = (14) y x+ = (7)
//si b>0: x- = (6) y x+ = (15)
// p es un polinomio 
function r = raices(p)
    c = coeff(p, 0);
    b = coeff(p, 1);
    a = coeff(p, 2);
    disc = (b^2-4*a*c);
    if disc > 0 then
        if b < 0 then
            r(1) = 2*c/(-b+sqrt(disc));
            r(2) = (-b+sqrt(disc))/(2*a);
        else
            r(1) = (-b-sqrt(disc))/(2*a);
            r(2) = 2*c/(-b-sqrt(disc));
        end
    elseif disc == 0 then
        r(1) = -(b/(2*a));
        r(2) = r(1);
    else
        r(1) = complex(-b/(2*a), sqrt(-disc)/(2*a));
        r(2) = complex(-b/(2*a), -sqrt(-disc)/(2*a));
    end
endfunction

function b = Horner(a, x_0)
    n = length(a);
    b = a(n);
    for i = 1:(n-1) do
        b = a(n-i) + x_0*b
    end
endfunction

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

//toma funcion(string), valor, orden, h
function r = derivar(f, v, n, h)
    deff("y=DF0(x)","y="+f);
    if n == 0 then r = DF0(v);
    else
        for i = 1:(n-1) do
            deff("y=DF"+string(i)+"(x)", "y=(DF"+string(i-1)+"(x+"+string(h)+")-DF"+string(i-1)+"(x))/"+string(h));
        end
        deff("y=DFn(x)", "y=(DF"+string(n-1)+"(x+"+string(h)+")-DF"+string(n-1)+"(x))/"+string(h));
        r = DFn(v); 
    end
endfunction
//toma funcion(string), valor, orden, h
function r = derivarNum(f, v, n, h)
    deff("y=DF0(x)", "y="+f);
    if n == 0 then r = DF0(v);
    else
        for i = 1:(n-1) do
            deff("y=DF"+string(i)+"(x)", "y=numderivative(DF"+string(i-1)+",x,"+string(h)+",4)");
        end
        deff("y=DFn(x)", "y=numderivative(DF"+string(n-1)+",x,"+string(h)+",4)");
        r = DFn(v);
    end
endfunction

function r = taylor(f, v, n, a, h)
    coef = [0:n];
    coef = 1./(factorial(coef));
    deff("y=DF0(x)", "y="+f);
    for i = 0:n do
        if i == 0 then deff("y=F(x)","y=DF0(x)");
        else
            deff("y=DF"+string(i)+"(x)", "y=numderivative(DF"+string(i-1)+",x,"+string(h)+",4)");
            deff("y=F(x)","y=DF"+string(i)+"(x)");
        end
        coef(i+1) = coef(i+1)*F(a);printf("Esperado : %e\n", e1);
printf("misraices (nuestro) : %e (error= %e)\n", r2, error1);
    end
    p = poly(coef, "x", "coeff");
    r = horner(p, v-a);
endfunction

function r = TaylorFun(f, v, n, a, h)
    deff("y=DF0(x)", "y="+f);
    deff("y=Tf0(x)", "y=DF0(a)");
    if n == 0 then r = Tf0(v);
    else
        for i = 1:(n-1) do
            deff("y=DF"+string(i)+"(x)", "y=numderivative(DF"+string(i-1)+",x,"+string(h)+",4)");
            deff("y=Tf"+string(i)+"(x)", "y=Tf"+string(i-1)+"(x)+DF"+string(i)+"(a)*(x-a)^("+string(i)+")/factorial("+string(i)+")");
        end
        deff("y=DFn(x)", "y=numderivative(DF"+string(n-1)+",x,"+string(h)+",4)");
        deff("y=Tfn(x)", "y=Tf"+string(n-1)+"(x)+DFn(a)*(x-a)^(n)/factorial(n)");
    end
endfunction

//Ejercicio 1
p = poly([-0.0001 10000.0 0.0001],"x","coeff");
e1 = 1e-8;
roots1 = raices(p);
r1 = roots1(1);
roots2 = roots(p);
r2 = roots2(2);
error1 = abs(r2-e1)/e1;
error2 = abs(r2-e1)/e1;

printf("roots (Scilab) : %e (error= %e)\n", r2, error2);

//Esperado : 1.000000e-08
//misraices (nuestro) : 1.000000e-08 (error= 0.000000e+00)
//roots (Scilab) : 1.000000e-08 (error= 0.000000e+00)


//ejercicio 3
a = [1 , -2 ,3 ,1];
x_0 = 2;
p_x0 = Horner (a , x_0 );
p = poly (a ,"x", "coeff");
p_x0_scilab = horner (p , x_0 ) ;
disp (" p_x0 = "+ string ( p_x0 ))
disp (" p_x0_scilab = "+ string ( p_x0_scilab ))
//parte d
[p_x0, dp_x0] = gHorner(a, x_0);
disp("p_x0 = " + string(p_x0));
disp("dp_x0 = " + string(dp_x0));
//p_x0 = 17
//p_x0_scilab = 17
//p_x0 = 17
//dp_x0 = 22


//ejercicio 4
f = "%e^x";
v = 2;
h = 1e-2;
n = 3;
r0 = derivar(f,v,n,h);
r1 = derivarNum(f,v,n,h);
v_e = 6;
error0 = abs(r0-v_e)/v_e;
error1 = abs(r1-v_e)/v_e;
disp("r0="+string(r0)+" error="+string(error0));
disp("r1="+string(r1)+" error="+string(error1));
// r0=7.5008211 error=0.2501369
//r1=7.3890561 error=0.2315093

//ejercicio 5
funcprot(0);
r = taylor("%e^x", -2, 4, 0, 1e-2);
disp(string(r));
disp(string((r-(%e^-2))/(%e^-2)))
//0.3333333
//1.4630187

