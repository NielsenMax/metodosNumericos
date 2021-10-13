function x = biseccion(fun, a, b, tol)
    /*
    La funcion toma un string que sera una funcion y tres reales.
    Luego devolvera un valor de x, entre a y b, donde la funcion fun aplicada en este sera 0 dentro de la tolerancia tol pasada.
    sera 
    */
    deff("y=f(x)","y="+fun);
    if f(a)*f(b) >= 0 then
        error('biseccion - El signo en los extremos debe ser opuesto');
        abort;
    end
    c = (a+b)/2;
    while abs(b-c) > tol
        if f(c)*f(a)<0 then
            b = c;
        else
            a = c;
        end
        c = (a+b)/2;
    end
    x=c;
endfunction

fun = "x^3+x^2-x-1";
deff("y=f(x)","y="+fun);
a = 0;
b = 1.5;
tol = 10^-6;
x = biseccion(fun,a,b,tol);
disp("El punto es x="+string(x)+". EL valor de es f(x)="+string(f(x)));
