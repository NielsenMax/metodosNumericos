function x = secante(fun, x0, x1, tol, iter)
    /*
    Toma una funcion en forma de string y 4 reales.
    Devuelve una raiz de la funcion  fun, utilizando el metodo de la secante partiendo de los puntos x0 y x1.
    El resultado sera dentro de la tolerancia tol y realizando menos iteraciones que iter
    */
    deff("y=f(x)","y="+fun);
    i = 0;
    x2 = x1 - f(x1)*((x1-x0)/(f(x1)-f(x0)));
    while abs(x2-x1)>tol && i < iter
        i = i+1;
        x0 = x1;
        x1 = x2;
        x2 = x1 - f(x1)*((x1-x0)/(f(x1)-f(x0)));
    end
    if (abs(x2-x1)> tol) then disp('Se alcanzo el máximo de iteraciones'); end
    x = x2;
endfunction
fun = "(x^2/4) - sin(x)";
deff("y=f(x)", "y="+fun);
x0 = 1;
x1 = 2;
tol = 10^-6;
iter = 100;
x = secante(fun, x0, x1, tol, iter);
disp("El punto es x="+string(x)+". EL valor de es f(x)="+string(f(x)));
